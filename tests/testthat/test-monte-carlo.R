make_small_data <- function(seed = 42, n = 60, binary_td = FALSE) {
  set.seed(seed)
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 6)) + 0.01, n, 6)
  X <- matrix(rnorm(n * 3), n, 3)
  if (binary_td) X[, 1:2] <- rbinom(n * 2, 1, 0.5)
  prepare_gbkmr_data(Y, Z, X,
    time_points = 3, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    td_covariate_names = "waist", log_transform_mixtures = TRUE)
}

# n = 30 with 10 knots makes fields::cover.design() warn about the number of
# nearest neighbours; that warning is not the subject of these tests.
run_small <- function(dat, n_cores = 1, n_knots = 10, K = 4) {
  suppressWarnings(gbkmr_run(
    data = dat, outcome = "Y", outcome_type = "continuous",
    time_points = 3, iter = 200, n = 30, K = K,
    n_knots = n_knots, n_cores = n_cores, engine = "bkmr", verbose = FALSE))
}

test_that("continuous confounder draws vary across Monte Carlo samples", {
  skip_if_not_installed("bkmr")
  res <- run_small(make_small_data())
  L1a <- res$raw_results$L_samp_a[[1]][[1]]   # draws x K at t = 1
  expect_equal(dim(L1a), c(length(res$raw_results$meta$sel), 4))
  # Before the fix every row was constant (same seed, same kernel inputs).
  expect_true(all(apply(L1a, 1, function(r) length(unique(r)) > 1)))
  expect_equal(res$raw_results$meta$confounder_noise, "residual")
})

test_that("binary confounder draws are Bernoulli and vary across samples", {
  skip_if_not_installed("bkmr")
  res <- run_small(make_small_data(binary_td = TRUE), K = 40)
  L1a <- res$raw_results$L_samp_a[[1]][[1]]
  expect_true(all(L1a %in% c(0, 1)))
  expect_equal(res$call_info$confounder_types[["waist"]], "binary")
  expect_true(any(apply(L1a, 1, function(r) length(unique(r)) > 1)))
})

test_that("results do not depend on the number of workers", {
  skip_if_not_installed("bkmr")
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2)
  dat <- make_small_data()
  r1 <- run_small(dat, n_cores = 1)
  r2 <- run_small(dat, n_cores = 2)
  expect_equal(r1$causal_effect$estimate, r2$causal_effect$estimate)
  expect_equal(r1$raw_results$Ya_mat, r2$raw_results$Ya_mat)
  expect_equal(r1$raw_results$L_samp_a, r2$raw_results$L_samp_a)
  expect_equal(r1$raw_results$fit_y$beta, r2$raw_results$fit_y$beta)
  expect_equal(r2$raw_results$meta$mc_cores, 2L)
})

test_that("n_knots = NULL fits the exact Gaussian process", {
  skip_if_not_installed("bkmr")
  dat <- make_small_data()
  res <- run_small(dat, n_knots = NULL)
  expect_null(res$raw_results$meta$n_knots)
  # bkmr stores the predictive-process knots in fit$data.comps$knots.
  expect_null(res$raw_results$fit_y$data.comps$knots)
  res10 <- run_small(dat, n_knots = 10)   # positive control
  expect_equal(res10$raw_results$meta$n_knots, 10)
  expect_equal(nrow(res10$raw_results$fit_y$data.comps$knots), 10)
  res0 <- run_small(dat, n_knots = 0)
  expect_null(res0$raw_results$meta$n_knots)
})

test_that("Monte Carlo seeds are distinct across stages, draws, times and confounders", {
  grid <- expand.grid(stage = 1:4, t = 0:5, li = 0:5, j = 1:600)
  seeds <- mapply(.gbkmr_seed, 7L, grid$stage, grid$t, grid$li, grid$j)
  expect_equal(anyDuplicated(seeds), 0L)
  expect_true(all(is.finite(seeds)))
})

test_that(".gbkmr_lapply_mc reports killed workers instead of dropping results", {
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2)
  expect_equal(.gbkmr_lapply_mc(1:3, function(i) i * 2, 1L), list(2, 4, 6))
  expect_equal(.gbkmr_lapply_mc(1:3, function(i) i * 2, 2L), list(2, 4, 6))
  expect_error(.gbkmr_lapply_mc(1:2, function(i) stop("boom"), 2L), "boom")
  expect_error(
    suppressWarnings(.gbkmr_lapply_mc(1:2, function(i) if (i == 2) quit(save = "no") else i, 2L)),
    "returned no result"
  )
})

test_that("causal plots are reproducible for a given seed", {
  skip_if_not_installed("bkmr")
  res <- run_small(make_small_data())
  sel2 <- tail(res$raw_results$meta$sel, 2)
  p1 <- gbkmr_causal_overall(res, quantiles = c(0.25, 0.75), K = 3, sel = sel2, seed = 8)
  p2 <- gbkmr_causal_overall(res, quantiles = c(0.25, 0.75), K = 3, sel = sel2, seed = 8)
  expect_equal(p1$draws, p2$draws)
  p3 <- gbkmr_causal_overall(res, quantiles = c(0.25, 0.75), K = 3, sel = sel2, seed = 9)
  expect_false(isTRUE(all.equal(p1$draws, p3$draws)))
  # Serial and forked evaluation agree.
  withr::local_options(causalBKMR.cores = 1)
  p4 <- gbkmr_causal_overall(res, quantiles = c(0.25, 0.75), K = 3, sel = sel2, seed = 8)
  expect_equal(p1$draws, p4$draws)
})

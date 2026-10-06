# Core g-BKMR analysis implementation

# Number of forked workers to use for the standard-BKMR path. Forking is not
# available on Windows, so the work runs serially there.
.gbkmr_mc_cores <- function(n_cores) {
  if (is.null(n_cores) || !is.numeric(n_cores) || length(n_cores) != 1L ||
      !is.finite(n_cores) || n_cores < 1) {
    return(1L)
  }
  if (.Platform$OS.type == "windows") return(1L)
  available <- suppressWarnings(parallel::detectCores())
  if (is.na(available) || available < 1) available <- 1L
  as.integer(max(1L, min(n_cores, available)))
}

# lapply() over posterior draws / model fits, forked when mc_cores > 1. Every
# FUN call seeds its own RNG, so the result is the same for any mc_cores.
# mclapply() reports an R error as a "try-error" element and a killed worker
# (e.g. by the out-of-memory killer) as NULL plus a warning; both are turned
# into errors here because downstream code would otherwise silently drop the
# missing draws.
.gbkmr_lapply_mc <- function(X, FUN, mc_cores) {
  if (mc_cores <= 1L) return(lapply(X, FUN))
  out <- parallel::mclapply(X, FUN, mc.cores = mc_cores, mc.preschedule = TRUE)
  failed <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    e <- out[[which(failed)[1]]]
    cond <- attr(e, "condition")
    stop("A parallel worker failed: ",
         if (is.null(cond)) as.character(e) else conditionMessage(cond),
         call. = FALSE)
  }
  missing <- vapply(out, is.null, logical(1))
  if (length(out) != length(X) || any(missing)) {
    stop(sum(missing), " of ", length(X), " parallel tasks returned no result ",
         "(worker killed, probably out of memory). Re-run with a smaller ",
         "n_cores or more memory.", call. = FALSE)
  }
  out
}

# Seeds for the Monte Carlo g-computation. Stages: 1 = posterior function draw
# for a confounder model, 2 = posterior function draw for the outcome model,
# 3 = residual noise / Bernoulli uniforms for a confounder, 4 = BKMR model fit.
# The multipliers keep every (stage, t, li, j) combination distinct for
# t, li <= 5 and j <= 9000, and `base` (currind, or the seed argument of the
# plotting functions) separates analyses.
.gbkmr_seed <- function(seed, stage, t = 0L, variable = 0L, draw = 0L, mc = 0L) {
  value <- seed + stage * 1000003 + t * 10007 + variable * 1009 +
    draw * 101 + mc
  as.integer(value %% .Machine$integer.max)
}

#' Run g-BKMR panel analysis
#'
#' @param sim_popn Data frame in g-BKMR format (see \code{\link{prepare_gbkmr_data}}).
#' @param T Integer. Number of time points (including t=0).
#' @param p Integer. Number of exposures per time point.
#' @param confounder_basenames Character vector. Base names for time-dependent confounders.
#' @param confounder_types Character vector. "continuous" or "binary" for each
#'   time-dependent confounder. Binary time-varying confounders use probit BKMR and Bernoulli
#'   Monte Carlo sampling under engine="bkmr".
#' @param common_covariates Character vector. Baseline covariate names.
#' @param currind Integer. Random seed.
#' @param n Integer or NULL. Sample size for analysis. If NULL, all rows are used.
#' @param K Integer. Monte Carlo samples for g-computation.
#' @param sel Numeric vector. Post-burn-in MCMC indices for inference.
#' @param iter Integer. MCMC iterations for time-varying confounder models.
#' @param n_iter Integer or NULL. MCMC iterations for outcome model (default: iter).
#' @param n_knots Integer or NULL. Number of knots for the BKMR predictive-process
#'   approximation under engine="bkmr". NULL (or 0) fits the exact Gaussian
#'   process, as in Chai et al.; this is slower but gives wider, better
#'   calibrated posterior intervals than the knot approximation.
#' @param engine Character. Fitting engine: "bkmr" or "fastbkmr".
#' @param n_subset Integer. Number of subsets for fastBKMR.
#' @param n_cores Integer or NULL. Under engine="fastbkmr", cores for the
#'   subset fits (default 10). Under engine="bkmr", forked workers
#'   (\code{parallel::mclapply}) used to fit the BKMR models and to run the
#'   Monte Carlo g-computation across posterior draws (default 1 = serial;
#'   always 1 on Windows). Results do not depend on the number of workers.
#' @param outcome_type Character. "continuous" (Gaussian BKMR, default) or
#'   "binary" (probit BKMR via family="binomial"). Binary outcome requires
#'   engine="bkmr" (fastBKMR does not yet support non-Gaussian outcomes).
#' @param a_probs Numeric vector of length 2. Quantile probabilities for
#'   intervention levels (default: c(0.25, 0.75)).
#' @param a_vals Named numeric vector or NULL. Custom intervention values for
#'   low-exposure scenario. Overrides a_probs if provided.
#' @param astar_vals Named numeric vector or NULL. Custom intervention values for
#'   high-exposure scenario. Overrides a_probs if provided.
#' @param verbose_every Integer. Kept for backward compatibility; progress is
#'   now reported once per model and per sampling stage.
#'
#' @return A list with causal effect estimate and model fits.
#' @importFrom stats complete.cases quantile rnorm runif
#' @export
run_gbkmr_panel <- function(
    sim_popn,
    T = 5,
    p = 3,
    confounder_basenames = c("waist"),
    confounder_types = NULL,
    common_covariates = "baseline_1",
    currind = 1,
    n = NULL,
    K = 1000,
    sel = seq(22000, 24000, by = 25),
    iter = 24000,
    n_iter = NULL,
    n_knots = 50,
    engine = c("bkmr", "fastbkmr"),
    n_subset = 10,
    n_cores = NULL,
    outcome_type = c("continuous", "binary"),
    a_probs = c(0.25, 0.75),
    a_vals = NULL,
    astar_vals = NULL,
    verbose_every = 50) {

  engine <- match.arg(engine)
  outcome_type <- match.arg(outcome_type)

  if (is.null(n_iter)) n_iter <- iter
  if (is.null(n_cores)) n_cores <- if (engine == "fastbkmr") 10L else 1L
  use_knots <- !is.null(n_knots) && is.finite(n_knots) && n_knots > 0
  if (!use_knots) n_knots <- NULL
  # Forked workers for the standard-BKMR path (fits + Monte Carlo). Every
  # worker sets its own seeds, so results are identical for any mc_cores.
  mc_cores <- if (engine == "bkmr") .gbkmr_mc_cores(n_cores) else 1L
  if (!"Y" %in% names(sim_popn)) stop("Data must contain outcome variable 'Y'")
  if (!"id" %in% names(sim_popn)) stop("Data must contain 'id' column")
  if (is.null(n)) n <- nrow(sim_popn)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1) {
    stop("n must be a positive integer or NULL.")
  }
  n <- as.integer(n)
  if (n > nrow(sim_popn)) {
    warning("n larger than data; using all rows.")
    n <- nrow(sim_popn)
  }
  if (max(sel) > max(iter, n_iter)) stop("sel contains indices beyond total MCMC iterations!")

  if (is.null(confounder_types)) {
    confounder_types <- rep("continuous", length(confounder_basenames))
  }
  if (length(confounder_types) != length(confounder_basenames)) {
    stop("confounder_types must have the same length as confounder_basenames")
  }
  if (length(confounder_basenames) == 0L) {
    confounder_types <- character(0)
  } else {
    confounder_types <- match.arg(confounder_types, c("continuous", "binary"),
                                several.ok = TRUE)
    names(confounder_types) <- confounder_basenames
  }

  # fastBKMR's public skmbayes() path is Gaussian-only in the current
  # fbkmr package, so binary outcomes or binary time-varying confounders must use standard BKMR.
  if (engine == "fastbkmr" &&
      (outcome_type == "binary" || any(confounder_types == "binary"))) {
    stop("Binary outcomes or binary time-varying confounders are not supported ",
         "with engine='fastbkmr'.\n",
         "  fbkmr::skmbayes() uses the Gaussian fast path in the current package.\n",
         "  Use engine='bkmr' for probit BKMR.")
  }

  # --- Internal helpers ---
  .fit_model <- function(y, Z_sc, X, it, knots = NULL, family = "gaussian") {
    if (engine == "fastbkmr") {
      if (!requireNamespace("fbkmr", quietly = TRUE))
        stop("Package 'fbkmr' is required for engine='fastbkmr'.\n",
             "Install with: remotes::install_github('junwei-lu/fbkmr')")
      nc <- min(n_subset, parallel::detectCores() - 1, n_cores)
      use_parallel <- nc > 1
      tryCatch(
        fbkmr::skmbayes(Z = Z_sc, X = X, y = y,
                         n_subset = n_subset, n_samp = 200,
                         iter = it, varsel = TRUE, est.h = FALSE,
                         parallel = use_parallel, n_cores = nc),
        error = function(e) {
          if (use_parallel && grepl("doSnowGlobals|parallel|snow", e$message, ignore.case = TRUE)) {
            warning("Parallel execution failed, falling back to sequential mode.\n",
                    "  Original error: ", e$message, call. = FALSE)
            fbkmr::skmbayes(Z = Z_sc, X = X, y = y,
                             n_subset = n_subset, n_samp = 200,
                             iter = it, varsel = TRUE, est.h = FALSE,
                             parallel = FALSE, n_cores = 1)
          } else {
            stop(e)
          }
        }
      )
    } else {
      bkmr::kmbayes(y = y, Z = Z_sc, X = X, iter = it,
                     family = family,
                     varsel = TRUE, verbose = FALSE, knots = knots)
    }
  }

  # Posterior-predictive draw of h(z) + x'beta at MCMC iteration sel_j for the
  # paired rows (Znew_a[k, ], Znew_astar[k, ]), k = 1..K, with common random
  # numbers across k (see .gbkmr_predict_blocks). Also returns the residual SD
  # of the same posterior draw from the same fit.
  .predict_pairs <- function(fit, Znew_a, Znew_astar, Xnew, sel_j, seed,
                             type = "link") {
    blocks <- lapply(seq_len(nrow(Znew_a)), function(k) {
      rbind(Znew_a[k, ], Znew_astar[k, ])
    })
    res <- .gbkmr_predict_blocks(fit, blocks, Xnew, sel_j, seed, type)
    list(a = res$pred[, 1], astar = res$pred[, 2], sigma = res$sigma)
  }

  .lapply_mc <- function(X, FUN) .gbkmr_lapply_mc(X, FUN, mc_cores)

  .extract_beta <- function(fit, sel_idx) {
    if (is.list(fit) && !inherits(fit, "bkmrfit")) {
      betas <- lapply(fit, function(f) colMeans(f$beta[sel_idx, , drop = FALSE]))
      colMeans(do.call(rbind, betas))
    } else {
      colMeans(fit$beta[sel_idx, , drop = FALSE])
    }
  }

  # --- Sampling and naming ---
  set.seed(currind)
  dat_sim <- sim_popn[sample(seq_len(nrow(sim_popn)), n, replace = FALSE), ]

  exposure_times <- 0:(T - 1)
  confounder_times <- if (T > 1) seq_len(T - 1L) else integer(0)
  exposure_names_at_t <- function(t) paste0("logM", seq_len(p), "_", t)
  all_exposure_names <- unlist(lapply(exposure_times, exposure_names_at_t),
                               use.names = FALSE)
  confounder_names_at_t <- function(t) paste0(confounder_basenames, "_", t)
  all_confounder_names <- unlist(lapply(confounder_times, confounder_names_at_t),
                               use.names = FALSE)

  needed_cols <- c("Y", "id", common_covariates, all_exposure_names, all_confounder_names)
  miss <- setdiff(needed_cols, names(dat_sim))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  X_common <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(common_covariates)))
  X_predict_common <- matrix(colMeans(X_common), nrow = 1)
  fitkm_list <- vector("list", length(confounder_times))
  scaleinfo_list <- vector("list", length(confounder_times))
  if (length(confounder_times) > 0L) {
    names(fitkm_list) <- paste0("L", confounder_times)
    names(scaleinfo_list) <- names(fitkm_list)
  }

  # --- Knot helper (only used for engine == "bkmr" with knots requested) ---
  .compute_knots <- function(Z_sc, n_knots) {
    n_unique <- nrow(unique(round(Z_sc, 10)))
    if (n_unique < nrow(Z_sc)) {
      Z_sc <- Z_sc + matrix(rnorm(length(Z_sc), 0, 1e-6), nrow = nrow(Z_sc))
      n_unique <- nrow(unique(round(Z_sc, 10)))
    }
    nd <- min(n_knots, n_unique - 1, floor(n_unique * 0.9))
    if (nd < 2) stop("Not enough unique rows to place knots")
    tryCatch(
      fields::cover.design(Z_sc, nd = nd)$design,
      error = function(e) Z_sc[sample(nrow(Z_sc), nd), , drop = FALSE]
    )
  }

  # =========================================================================
  # 1) Fit time-varying confounder models
  # =========================================================================
  # The confounder models and the outcome model are fit independently of one
  # another, so they are collected as jobs and fit together (in parallel when
  # mc_cores > 1). Each job seeds its own RNG from currind so the fits do not
  # depend on the number of workers or on the order of execution.
  fit_jobs <- list()
  .add_job <- function(name, y, Z_sc, X, it, knots, family) {
    fit_jobs[[length(fit_jobs) + 1L]] <<- list(
      name = name, y = y, Z_sc = Z_sc, X = X, it = it, knots = knots,
      family = family,
      seed = .gbkmr_seed(currind, 4L, draw = length(fit_jobs) + 1L)
    )
  }
  message("Preparing time-varying confounder models ...")

  for (t in confounder_times) {
    y_cols <- confounder_names_at_t(t)
    y_mat  <- as.matrix(dat_sim[, y_cols, drop = FALSE])
    colnames(y_mat) <- y_cols

    # Z: exposures 0..t-1 + confounders 1..t-1
    Z_names <- unlist(lapply(0:(t - 1), exposure_names_at_t))
    if (t > 1) Z_names <- c(Z_names, unlist(lapply(1:(t - 1), confounder_names_at_t)))

    Z_raw <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(Z_names)))
    rows_ok_ZX <- complete.cases(Z_raw, X_common)
    if (sum(rows_ok_ZX) < 3) stop("Not enough complete rows for time-varying confounder at t=", t)

    Z_sc <- scale(Z_raw[rows_ok_ZX, , drop = FALSE])
    sc_center <- attr(Z_sc, "scaled:center")
    sc_scale  <- attr(Z_sc, "scaled:scale")
    scaleinfo_list[[t]] <- list(center = sc_center, scale = sc_scale)

    knots_t <- if (engine == "bkmr" && use_knots) .compute_knots(Z_sc, n_knots) else NULL

    for (li in seq_len(ncol(y_mat))) {
      y_vec <- y_mat[, li]
      y_ok  <- y_vec[rows_ok_ZX]
      mask_y <- !is.na(y_ok)
      if (sum(mask_y) < 3) {
        stop(sprintf("Not enough complete rows for time-varying confounder %s at t=%d",
                     colnames(y_mat)[li], t))
      }
      valid_idx <- which(rows_ok_ZX)[mask_y]
      Z_sc_fit <- scale(Z_raw[valid_idx, , drop = FALSE],
                        center = sc_center, scale = sc_scale)
      X_common_fit <- X_common[valid_idx, , drop = FALSE]
      y_vec_fit <- y_vec[valid_idx]

      confounder_family <- if (confounder_types[[li]] == "binary") "binomial" else "gaussian"
      message(sprintf("  L%d: %s [engine=%s, Z=%d cols, n=%d, family=%s]",
                      t, colnames(y_mat)[li], engine, ncol(Z_sc_fit),
                      length(y_vec_fit), confounder_family))

      .add_job(paste0("L", t, "_", li), y_vec_fit, Z_sc_fit, X_common_fit,
               iter, knots_t, confounder_family)
    }
  }

  # =========================================================================
  # 2) Fit outcome model
  # =========================================================================
  message("Preparing outcome model Y ...")
  Y <- dat_sim$Y
  valid_y <- !is.na(Y)
  if (any(!valid_y)) {
    message(sprintf("  Removing %d NA outcomes", sum(!valid_y)))
    dat_sim <- dat_sim[valid_y, ]
    Y <- Y[valid_y]
    X_common <- X_common[valid_y, , drop = FALSE]
    X_predict_common <- matrix(colMeans(X_common), nrow = 1)
  }

  Zy_names <- c(all_exposure_names, all_confounder_names)
  Zy_raw   <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(Zy_names)))
  Zy_sc    <- scale(Zy_raw)
  scale_info_y <- list(center = attr(Zy_sc, "scaled:center"),
                       scale  = attr(Zy_sc, "scaled:scale"))

  knots_y <- if (engine == "bkmr" && use_knots) .compute_knots(Zy_sc, n_knots) else NULL

  y_family <- if (outcome_type == "binary") "binomial" else "gaussian"
  message(sprintf("  Y [engine=%s, Z=%d cols, n=%d, family=%s]",
                  engine, ncol(Zy_sc), length(Y), y_family))
  .add_job("Y", Y, Zy_sc, X_common, n_iter, knots_y, y_family)

  message(sprintf("Fitting %d BKMR model(s) [knots=%s, workers=%d] ...",
                  length(fit_jobs), if (use_knots) n_knots else "none", mc_cores))
  start_time_fit <- proc.time()
  fits <- .lapply_mc(fit_jobs, function(job) {
    set.seed(job$seed)
    .fit_model(job$y, job$Z_sc, job$X, job$it, job$knots, family = job$family)
  })
  names(fits) <- vapply(fit_jobs, `[[`, character(1), "name")
  message(sprintf("  fits done in %.1f min",
                  (proc.time() - start_time_fit)["elapsed"] / 60))

  for (t in confounder_times) {
    fitkm_list[[t]] <- lapply(seq_along(confounder_basenames), function(li) {
      fits[[paste0("L", t, "_", li)]]
    })
  }
  fit_y <- fits[["Y"]]

  # =========================================================================
  # 3) Intervention levels (a / a*)
  # =========================================================================
  A_all <- as.matrix(dplyr::select(dat_sim, dplyr::all_of(all_exposure_names)))
  if (!is.null(a_vals) && !is.null(astar_vals)) {
    a_vec     <- a_vals[all_exposure_names]
    astar_vec <- astar_vals[all_exposure_names]
  } else {
    a_vec     <- apply(A_all, 2, quantile, probs = a_probs[1])
    astar_vec <- apply(A_all, 2, quantile, probs = a_probs[2])
  }

  # Containers
  L_samp_a     <- vector("list", length(confounder_times))
  L_samp_astar <- vector("list", length(confounder_times))

  scale_like <- function(newZ, center, sc) scale(newZ, center = center, scale = sc)

  # =========================================================================
  # 4) Sequential time-varying confounder sampling
  # =========================================================================
  # For posterior draw j and Monte Carlo sample k, a continuous confounder is
  # drawn as  L = h_j(z_k) + x'beta_j + sigma_j * eps_jk,  eps_jk ~ N(0, 1),
  # i.e. from the fitted conditional distribution, not at its mean (Chai et
  # al.); a binary confounder is Bernoulli(Phi(h_j(z_k) + x'beta_j)). The
  # posterior-function draws use common random numbers across k (see
  # .gbkmr_predict_blocks); eps_jk and the Bernoulli uniforms are seeded per
  # (analysis, draw, time, confounder) so the K values are distinct and the
  # result is independent of mc_cores.
  message("\n=== Sampling time-varying confounders sequentially ===")
  start_time_global <- proc.time()

  for (t in confounder_times) {
    message(sprintf("\n--- Time point t=%d ---", t))
    start_time_t <- proc.time()

    a_exp_t     <- unlist(lapply(0:(t - 1), function(s) a_vec[exposure_names_at_t(s)]))
    astar_exp_t <- unlist(lapply(0:(t - 1), function(s) astar_vec[exposure_names_at_t(s)]))
    Za_exp_mat     <- matrix(a_exp_t,     nrow = K, ncol = length(a_exp_t),     byrow = TRUE)
    Zastar_exp_mat <- matrix(astar_exp_t, nrow = K, ncol = length(astar_exp_t), byrow = TRUE)

    L_samp_a_t     <- vector("list", length(confounder_basenames))
    L_samp_astar_t <- vector("list", length(confounder_basenames))

    for (li in seq_along(confounder_basenames)) {
      message(sprintf("  Sampling %s at t=%d [%d draws x K=%d, workers=%d]",
                      confounder_basenames[li], t, length(sel), K, mc_cores))

      fit_li    <- fitkm_list[[t]][[li]]
      scinfo_t  <- scaleinfo_list[[t]]
      is_binary <- confounder_types[[li]] == "binary"

      draws_j <- .lapply_mc(seq_along(sel), function(j) {
        # Historical confounder block for this posterior draw
        if (t == 1) {
          aL_a_j         <- Za_exp_mat
          astarL_astar_j <- Zastar_exp_mat
        } else {
          L_hist_a_blocks <- L_hist_astar_blocks <- list()
          for (tt in seq_len(t - 1L)) {
            for (lj in seq_along(confounder_basenames)) {
              L_hist_a_blocks[[length(L_hist_a_blocks) + 1]] <- L_samp_a[[tt]][[lj]][j, ]
              L_hist_astar_blocks[[length(L_hist_astar_blocks) + 1]] <- L_samp_astar[[tt]][[lj]][j, ]
            }
          }
          aL_a_j         <- cbind(Za_exp_mat,     do.call(cbind, L_hist_a_blocks))
          astarL_astar_j <- cbind(Zastar_exp_mat, do.call(cbind, L_hist_astar_blocks))
        }

        pred <- .predict_pairs(
          fit_li,
          Znew_a     = scale_like(aL_a_j,         scinfo_t$center, scinfo_t$scale),
          Znew_astar = scale_like(astarL_astar_j, scinfo_t$center, scinfo_t$scale),
          Xnew = X_predict_common, sel_j = sel[j],
          seed = .gbkmr_seed(currind, 1L, t, li, j),
          type = if (is_binary) "response" else "link"
        )

        set.seed(.gbkmr_seed(currind, 3L, t, li, j))
        if (is_binary) {
          u_a     <- stats::runif(K)
          u_astar <- stats::runif(K)
          list(a     = as.numeric(u_a     < pmin(pmax(pred$a,     0), 1)),
               astar = as.numeric(u_astar < pmin(pmax(pred$astar, 0), 1)))
        } else {
          eps_a     <- stats::rnorm(K)
          eps_astar <- stats::rnorm(K)
          list(a     = pred$a     + pred$sigma * eps_a,
               astar = pred$astar + pred$sigma * eps_astar)
        }
      })

      L_samp_a_t[[li]]     <- do.call(rbind, lapply(draws_j, `[[`, "a"))
      L_samp_astar_t[[li]] <- do.call(rbind, lapply(draws_j, `[[`, "astar"))

      message(sprintf("    done | %.2f min",
                      (proc.time() - start_time_t)["elapsed"] / 60))
    }

    L_samp_a[[t]]     <- L_samp_a_t
    L_samp_astar[[t]] <- L_samp_astar_t
  }

  # =========================================================================
  # 5) Sample outcome Y
  # =========================================================================
  # E[Y | a, L_k] under posterior draw j: h_j(z_k) + x'beta_j (continuous) or
  # Phi(h_j(z_k) + x'beta_j) (binary, probit BKMR). No residual noise is added
  # because the estimand is a mean.
  message(sprintf("\n=== Sampling outcome Y [%d draws x K=%d, workers=%d] ===",
                  length(sel), K, mc_cores))
  start_time_y <- proc.time()

  pT <- length(all_exposure_names)
  exp_a_block     <- matrix(a_vec[all_exposure_names],     nrow = K, ncol = pT, byrow = TRUE)
  exp_astar_block <- matrix(astar_vec[all_exposure_names], nrow = K, ncol = pT, byrow = TRUE)

  y_draws_j <- .lapply_mc(seq_along(sel), function(j) {
    if (length(confounder_times) > 0 && length(confounder_basenames) > 0) {
      L_a_blocks <- L_astar_blocks <- list()
      for (t in confounder_times) {
        for (li in seq_along(confounder_basenames)) {
          L_a_blocks[[length(L_a_blocks) + 1]]         <- L_samp_a[[t]][[li]][j, ]
          L_astar_blocks[[length(L_astar_blocks) + 1]] <- L_samp_astar[[t]][[li]][j, ]
        }
      }
      aL_a_j         <- cbind(exp_a_block,     do.call(cbind, L_a_blocks))
      astarL_astar_j <- cbind(exp_astar_block, do.call(cbind, L_astar_blocks))
    } else {
      aL_a_j         <- exp_a_block
      astarL_astar_j <- exp_astar_block
    }

    pred <- .predict_pairs(
      fit_y,
      Znew_a     = scale_like(aL_a_j,         scale_info_y$center, scale_info_y$scale),
      Znew_astar = scale_like(astarL_astar_j, scale_info_y$center, scale_info_y$scale),
      Xnew = X_predict_common, sel_j = sel[j],
      seed = .gbkmr_seed(currind, 2L, draw = j)
    )
    if (outcome_type == "binary") {
      list(a = stats::pnorm(pred$a), astar = stats::pnorm(pred$astar))
    } else {
      list(a = pred$a, astar = pred$astar)
    }
  })

  Ya_mat     <- do.call(rbind, lapply(y_draws_j, `[[`, "a"))
  Yastar_mat <- do.call(rbind, lapply(y_draws_j, `[[`, "astar"))
  stopifnot(nrow(Ya_mat) == length(sel), nrow(Yastar_mat) == length(sel))
  message(sprintf("  done | %.2f min", (proc.time() - start_time_y)["elapsed"] / 60))
  # =========================================================================
  # 6) Aggregate results
  # =========================================================================
  Ya     <- rowMeans(Ya_mat)
  Yastar <- rowMeans(Yastar_mat)
  diff_gBKMR <- mean(Yastar) - mean(Ya)

  beta_L <- lapply(seq_along(fitkm_list), function(t) {
    lapply(seq_along(fitkm_list[[t]]), function(li) {
      .extract_beta(fitkm_list[[t]][[li]], sel)
    })
  })
  beta_y <- .extract_beta(fit_y, sel)

  end_time_global <- proc.time()
  total_time_min <- round((end_time_global - start_time_global)["elapsed"] / 60, 2)
  message(sprintf("\n=== Total time: %s minutes ===", total_time_min))

  list(
    diff_gBKMR = diff_gBKMR,
    Ya = Ya,
    Yastar = Yastar,
    Ya_mat = Ya_mat,
    Yastar_mat = Yastar_mat,
    L_samp_a = L_samp_a,
    L_samp_astar = L_samp_astar,
    fit_confounders = fitkm_list,
    fit_y = fit_y,
    beta_L = beta_L,
    beta_y = beta_y,
    gcomp_state = list(
      exposure_data = A_all,
      confounder_scale_info = scaleinfo_list,
      outcome_scale_info = scale_info_y,
      baseline_predictors = X_predict_common
    ),
    meta = list(
      T = T, p = p,
      confounder_basenames = confounder_basenames,
      confounder_types = confounder_types,
      common_covariates = common_covariates,
      K = K, sel = sel, iter = iter, n_iter = n_iter,
      n_knots = n_knots, n = n,
      currind = currind,
      engine = engine, n_subset = n_subset, n_cores = n_cores,
      mc_cores = mc_cores,
      confounder_noise = "residual",   # L drawn from fitted conditional, not its mean
      outcome_type = outcome_type,
      a_probs = a_probs, a_vals = a_vals, astar_vals = astar_vals,
      a_vec = a_vec, astar_vec = astar_vec,
      exposure_names = all_exposure_names,
      confounder_names = all_confounder_names,
      total_time_minutes = total_time_min
    )
  )
}

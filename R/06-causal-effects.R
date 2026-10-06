# Causal effect summaries and plots

utils::globalVariables(".data")

.gbkmr_or <- function(x, y) {
  if (is.null(x)) y else x
}

.gbkmr_get_raw <- function(object) {
  if (inherits(object, "gbkmr_results")) {
    raw <- object$raw_results
  } else if (is.list(object) && all(c("fit_y", "meta") %in% names(object))) {
    raw <- object
  } else {
    stop("object must be a gbkmr_results object returned by gbkmr_run().")
  }

  if (is.null(raw$gcomp_state)) {
    stop(
      "This fitted object predates the causal plotting API and does not contain ",
      "the saved g-computation state. Refit the model with the current version ",
      "of causalBKMR."
    )
  }
  raw
}

.gbkmr_validate_probabilities <- function(x, name, min_length = 1L) {
  if (!is.numeric(x) || length(x) < min_length || any(!is.finite(x)) ||
      any(x < 0 | x > 1)) {
    stop(name, " must contain finite probabilities between 0 and 1.")
  }
  unique(as.numeric(x))
}

.gbkmr_validate_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("credible_level must be a single number between 0 and 1.")
  }
  as.numeric(level)
}

.gbkmr_select_exposures <- function(raw, exposures = NULL, time_points = NULL) {
  available <- raw$meta$exposure_names
  if (is.null(exposures)) {
    selected <- available
  } else {
    if (!is.character(exposures) || length(exposures) == 0L) {
      stop("exposures must be a non-empty character vector or NULL.")
    }
    missing_exposures <- setdiff(exposures, available)
    if (length(missing_exposures) > 0L) {
      stop("Unknown exposure columns: ", paste(missing_exposures, collapse = ", "))
    }
    selected <- available[available %in% unique(exposures)]
  }

  if (!is.null(time_points)) {
    if (!is.numeric(time_points) || any(!is.finite(time_points)) ||
        any(time_points < 0 | time_points != floor(time_points))) {
      stop("time_points must contain non-negative integer visit indices.")
    }
    exposure_times <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", selected)))
    selected <- selected[exposure_times %in% as.integer(time_points)]
  }

  if (length(selected) == 0L) {
    stop("No exposure columns remain after applying the requested filters.")
  }
  selected
}

.gbkmr_exposure_quantiles <- function(raw, probs) {
  exposure_names <- raw$meta$exposure_names
  exposure_data <- raw$gcomp_state$exposure_data
  if (is.null(colnames(exposure_data))) colnames(exposure_data) <- exposure_names
  exposure_data <- exposure_data[, exposure_names, drop = FALSE]

  values <- unlist(lapply(exposure_names, function(nm) {
    stats::quantile(exposure_data[, nm], probs = probs, na.rm = TRUE, names = FALSE)
  }), use.names = FALSE)
  out <- matrix(values, nrow = length(probs), ncol = length(exposure_names))
  colnames(out) <- exposure_names
  rownames(out) <- format(probs, trim = TRUE, scientific = FALSE)
  if (any(!is.finite(out))) stop("Exposure quantiles contain non-finite values.")
  out
}

.gbkmr_scale_matrix <- function(x, scale_info) {
  out <- sweep(x, 2L, scale_info$center, FUN = "-")
  sweep(out, 2L, scale_info$scale, FUN = "/")
}

# Posterior-predictive draw of h(z) + x'beta at MCMC iteration sel_j for each of
# K blocks of new rows (one block per Monte Carlo sample, n_regimes rows each).
# The seed is reset before every SamplePred() call, so the K blocks reuse the
# same standard-normal innovations (common random numbers, as in Chai et al.'s
# per-sample seed reset); this couples the posterior-uncertainty component
# across the Monte Carlo sample instead of averaging it away with K
# independent draws. Identical blocks are evaluated once. Returns a
# K x n_regimes matrix and the residual SD of the same posterior draw; when
# fastBKMR supplies a list of subset fits, one subset is chosen
# (deterministically from the seed) and used for both.
.gbkmr_predict_blocks <- function(fit, blocks, Xnew, sel_j, seed, type = "link") {
  set.seed(seed)
  if (is.list(fit) && !inherits(fit, "bkmrfit")) {
    fit <- fit[[sample.int(length(fit), 1L)]]
  }
  keys <- vapply(blocks, function(b) paste(b, collapse = "\r"), character(1))
  unique_pos <- which(!duplicated(keys))
  group <- match(keys, keys[unique_pos])
  out <- matrix(NA_real_, nrow = length(blocks), ncol = nrow(blocks[[1]]))
  for (u in seq_along(unique_pos)) {
    set.seed(seed)
    pr <- bkmr::SamplePred(fit, Znew = blocks[[unique_pos[u]]], Xnew = Xnew,
                           sel = sel_j, type = type)
    hit <- which(group == u)
    out[hit, ] <- matrix(as.numeric(pr), nrow = length(hit), ncol = ncol(out),
                         byrow = TRUE)
  }
  list(pred = out, sigma = sqrt(fit$sigsq.eps[sel_j]))
}

# .gbkmr_seed() is defined in 04-core-analysis.R and shared with gbkmr_run().

.gbkmr_validate_interventions <- function(raw, interventions) {
  if (is.data.frame(interventions)) interventions <- as.matrix(interventions)
  if (is.null(dim(interventions))) interventions <- matrix(interventions, nrow = 1L)
  if (!is.numeric(interventions)) stop("interventions must be numeric.")

  exposure_names <- raw$meta$exposure_names
  if (is.null(colnames(interventions))) {
    stop("interventions must have exposure column names.")
  }
  missing_columns <- setdiff(exposure_names, colnames(interventions))
  if (length(missing_columns) > 0L) {
    stop("Interventions are missing columns: ", paste(missing_columns, collapse = ", "))
  }
  interventions <- interventions[, exposure_names, drop = FALSE]
  storage.mode(interventions) <- "double"
  if (any(!is.finite(interventions))) {
    stop("interventions must contain only finite values.")
  }
  if (is.null(rownames(interventions))) {
    rownames(interventions) <- paste0("regime_", seq_len(nrow(interventions)))
  }
  interventions
}

.gbkmr_gcompute <- function(object, interventions, K = NULL, sel = NULL,
                             seed = 1L, verbose = FALSE) {
  raw <- .gbkmr_get_raw(object)
  interventions <- .gbkmr_validate_interventions(raw, interventions)
  meta <- raw$meta
  state <- raw$gcomp_state

  if (is.null(K)) K <- meta$K
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 1 ||
      K != floor(K)) {
    stop("K must be a positive integer or NULL.")
  }
  K <- as.integer(K)

  if (is.null(sel)) sel <- meta$sel
  if (!is.numeric(sel) || length(sel) == 0L || any(!is.finite(sel)) ||
      any(sel < 1 | sel != floor(sel))) {
    stop("sel must contain positive MCMC iteration indices.")
  }
  sel <- as.integer(sel)

  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) {
    stop("seed must be a non-negative integer.")
  }
  seed <- as.integer(seed)

  n_regimes <- nrow(interventions)
  n_draws <- length(sel)
  confounder_times <- if (meta$T > 1L) seq_len(meta$T - 1L) else integer(0)
  n_confounders <- length(meta$confounder_basenames)
  confounder_types <- meta$confounder_types
  if (is.null(confounder_types)) {
    confounder_types <- rep("continuous", n_confounders)
  }

  # Same Monte Carlo scheme as run_gbkmr_panel(): common random numbers for
  # the posterior-function draw across the K Monte Carlo samples of each
  # posterior draw, and continuous confounders drawn from the fitted
  # conditional distribution (mean + sigma_j * eps). Draws are seeded per
  # (seed, draw, time, confounder), so results do not depend on the number of
  # workers. Workers default to the fit's n_cores; override with
  # options(causalBKMR.cores = 1) e.g. when forking is unavailable.
  mc_cores <- .gbkmr_mc_cores(getOption("causalBKMR.cores", .gbkmr_or(meta$mc_cores, 1L)))
  lapply_mc <- function(X, FUN) .gbkmr_lapply_mc(X, FUN, mc_cores)

  L_samples <- vector("list", length(confounder_times))
  if (length(confounder_times) > 0L) names(L_samples) <- paste0("L", confounder_times)

  for (t in confounder_times) {
    L_samples_t <- vector("list", n_confounders)
    exposure_history <- interventions[, seq_len(t * meta$p), drop = FALSE]

    for (li in seq_len(n_confounders)) {
      fit_li <- raw$fit_confounders[[t]][[li]]
      scale_info <- state$confounder_scale_info[[t]]
      is_binary <- confounder_types[[li]] == "binary"

      per_draw <- lapply_mc(seq_len(n_draws), function(j) {
        blocks <- lapply(seq_len(K), function(k) {
          if (t > 1L) {
            history <- do.call(cbind, lapply(seq_len(t - 1L), function(tt) {
              do.call(cbind, lapply(seq_len(n_confounders), function(lj) {
                L_samples[[tt]][[lj]][j, k, ]
              }))
            }))
            newz <- cbind(exposure_history, history)
          } else {
            newz <- exposure_history
          }
          .gbkmr_scale_matrix(newz, scale_info)
        })

        pred <- .gbkmr_predict_blocks(
          fit_li, blocks, state$baseline_predictors, sel[j],
          seed = .gbkmr_seed(seed, 1L, t, li, j),
          type = if (is_binary) "response" else "link"
        )

        # The K residual draws (or Bernoulli uniforms) are shared by all
        # regimes (common random numbers): a regime's confounder draw then
        # depends only on its own kernel inputs, not on which other regimes
        # are evaluated in the same call, and contrasts between regimes carry
        # less Monte Carlo noise. gbkmr_run() instead uses independent noise
        # for a and a*, as in Chai et al.
        set.seed(.gbkmr_seed(seed, 3L, t, li, j))
        if (is_binary) {
          u <- stats::runif(K)
          (u < pmin(pmax(pred$pred, 0), 1)) * 1
        } else {
          pred$pred + pred$sigma * stats::rnorm(K)
        }
      })

      sampled <- array(NA_real_, dim = c(n_draws, K, n_regimes))
      for (j in seq_len(n_draws)) sampled[j, , ] <- per_draw[[j]]
      L_samples_t[[li]] <- sampled
    }
    L_samples[[t]] <- L_samples_t
    if (isTRUE(verbose)) message("Generated confounders for t = ", t)
  }

  per_draw_y <- lapply_mc(seq_len(n_draws), function(j) {
    blocks <- lapply(seq_len(K), function(k) {
      if (length(confounder_times) > 0L && n_confounders > 0L) {
        confounder_history <- do.call(cbind, lapply(confounder_times, function(t) {
          do.call(cbind, lapply(seq_len(n_confounders), function(li) {
            L_samples[[t]][[li]][j, k, ]
          }))
        }))
        newz <- cbind(interventions, confounder_history)
      } else {
        newz <- interventions
      }
      .gbkmr_scale_matrix(newz, state$outcome_scale_info)
    })
    pred <- .gbkmr_predict_blocks(
      raw$fit_y, blocks, state$baseline_predictors, sel[j],
      seed = .gbkmr_seed(seed, 2L, draw = j)
    )$pred
    if (meta$outcome_type == "binary") pred <- stats::pnorm(pred)
    pred
  })
  if (isTRUE(verbose)) message("Completed ", n_draws, " posterior draws")

  outcome_samples <- array(NA_real_, dim = c(n_draws, K, n_regimes))
  for (j in seq_len(n_draws)) outcome_samples[j, , ] <- per_draw_y[[j]]

  draw_means <- matrix(NA_real_, nrow = n_draws, ncol = n_regimes)
  for (r in seq_len(n_regimes)) {
    draw_means[, r] <- rowMeans(outcome_samples[, , r, drop = FALSE],
                                dims = 1L, na.rm = TRUE)
  }
  colnames(draw_means) <- rownames(interventions)

  list(
    draws = draw_means,
    iterations = sel,
    interventions = interventions,
    K = K
  )
}

.gbkmr_draw_interval <- function(draws, credible_level) {
  alpha <- (1 - credible_level) / 2
  c(
    estimate = mean(draws, na.rm = TRUE),
    lower = unname(stats::quantile(draws, alpha, na.rm = TRUE)),
    upper = unname(stats::quantile(draws, 1 - alpha, na.rm = TRUE))
  )
}

.gbkmr_new_causal_result <- function(type, summary, draws, interventions,
                                     plot, estimand, settings) {
  structure(
    list(
      summary = summary,
      draws = draws,
      interventions = interventions,
      plot = plot,
      estimand = estimand,
      settings = settings
    ),
    class = c(paste0("gbkmr_causal_", type), "gbkmr_causal_effects")
  )
}

.gbkmr_plot_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70")
    )
}

#' Overall mixture causal dose-response
#'
#' Estimates the full g-formula contrast between jointly setting every mixture
#' component to each requested marginal quantile and jointly setting every
#' component to the reference quantile. Time-dependent confounders are generated
#' sequentially under every intervention.
#'
#' @param object A fitted `gbkmr_results` object.
#' @param quantiles Joint mixture quantiles to evaluate.
#' @param reference Reference quantile for every mixture component.
#' @param K Number of Monte Carlo trajectories per posterior draw. `NULL` uses
#'   the value from the fitted model.
#' @param sel Posterior iteration indices. `NULL` uses the fitted model's
#'   post-burn-in indices.
#' @param seed Random seed for g-computation.
#' @param credible_level Posterior credible interval level.
#' @param verbose Logical; print g-computation progress.
#'
#' @return A `gbkmr_causal_effects` object containing tidy `summary` and
#'   `draws` data frames, the intervention matrix, estimand metadata, and a
#'   ggplot object in `plot`.
#' @export
gbkmr_causal_overall <- function(
    object,
    quantiles = seq(0.30, 0.90, by = 0.05),
    reference = 0.25,
    K = NULL,
    sel = NULL,
    seed = 1L,
    credible_level = 0.95,
    verbose = FALSE) {
  raw <- .gbkmr_get_raw(object)
  quantiles <- .gbkmr_validate_probabilities(quantiles, "quantiles")
  reference <- .gbkmr_validate_probabilities(reference, "reference")
  if (length(reference) != 1L) stop("reference must be a single probability.")
  credible_level <- .gbkmr_validate_level(credible_level)

  evaluated_quantiles <- unique(c(reference, quantiles))
  interventions <- .gbkmr_exposure_quantiles(raw, evaluated_quantiles)
  rownames(interventions) <- paste0("overall_q", seq_along(evaluated_quantiles))
  evaluated <- .gbkmr_gcompute(object, interventions, K, sel, seed, verbose)
  reference_index <- match(reference, evaluated_quantiles)

  summaries <- vector("list", length(quantiles))
  draw_tables <- vector("list", length(quantiles))
  for (i in seq_along(quantiles)) {
    regime_index <- match(quantiles[i], evaluated_quantiles)
    effect <- evaluated$draws[, regime_index] - evaluated$draws[, reference_index]
    interval <- .gbkmr_draw_interval(effect, credible_level)
    summaries[[i]] <- data.frame(
      quantile = quantiles[i],
      counterfactual_mean = mean(evaluated$draws[, regime_index]),
      estimate = interval[["estimate"]],
      lower = interval[["lower"]],
      upper = interval[["upper"]],
      reference_quantile = reference
    )
    draw_tables[[i]] <- data.frame(
      .draw = seq_len(nrow(evaluated$draws)),
      iteration = evaluated$iterations,
      quantile = quantiles[i],
      counterfactual_mean = evaluated$draws[, regime_index],
      reference_mean = evaluated$draws[, reference_index],
      effect = effect
    )
  }
  summary <- do.call(rbind, summaries)
  draws <- do.call(rbind, draw_tables)

  plot <- ggplot2::ggplot(summary, ggplot2::aes(
    x = .data$quantile, y = .data$estimate
  )) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "#5B8DB8", alpha = 0.2) +
    ggplot2::geom_line(color = "#215A7A", linewidth = 0.8) +
    ggplot2::geom_point(color = "#215A7A", size = 1.8) +
    ggplot2::labs(
      x = "Joint mixture quantile",
      y = paste0("Causal effect vs q = ", format(reference)),
      title = "Overall mixture causal effect"
    ) +
    .gbkmr_plot_theme()

  intervention_table <- data.frame(
    regime = rownames(interventions),
    quantile = evaluated_quantiles,
    as.data.frame(interventions, check.names = FALSE),
    check.names = FALSE
  )
  .gbkmr_new_causal_result(
    "overall", summary, draws, intervention_table, plot,
    paste0("E[Y(A=q)] - E[Y(A=", reference,
           ")], with every time-dependent confounder generated under its intervention history."),
    list(reference = reference, quantiles = quantiles, K = evaluated$K,
         sel = evaluated$iterations, credible_level = credible_level, seed = seed)
  )
}

#' Univariate causal dose-response curves
#'
#' Varies one exposure at a time, fixes all remaining exposures at a common
#' background quantile, and computes a full g-formula contrast against the
#' focal exposure's reference quantile.
#'
#' @inheritParams gbkmr_causal_overall
#' @param exposures Exposure column names to include. `NULL` uses all exposures.
#' @param time_points Optional integer visit indices used to filter exposures.
#' @param background Quantile at which non-focal exposures are fixed.
#'
#' @return A `gbkmr_causal_effects` object; see
#'   [gbkmr_causal_overall()].
#' @export
gbkmr_causal_univariate <- function(
    object,
    exposures = NULL,
    time_points = NULL,
    quantiles = seq(0.10, 0.90, by = 0.05),
    reference = 0.25,
    background = 0.50,
    K = NULL,
    sel = NULL,
    seed = 1L,
    credible_level = 0.95,
    verbose = FALSE) {
  raw <- .gbkmr_get_raw(object)
  selected <- .gbkmr_select_exposures(raw, exposures, time_points)
  quantiles <- .gbkmr_validate_probabilities(quantiles, "quantiles")
  reference <- .gbkmr_validate_probabilities(reference, "reference")
  background <- .gbkmr_validate_probabilities(background, "background")
  if (length(reference) != 1L || length(background) != 1L) {
    stop("reference and background must each be a single probability.")
  }
  credible_level <- .gbkmr_validate_level(credible_level)

  all_exposure_names <- raw$meta$exposure_names
  background_values <- .gbkmr_exposure_quantiles(raw, background)[1L, ]
  evaluated_quantiles <- unique(c(reference, quantiles))
  summary_tables <- draw_tables <- intervention_tables <- vector("list", length(selected))

  for (e in seq_along(selected)) {
    exposure <- selected[e]
    focal_values <- .gbkmr_exposure_quantiles(raw, evaluated_quantiles)[, exposure]
    interventions <- matrix(
      rep(background_values, each = length(evaluated_quantiles)),
      nrow = length(evaluated_quantiles),
      dimnames = list(NULL, all_exposure_names)
    )
    interventions[, exposure] <- focal_values
    rownames(interventions) <- paste0("univariate_", e, "_q", seq_along(evaluated_quantiles))
    evaluated <- .gbkmr_gcompute(object, interventions, K, sel, seed, verbose)
    reference_index <- match(reference, evaluated_quantiles)

    exposure_summaries <- exposure_draws <- vector("list", length(quantiles))
    for (i in seq_along(quantiles)) {
      regime_index <- match(quantiles[i], evaluated_quantiles)
      effect <- evaluated$draws[, regime_index] - evaluated$draws[, reference_index]
      interval <- .gbkmr_draw_interval(effect, credible_level)
      exposure_summaries[[i]] <- data.frame(
        exposure = exposure,
        quantile = quantiles[i],
        exposure_value = interventions[regime_index, exposure],
        counterfactual_mean = mean(evaluated$draws[, regime_index]),
        estimate = interval[["estimate"]],
        lower = interval[["lower"]],
        upper = interval[["upper"]],
        reference_quantile = reference,
        background_quantile = background
      )
      exposure_draws[[i]] <- data.frame(
        .draw = seq_len(nrow(evaluated$draws)),
        iteration = evaluated$iterations,
        exposure = exposure,
        quantile = quantiles[i],
        exposure_value = interventions[regime_index, exposure],
        counterfactual_mean = evaluated$draws[, regime_index],
        reference_mean = evaluated$draws[, reference_index],
        effect = effect
      )
    }
    summary_tables[[e]] <- do.call(rbind, exposure_summaries)
    draw_tables[[e]] <- do.call(rbind, exposure_draws)
    intervention_tables[[e]] <- data.frame(
      exposure = exposure,
      regime = rownames(interventions),
      quantile = evaluated_quantiles,
      as.data.frame(interventions, check.names = FALSE),
      check.names = FALSE
    )
  }

  summary <- do.call(rbind, summary_tables)
  draws <- do.call(rbind, draw_tables)
  intervention_table <- do.call(rbind, intervention_tables)
  plot <- ggplot2::ggplot(summary, ggplot2::aes(
    x = .data$exposure_value, y = .data$estimate
  )) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "#5B8DB8", alpha = 0.2) +
    ggplot2::geom_line(color = "#215A7A", linewidth = 0.8) +
    ggplot2::facet_wrap(ggplot2::vars(.data$exposure), scales = "free_x") +
    ggplot2::labs(
      x = "Exposure value",
      y = paste0("Causal effect vs focal q = ", format(reference)),
      title = "Univariate causal dose-response"
    ) +
    .gbkmr_plot_theme()

  .gbkmr_new_causal_result(
    "univariate", summary, draws, intervention_table, plot,
    paste0("E[Y(A_j=x, A_-j=q", background,
           ")] - E[Y(A_j=q", reference, ", A_-j=q", background,
           ")], using sequential g-computation."),
    list(exposures = selected, reference = reference, background = background,
         quantiles = quantiles, K = .gbkmr_or(K, raw$meta$K),
         sel = .gbkmr_or(sel, raw$meta$sel),
         credible_level = credible_level, seed = seed)
  )
}

.gbkmr_pairs <- function(selected, pairs, ordered = FALSE) {
  if (is.null(pairs)) {
    if (length(selected) < 2L) stop("At least two exposures are required.")
    if (ordered) {
      out <- expand.grid(
        focal = selected,
        conditional = selected,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      out <- as.matrix(out[out$focal != out$conditional, , drop = FALSE])
    } else {
      out <- t(utils::combn(selected, 2L))
    }
  } else if (is.data.frame(pairs) || is.matrix(pairs)) {
    out <- as.matrix(pairs)
  } else if (is.list(pairs)) {
    if (length(pairs) == 0L || any(lengths(pairs) != 2L)) {
      stop("Each element of pairs must contain two exposure names.")
    }
    out <- do.call(rbind, pairs)
  } else {
    stop("pairs must be a two-column matrix/data frame, a list, or NULL.")
  }
  if (ncol(out) != 2L) stop("pairs must have exactly two columns.")
  storage.mode(out) <- "character"
  colnames(out) <- c("focal", "conditional")
  unknown <- setdiff(unique(as.vector(out)), selected)
  if (length(unknown) > 0L) {
    stop("Pair exposures are not in the selected exposure set: ",
         paste(unknown, collapse = ", "))
  }
  if (any(out[, 1L] == out[, 2L])) stop("A bivariate pair must contain two different exposures.")
  unique(as.data.frame(out, stringsAsFactors = FALSE))
}

#' Conditional bivariate causal dose-response curves
#'
#' Varies a focal exposure while fixing a second exposure at requested
#' conditional quantiles. All remaining exposures are fixed at the background
#' quantile. Every curve is contrasted with one common intervention that fixes
#' the focal and conditional exposures at their respective reference quantiles.
#'
#' @inheritParams gbkmr_causal_univariate
#' @param pairs Exposure pairs. The first column or element is the focal
#'   exposure and the second is the conditional exposure. `NULL` uses all
#'   ordered pairs when `layout = "matrix"` and all unordered pairs when
#'   `layout = "wrap"`.
#' @param layout Plot layout. `matrix` uses focal exposures as columns and
#'   conditional exposures as rows; `wrap` uses one panel per pair.
#' @param conditional_quantiles Quantiles for the conditional exposure.
#' @param reference Focal exposure quantile in the common reference intervention.
#' @param conditional_reference Conditional exposure quantile in the common
#'   reference intervention.
#'
#' @return A `gbkmr_causal_effects` object; see
#'   [gbkmr_causal_overall()].
#' @export
gbkmr_causal_bivariate <- function(
    object,
    exposures = NULL,
    time_points = NULL,
    pairs = NULL,
    layout = c("matrix", "wrap"),
    quantiles = seq(0.10, 0.90, by = 0.05),
    conditional_quantiles = c(0.25, 0.50, 0.75),
    reference = 0.25,
    background = 0.50,
    K = NULL,
    sel = NULL,
    seed = 1L,
    credible_level = 0.95,
    verbose = FALSE,
    conditional_reference = 0.50) {
  raw <- .gbkmr_get_raw(object)
  selected <- .gbkmr_select_exposures(raw, exposures, time_points)
  layout <- match.arg(layout)
  pair_table <- .gbkmr_pairs(
    selected, pairs, ordered = is.null(pairs) && layout == "matrix"
  )
  quantiles <- .gbkmr_validate_probabilities(quantiles, "quantiles")
  conditional_quantiles <- .gbkmr_validate_probabilities(
    conditional_quantiles, "conditional_quantiles"
  )
  reference <- .gbkmr_validate_probabilities(reference, "reference")
  background <- .gbkmr_validate_probabilities(background, "background")
  conditional_reference <- .gbkmr_validate_probabilities(
    conditional_reference, "conditional_reference"
  )
  if (length(reference) != 1L || length(background) != 1L ||
      length(conditional_reference) != 1L) {
    stop(paste(
      "reference, background, and conditional_reference must each be",
      "a single probability."
    ))
  }
  credible_level <- .gbkmr_validate_level(credible_level)
  if (nrow(pair_table) > 20L) {
    warning("More than 20 bivariate panels requested. Use exposures, time_points, or pairs to reduce computation.")
  }

  all_names <- raw$meta$exposure_names
  background_values <- .gbkmr_exposure_quantiles(raw, background)[1L, ]
  evaluated_quantiles <- unique(c(reference, quantiles))
  evaluated_conditions <- unique(c(conditional_reference, conditional_quantiles))
  summary_tables <- draw_tables <- intervention_tables <- vector("list", nrow(pair_table))

  for (pair_index in seq_len(nrow(pair_table))) {
    focal <- pair_table$focal[pair_index]
    conditional <- pair_table$conditional[pair_index]
    specification <- expand.grid(
      focal_quantile = evaluated_quantiles,
      conditional_quantile = evaluated_conditions,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    interventions <- matrix(
      rep(background_values, each = nrow(specification)),
      nrow = nrow(specification),
      dimnames = list(NULL, all_names)
    )
    interventions[, focal] <- .gbkmr_exposure_quantiles(
      raw, specification$focal_quantile
    )[, focal]
    interventions[, conditional] <- .gbkmr_exposure_quantiles(
      raw, specification$conditional_quantile
    )[, conditional]
    rownames(interventions) <- paste0("bivariate_", pair_index, "_", seq_len(nrow(specification)))
    evaluated <- .gbkmr_gcompute(object, interventions, K, sel, seed, verbose)
    pair_label <- paste0(focal, " | ", conditional)
    reference_index <- which(
      specification$focal_quantile == reference &
        specification$conditional_quantile == conditional_reference
    )[1L]
    reference_focal_value <- interventions[reference_index, focal]
    reference_conditional_value <- interventions[reference_index, conditional]

    pair_summaries <- pair_draws <- list()
    output_index <- 0L
    for (condition in conditional_quantiles) {
      for (q in quantiles) {
        output_index <- output_index + 1L
        regime_index <- which(
          specification$focal_quantile == q &
            specification$conditional_quantile == condition
        )[1L]
        effect <- evaluated$draws[, regime_index] - evaluated$draws[, reference_index]
        interval <- .gbkmr_draw_interval(effect, credible_level)
        pair_summaries[[output_index]] <- data.frame(
          pair = pair_label,
          focal = focal,
          conditional = conditional,
          focal_quantile = q,
          focal_value = interventions[regime_index, focal],
          conditional_quantile = condition,
          conditional_value = interventions[regime_index, conditional],
          counterfactual_mean = mean(evaluated$draws[, regime_index]),
          estimate = interval[["estimate"]],
          lower = interval[["lower"]],
          upper = interval[["upper"]],
          reference_quantile = reference,
          conditional_reference_quantile = conditional_reference,
          reference_focal_value = reference_focal_value,
          reference_conditional_value = reference_conditional_value,
          background_quantile = background
        )
        pair_draws[[output_index]] <- data.frame(
          .draw = seq_len(nrow(evaluated$draws)),
          iteration = evaluated$iterations,
          pair = pair_label,
          focal = focal,
          conditional = conditional,
          focal_quantile = q,
          focal_value = interventions[regime_index, focal],
          conditional_quantile = condition,
          conditional_value = interventions[regime_index, conditional],
          counterfactual_mean = evaluated$draws[, regime_index],
          reference_mean = evaluated$draws[, reference_index],
          reference_quantile = reference,
          conditional_reference_quantile = conditional_reference,
          reference_focal_value = reference_focal_value,
          reference_conditional_value = reference_conditional_value,
          effect = effect
        )
      }
    }
    summary_tables[[pair_index]] <- do.call(rbind, pair_summaries)
    draw_tables[[pair_index]] <- do.call(rbind, pair_draws)
    intervention_tables[[pair_index]] <- data.frame(
      pair = pair_label,
      regime = rownames(interventions),
      is_reference = specification$focal_quantile == reference &
        specification$conditional_quantile == conditional_reference,
      specification,
      as.data.frame(interventions, check.names = FALSE),
      check.names = FALSE
    )
  }

  summary <- do.call(rbind, summary_tables)
  draws <- do.call(rbind, draw_tables)
  intervention_table <- do.call(rbind, intervention_tables)
  summary$condition_level <- factor(
    summary$conditional_quantile,
    levels = conditional_quantiles,
    labels = formatC(conditional_quantiles, format = "f", digits = 2)
  )
  plot_data <- summary
  plot_data$focal <- factor(plot_data$focal, levels = selected)
  plot_data$conditional <- factor(plot_data$conditional, levels = selected)
  conditional_colors <- if (length(conditional_quantiles) == 3L) {
    c("#D55E00", "#009E73", "#0072B2")
  } else {
    grDevices::hcl(
      seq(15, 375, length.out = length(conditional_quantiles) + 1L)[
        seq_along(conditional_quantiles)
      ],
      c = 80, l = 55
    )
  }

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$focal_value, y = .data$estimate,
    color = .data$condition_level, fill = .data$condition_level,
    group = .data$condition_level
  )) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      alpha = 0.20, color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_color_manual(values = conditional_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = conditional_colors, drop = FALSE)

  if (layout == "matrix") {
    plot <- plot +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$conditional),
        cols = ggplot2::vars(.data$focal),
        scales = "free_x",
        drop = FALSE
      )
  } else {
    plot <- plot +
      ggplot2::facet_wrap(ggplot2::vars(.data$pair), scales = "free_x")
  }

  plot <- plot +
    ggplot2::labs(
      x = "Focal exposure value",
      y = "Causal effect vs common reference",
      color = "Conditional quantile",
      fill = "Conditional quantile",
      title = "Conditional bivariate causal dose-response",
      subtitle = if (layout == "matrix") {
        "Columns: focal exposure; rows: conditional exposure"
      } else {
        NULL
      }
    ) +
    .gbkmr_plot_theme()

  if (layout == "matrix") {
    plot <- plot +
      ggplot2::theme(
        legend.position = "bottom",
        strip.text = ggplot2::element_text(size = 8)
      )
  }

  .gbkmr_new_causal_result(
    "bivariate", summary, draws, intervention_table, plot,
    paste0("E[Y(A_j=x, A_k=c, A_-(j,k)=q", background,
           ")] - E[Y(A_j=q", reference, ", A_k=q",
           conditional_reference, ", A_-(j,k)=q", background,
           ")], with time-dependent confounders generated sequentially."),
    list(pairs = pair_table, layout = layout, reference = reference,
         background = background, conditional_reference = conditional_reference,
         quantiles = quantiles, conditional_quantiles = conditional_quantiles,
         K = .gbkmr_or(K, raw$meta$K), sel = .gbkmr_or(sel, raw$meta$sel),
         credible_level = credible_level, seed = seed)
  )
}

#' Single-exposure IQR causal effects
#'
#' For each focal exposure, estimates its interquartile-range causal effect
#' while jointly fixing all remaining exposures at each requested background
#' quantile.
#'
#' @inheritParams gbkmr_causal_univariate
#' @param contrast_quantiles Two quantiles defining the focal exposure contrast.
#' @param background_quantiles Quantiles at which all non-focal exposures are
#'   jointly fixed.
#'
#' @return A `gbkmr_causal_effects` object; see
#'   [gbkmr_causal_overall()].
#' @export
gbkmr_causal_iqr <- function(
    object,
    exposures = NULL,
    time_points = NULL,
    contrast_quantiles = c(0.25, 0.75),
    background_quantiles = c(0.25, 0.50, 0.75),
    K = NULL,
    sel = NULL,
    seed = 1L,
    credible_level = 0.95,
    verbose = FALSE) {
  raw <- .gbkmr_get_raw(object)
  selected <- .gbkmr_select_exposures(raw, exposures, time_points)
  contrast_quantiles <- .gbkmr_validate_probabilities(
    contrast_quantiles, "contrast_quantiles", min_length = 2L
  )
  if (length(contrast_quantiles) != 2L || contrast_quantiles[1L] >= contrast_quantiles[2L]) {
    stop("contrast_quantiles must contain two increasing probabilities.")
  }
  background_quantiles <- .gbkmr_validate_probabilities(
    background_quantiles, "background_quantiles"
  )
  credible_level <- .gbkmr_validate_level(credible_level)
  all_names <- raw$meta$exposure_names
  summary_tables <- draw_tables <- intervention_tables <- vector("list", length(selected))

  for (e in seq_along(selected)) {
    exposure <- selected[e]
    specification <- expand.grid(
      focal_quantile = contrast_quantiles,
      background_quantile = background_quantiles,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    interventions <- matrix(NA_real_, nrow(specification), length(all_names),
                            dimnames = list(NULL, all_names))
    for (i in seq_len(nrow(specification))) {
      interventions[i, ] <- .gbkmr_exposure_quantiles(
        raw, specification$background_quantile[i]
      )[1L, ]
      interventions[i, exposure] <- .gbkmr_exposure_quantiles(
        raw, specification$focal_quantile[i]
      )[1L, exposure]
    }
    rownames(interventions) <- paste0("iqr_", e, "_", seq_len(nrow(specification)))
    evaluated <- .gbkmr_gcompute(object, interventions, K, sel, seed, verbose)

    exposure_summaries <- exposure_draws <- vector("list", length(background_quantiles))
    for (i in seq_along(background_quantiles)) {
      background_q <- background_quantiles[i]
      low_index <- which(
        specification$focal_quantile == contrast_quantiles[1L] &
          specification$background_quantile == background_q
      )[1L]
      high_index <- which(
        specification$focal_quantile == contrast_quantiles[2L] &
          specification$background_quantile == background_q
      )[1L]
      effect <- evaluated$draws[, high_index] - evaluated$draws[, low_index]
      interval <- .gbkmr_draw_interval(effect, credible_level)
      exposure_summaries[[i]] <- data.frame(
        exposure = exposure,
        low_quantile = contrast_quantiles[1L],
        high_quantile = contrast_quantiles[2L],
        background_quantile = background_q,
        estimate = interval[["estimate"]],
        lower = interval[["lower"]],
        upper = interval[["upper"]]
      )
      exposure_draws[[i]] <- data.frame(
        .draw = seq_len(nrow(evaluated$draws)),
        iteration = evaluated$iterations,
        exposure = exposure,
        low_quantile = contrast_quantiles[1L],
        high_quantile = contrast_quantiles[2L],
        background_quantile = background_q,
        low_mean = evaluated$draws[, low_index],
        high_mean = evaluated$draws[, high_index],
        effect = effect
      )
    }
    summary_tables[[e]] <- do.call(rbind, exposure_summaries)
    draw_tables[[e]] <- do.call(rbind, exposure_draws)
    intervention_tables[[e]] <- data.frame(
      exposure = exposure,
      regime = rownames(interventions),
      specification,
      as.data.frame(interventions, check.names = FALSE),
      check.names = FALSE
    )
  }

  summary <- do.call(rbind, summary_tables)
  draws <- do.call(rbind, draw_tables)
  intervention_table <- do.call(rbind, intervention_tables)
  plot_data <- summary
  plot_data$exposure_plot <- factor(plot_data$exposure, levels = rev(selected))
  plot_data$background_level <- factor(
    plot_data$background_quantile,
    levels = background_quantiles,
    labels = paste0(round(100 * background_quantiles), "th")
  )
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(
    y = .data$exposure_plot, x = .data$estimate, color = .data$background_level
  )) +
    ggplot2::geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      orientation = "y", width = 0.25,
      position = ggplot2::position_dodge(width = 0.65)
    ) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.65)) +
    ggplot2::labs(
      x = paste0("IQR causal effect (posterior mean and ",
                 round(100 * credible_level), "% CrI)"),
      y = "Exposure",
      color = "Other exposures",
      title = "Single-exposure IQR causal effects"
    ) +
    .gbkmr_plot_theme()

  .gbkmr_new_causal_result(
    "iqr", summary, draws, intervention_table, plot,
    paste0("E[Y(A_j=q", contrast_quantiles[2L], ", A_-j=q_b)] - ",
           "E[Y(A_j=q", contrast_quantiles[1L], ", A_-j=q_b)], ",
           "with sequential g-computation at each background q_b."),
    list(exposures = selected, contrast_quantiles = contrast_quantiles,
         background_quantiles = background_quantiles,
         K = .gbkmr_or(K, raw$meta$K), sel = .gbkmr_or(sel, raw$meta$sel),
         credible_level = credible_level, seed = seed)
  )
}

#' High-versus-low interaction contrasts
#'
#' Computes the posterior difference between each focal exposure's IQR effect
#' when all remaining exposures are high and the corresponding IQR effect when
#' all remaining exposures are low.
#'
#' @inheritParams gbkmr_causal_iqr
#' @param low_background Background quantile defining the low-mixture setting.
#' @param high_background Background quantile defining the high-mixture setting.
#'
#' @return A `gbkmr_causal_effects` object; see
#'   [gbkmr_causal_overall()].
#' @export
gbkmr_causal_interaction <- function(
    object,
    exposures = NULL,
    time_points = NULL,
    contrast_quantiles = c(0.25, 0.75),
    low_background = 0.25,
    high_background = 0.75,
    K = NULL,
    sel = NULL,
    seed = 1L,
    credible_level = 0.95,
    verbose = FALSE) {
  low_background <- .gbkmr_validate_probabilities(low_background, "low_background")
  high_background <- .gbkmr_validate_probabilities(high_background, "high_background")
  if (length(low_background) != 1L || length(high_background) != 1L ||
      low_background >= high_background) {
    stop("low_background and high_background must be increasing scalar probabilities.")
  }

  iqr_result <- gbkmr_causal_iqr(
    object = object,
    exposures = exposures,
    time_points = time_points,
    contrast_quantiles = contrast_quantiles,
    background_quantiles = c(low_background, high_background),
    K = K,
    sel = sel,
    seed = seed,
    credible_level = credible_level,
    verbose = verbose
  )

  selected <- unique(iqr_result$summary$exposure)
  summary_tables <- draw_tables <- vector("list", length(selected))
  for (e in seq_along(selected)) {
    exposure <- selected[e]
    low_draws <- iqr_result$draws[
      iqr_result$draws$exposure == exposure &
        iqr_result$draws$background_quantile == low_background,
    ]
    high_draws <- iqr_result$draws[
      iqr_result$draws$exposure == exposure &
        iqr_result$draws$background_quantile == high_background,
    ]
    low_draws <- low_draws[order(low_draws$.draw), ]
    high_draws <- high_draws[order(high_draws$.draw), ]
    contrast <- high_draws$effect - low_draws$effect
    interval <- .gbkmr_draw_interval(contrast, credible_level)
    summary_tables[[e]] <- data.frame(
      exposure = exposure,
      low_background = low_background,
      high_background = high_background,
      estimate = interval[["estimate"]],
      lower = interval[["lower"]],
      upper = interval[["upper"]],
      excludes_zero = interval[["lower"]] > 0 | interval[["upper"]] < 0
    )
    draw_tables[[e]] <- data.frame(
      .draw = low_draws$.draw,
      iteration = low_draws$iteration,
      exposure = exposure,
      iqr_effect_low = low_draws$effect,
      iqr_effect_high = high_draws$effect,
      interaction_contrast = contrast
    )
  }
  summary <- do.call(rbind, summary_tables)
  draws <- do.call(rbind, draw_tables)
  plot_data <- summary
  plot_data$exposure_plot <- factor(plot_data$exposure, levels = rev(selected))
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(
    y = .data$exposure_plot, x = .data$estimate
  )) +
    ggplot2::geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      orientation = "y", width = 0.2, color = "#215A7A"
    ) +
    ggplot2::geom_point(size = 2, color = "#215A7A") +
    ggplot2::labs(
      x = "IQR effect at high mixture - IQR effect at low mixture",
      y = "Exposure",
      title = "Interaction contrasts"
    ) +
    .gbkmr_plot_theme()

  .gbkmr_new_causal_result(
    "interaction", summary, draws, iqr_result$interventions, plot,
    paste0("IQR effect with A_-j at q", high_background,
           " minus the IQR effect with A_-j at q", low_background, "."),
    list(exposures = selected, contrast_quantiles = contrast_quantiles,
         low_background = low_background, high_background = high_background,
         K = iqr_result$settings$K, sel = iqr_result$settings$sel,
         credible_level = credible_level, seed = seed)
  )
}

#' @export
print.gbkmr_causal_effects <- function(x, ...) {
  print(x$plot)
  invisible(x)
}

#' @export
plot.gbkmr_causal_effects <- function(x, ...) {
  print(x$plot)
  invisible(x$plot)
}

#' @export
summary.gbkmr_causal_effects <- function(object, ...) {
  object$summary
}

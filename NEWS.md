# causalBKMR 0.1.1

## Bug fixes

* Monte Carlo g-computation now draws continuous time-varying confounders from
  the fitted conditional distribution, `h_j(z) + x'beta_j + sigma_j * eps`,
  instead of at its conditional mean. The missing residual term made the K
  Monte Carlo samples of each confounder identical at the first post-baseline
  visit and understated confounder variability afterwards, which biased the
  ACE when the outcome surface is non-linear in the confounder and narrowed its
  credible interval (Chai et al. scenario 3, 50 paired replicates, 50 knots:
  coverage of the 95% interval 0.56 before the fix, 0.86 after).
  Affects `run_gbkmr_panel()`, `gbkmr_run()` and all `gbkmr_causal_*()` plots.
* Binary time-varying confounders are now drawn with K independent Bernoulli
  uniforms per posterior draw; previously the seed reset before each Monte
  Carlo sample made the K draws identical whenever the kernel inputs were.
* The `gbkmr_causal_*()` functions used an independent posterior-function draw
  for every Monte Carlo sample, averaging away posterior uncertainty. They now
  reset the seed before every `SamplePred()` call so the K samples of a
  posterior draw use common random numbers, as `gbkmr_run()` always did (and as
  in Chai et al.'s per-sample seed reset).
  Their confounder noise is likewise shared across the regimes evaluated in one
  call, so a regime's draw does not depend on which other regimes are
  requested (`gbkmr_causal_interaction()` therefore equals the difference of
  the matching `gbkmr_causal_iqr()` effects exactly).

## New features

* `n_knots = NULL` (or 0) fits the exact Gaussian process with no knots, as in
  Chai et al. The default remains 50 knots; note that the predictive-process
  approximation narrows posterior intervals.
* `n_cores` now also applies to the standard BKMR engine: the BKMR models are
  fit in parallel and the Monte Carlo g-computation is split across posterior
  draws with `parallel::mclapply()`. Results are identical for any `n_cores`.
  Default is 1 (serial); forking is unavailable on Windows.
* Identical kernel-input rows within a posterior draw are evaluated once, which
  removes K-1 redundant `SamplePred()` calls at the first post-baseline visit.

## Behaviour changes

* Numerical results differ from 0.1.0 because of the fixes above and because
  all Monte Carlo seeds (model fits, posterior-function draws, confounder
  noise) are now derived from `currind` through one collision-free scheme.
* `run_gbkmr_panel(n_cores = )` defaults to `NULL`: 1 for `engine = "bkmr"`
  (previously 10, but unused) and 10 for `engine = "fastbkmr"` (unchanged).
* `run_gbkmr_panel(verbose_every = )` is ignored; progress is reported once per
  model fit and per sampling stage.
* `gbkmr_causal_*()` fork as many workers as the fit used (`meta$mc_cores`);
  set `options(causalBKMR.cores = 1)` to run them serially.
* New fields: `raw_results$meta$mc_cores`, `raw_results$meta$confounder_noise`,
  `call_info$n_knots`, `call_info$n_cores`.

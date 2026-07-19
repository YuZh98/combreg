# combreg 0.2.0

- `crr_control()` gains an adaptive latent-utility block size. `zeta_block` now
  accepts `"adaptive"` (the new default) in addition to a fixed integer. The
  adaptive rule tunes the block size during warmup with a dual-averaging
  controller to hit `zeta_target_accept` (default 0.6), then freezes it for
  sampling. The optimal block size depends on the constraint geometry rather
  than the response dimension alone, so this removes a manual tuning step and
  improves mixing efficiency across problems (up to ~2x faster than the former
  fixed `min(d, 100)` on tightly constrained, high-dimensional responses) with
  no change to the posterior. Pass an integer `zeta_block` for the previous
  fixed behavior.
- **Reproducibility note:** because the default changed from a fixed
  `min(d, 100)` block to `"adaptive"`, a seeded `crr()` call with no explicit
  `zeta_block` now consumes a different RNG stream and returns different draws
  than in 0.1.0. Pass `zeta_block = 100` (or your prior value) to reproduce
  0.1.0 results exactly; `inst/scripts/reproduce-table2.R` already does this.
- The tuned block size is stored on the fit (`zeta_block_tuned`) and shown by
  `print()`.
- `accept_rate` (and the low-acceptance diagnostic) now averages over the
  post-warmup sampling phase only, so it reflects the frozen proposal rather
  than warmup adaptation transients.
- `crr()` accepts a one-sided formula and `data` in place of a covariate
  matrix, building the design matrix via `model.matrix()` (factor expansion,
  interactions, intercept). The matrix interface is unchanged.
- `crr_constraints()` gains a `dedup` argument (default `TRUE`) that removes
  redundant constraint rows: rows identical in `A` are collapsed to the one
  with the smallest `b`, since a tighter `A y <= b_1` implies any looser
  `A y <= b_2`. This leaves the feasible set and posterior unchanged, preserves
  total unimodularity, and speeds the sampler by dropping redundant dual
  dimensions. `random_constraints()` uses `dedup = FALSE` so simulated designs
  keep their exact `(d, m)`.

# combreg 0.1.0

Initial release.

- `crr()`: MH-within-Gibbs sampler for combinatorial response regression
  with dual-certificate augmentation (exponential and half-Gaussian dual
  kernels), plus an unconstrained Albert–Chib probit baseline.
- Constraint utilities: `crr_constraints()`, `is_tum()`, `is_feasible()`.
- Sampling primitives exported for custom models: `init_dual()`,
  `sample_dual()`, `draw_utility()`, `sample_utility()`, `dual_feasible()`,
  `coef_precompute()`, `update_coef()`.
- Data simulation: `simulate_crr()`, `random_constraints()`.
- Diagnostics: `crr_diagnostics()` structured MCMC + regression report
  (including posterior predictive goodness-of-fit checks via `crr_ppc()`
  and, with known truth, coefficient RMSE and interval coverage),
  `crr_ess()`, `crr_rhat()`, `fitted()` and `residuals()` methods.
- Benchmarking: `crr_benchmark()` compares methods on common data (timing,
  ESS per second, RMSE, and interval coverage against known truth).
- Plotting: `plot()` types `"trace"`, `"acf"`, `"violin"`, `"ess"`,
  `"ess_time"`, `"residual"` (heat map), and `"beta_diff"` (coefficient
  error heat map against known truth).
- `coda` and `posterior` converters; `summary()`, `coef()`, `predict()`
  methods.
- Vignettes: introduction (`combreg`) and diagnostics/benchmarking
  (`diagnostics`).
- `inst/scripts/reproduce-table2.R`: reproduces the paper's Table 2
  simulation study at the original settings.

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

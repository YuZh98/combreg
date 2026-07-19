# combreg

<!-- badges: start -->
[![R-CMD-check](https://github.com/YuZh98/combreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YuZh98/combreg/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Bayesian regression for **combinatorial response data**: each response is a
binary vector constrained to an integral polytope,
$\mathcal{Y} = \{\,y \in \{0,1\}^d : A y \le b\,\}$. Ignoring the
constraint and fitting an unconstrained model leads to biased estimation, and enumerating all feasible outcomes is computationally impractical in high-dimensional settings. In [Zheng, Ghosh & Duan (2026+)](https://arxiv.org/abs/2504.11630), an augmented likelihood that respects the combinatorial constraints is proposed, along with a Metropolis–Hastings-within-Gibbs (MH-Within-Gibbs) sampler that scales to high-dimensional problems.

The package implements the MH-Within-Gibbs sampler of [Zheng, Ghosh & Duan (2026+)](https://arxiv.org/abs/2504.11630), with hit-and-run
dual updates in C++ (OpenMP-parallel across observations), an unconstrained probit baseline, one-call benchmarking, MCMC and regression diagnostics, and utilities for constraint validation and data simulation.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("YuZh98/combreg")
```

Building from source requires a C++17 toolchain (the package links to
`RcppArmadillo`); OpenMP is used when available and falls back to a serial build
otherwise.

## Usage

```r
library(combreg)

# a small matching-type constraint: A y <= b on y in {0,1}^3
A <- rbind(c(1, 1, 0),
           c(0, 1, 1))
con <- crr_constraints(A, b = c(1, 1))

# simulate constrained responses, then fit
sim <- simulate_crr(n = 300, p = 2, constraints = con, seed = 1)
fit <- crr(sim$Y, sim$X, con, n_iter = 2000, warmup = 1000, chains = 2, seed = 1)
summary(fit)

# alternatively, predictors can also be given as a formula + data frame
df  <- data.frame(x1 = sim$X[, 1], x2 = sim$X[, 2])
fit <- crr(sim$Y, ~ 0 + x1 + x2, con, data = df,
           n_iter = 2000, warmup = 1000, chains = 2, seed = 1)

# diagnostics, plots, and interoperability
crr_diagnostics(fit)                 # structured MCMC + regression report
plot(fit, type = "trace")            # also: "acf", "violin", "ess", "ess_time", "residual"
coda::as.mcmc(fit)                   # -> coda ecosystem
posterior::as_draws(fit)             # -> posterior / bayesplot ecosystem
predict(fit, type = "response")

# benchmark the constrained sampler against the unconstrained probit baseline
crr_benchmark(sim$Y, sim$X, con, beta = sim$beta,
              n_iter = 2000, warmup = 1000, seed = 1)
```

## Overview

| Component | Functions |
|---|---|
| Model fitting | `crr()`, `crr_control()`, `crr_prior()` |
| Constraints | `crr_constraints()`, `is_tum()`, `is_feasible()` |
| Simulation | `simulate_crr()`, `random_constraints()` |
| Diagnostics | `crr_diagnostics()`, `crr_ess()`, `crr_rhat()`, `crr_ppc()` |
| Benchmarking | `crr_benchmark()` |
| Methods | `summary()`, `coef()`, `predict()`, `fitted()`, `residuals()`, `plot()`, `as.mcmc()`, `as_draws()` |
| Sampling primitives | `init_dual()`, `sample_dual()`, `draw_utility()`, `sample_utility()`, `dual_feasible()`, `coef_precompute()`, `update_coef()` |

Two vignettes and a R script cover the details:

- `vignette("combreg")`: the model, why constraint-ignoring baselines are
  biased, constraint setup, the adaptive block controller, and building custom models from the exported primitives.
- `vignette("diagnostics")`: the diagnostics and benchmarking workflow.
- `inst/scripts/reproduce-table2.R`: runs the simulation in [Zheng, Ghosh & Duan (2026+)](https://arxiv.org/abs/2504.11630) and reproduces its Table 2 (coefficient RMSE for varying `(d, m)`).



## Citation

```r
citation("combreg")
```

Zheng, Y., Ghosh, M., & Duan, L. (2026+). *Statistical Modeling of Combinatorial Response Data.* arXiv:2504.11630. <https://doi.org/10.48550/arXiv.2504.11630>

## License

MIT © Hugh Zheng

# combreg

Bayesian regression for combinatorial response data: responses are binary
vectors constrained to an integral polyhedron $\{y\in\{0,1\}^d:Ay\le b\}$.

The package implements the Metropolis–Hastings-within-Gibbs sampler with dual-certificate augmentation from [Statistical Modeling of Combinatorial Response Data](https://arxiv.org/abs/2504.11630), with hit-and-run dual updates in C++ (OpenMP-parallel across observations), an unconstrained probit baseline, one-call benchmarking, MCMC and regression diagnostics, and utilities for constraint validation and data simulation.

## Installation

```r
# install.packages("remotes")
remotes::install_github("YuZh98/combreg")
```

## Quick start

```r
library(combreg)

A <- rbind(c(1, 1, 0),
           c(0, 1, 1))
con <- crr_constraints(A, b = c(1, 1))

sim <- simulate_crr(n = 300, p = 2, constraints = con, seed = 1)

fit <- crr(sim$Y, sim$X, con, n_iter = 2000, warmup = 1000, chains = 2, seed = 1)
summary(fit)

# predictors can also be given as a formula + data frame:
df  <- data.frame(x1 = sim$X[, 1], x2 = sim$X[, 2])
fit <- crr(sim$Y, ~ 0 + x1 + x2, con, data = df,
           n_iter = 2000, warmup = 1000, chains = 2, seed = 1)
crr_diagnostics(fit)    # structured MCMC + regression diagnostics report

plot(fit, type = "trace")     # also: "acf", "violin", "ess", "ess_time",
plot(fit, type = "residual")  #       residual heat map

coda::as.mcmc(fit)      # further MCMC diagnostics via coda
predict(fit, type = "response")

# benchmark the paper's sampler against the unconstrained probit baseline
crr_benchmark(sim$Y, sim$X, con, beta = sim$beta,
              n_iter = 2000, warmup = 1000, seed = 1)
```

For a cross-package comparison against constraint-ignoring alternatives
(`bayesm::rmvpGibbs` multivariate probit, `MNP::mnp` on the simplex special
case) on RMSE, interval coverage, response feasibility, and sampling
efficiency, see `inst/benchmarks/benchmark_crosspkg.R`.

See `vignette("combreg")` for the model, the bias of constraint-ignoring
baselines, and how to build custom models from the exported sampling
primitives, and `vignette("diagnostics")` for the diagnostics and
benchmarking workflow.

## Reproducing the paper's simulation results

Table 2 of the paper (RMSE for estimating the coefficients across
dimensionalities `(d, m)`) can be reproduced with the shipped script, which
runs the paper's exact settings (`n = 1000`, `p = 5`, 50000 iterations,
warmup 5000, thinning 25, 100 inner hit-and-run steps):

```sh
Rscript "$(Rscript -e 'cat(system.file("scripts", "reproduce-table2.R", package = "combreg"))')" 2,1 5,1 10,5
```

Low-dimensional cells take minutes; the largest cells (`d = 1000`) take
hours to days, matching the timings reported in the paper.

## License

MIT © Hugh Zheng

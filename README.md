# combreg

Bayesian regression for combinatorial response data: responses are binary
vectors constrained to a polytope `{y ∈ {0,1}^d : Ay ≤ b}` with totally
unimodular `A` (matchings, assignments, and related structures).

The package implements the Metropolis–Hastings-within-Gibbs sampler with
dual-certificate augmentation from *Statistical Modeling of Combinatorial
Response Data*, with hit-and-run dual updates in C++ (OpenMP-parallel across
observations), an unconstrained probit baseline, and utilities for constraint
validation and data simulation.

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

fit <- crr(sim$Y, sim$X, con, n_iter = 2000, warmup = 1000, seed = 1)
summary(fit)
plot(fit)

coda::as.mcmc(fit)      # MCMC diagnostics via coda
predict(fit, type = "response")
```

See `vignette("combreg")` for the model, the bias of constraint-ignoring
baselines, and how to build custom models from the exported sampling
primitives.

## License

MIT © Yu Zheng

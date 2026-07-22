# combreg — code map

Bayesian regression for combinatorial response data. The package is organized as a three-layer stack: native samplers at the base, exported R primitives in the middle, and model orchestration plus diagnostics on top.

| | | |
|---|---|---|
| **v0.2.0** MH-within-Gibbs | **20** exported functions | **4** C++ entry points · **7** S3 classes |

## Architecture

```mermaid
flowchart TD
    subgraph L2["Layer 2 · orchestration & analysis"]
        crr["crr.R<br/>crr() · fit_mhwg · fit_unconstrained"]
        methods["methods.R<br/>print/summary/coef/plot/predict · as.mcmc/as_draws"]
        diag["diagnostics.R<br/>crr_diagnostics · crr_ess · crr_rhat · crr_ppc"]
        bench["benchmark.R<br/>crr_benchmark()"]
        sim["simulate.R<br/>simulate_crr · random_constraints · ilp_argmax"]
    end
    subgraph L1["Layer 1 · primitives & constructors"]
        prim["primitives.R<br/>init_dual · sample_dual · draw_utility<br/>dual_feasible · sample_utility · coef_precompute · update_coef"]
        cons["constraints.R<br/>crr_constraints · is_tum · is_feasible"]
        cfg["priors.R · control.R<br/>crr_prior · crr_control"]
    end
    subgraph L0["Layer 0 · native kernels (Rcpp / RcppArmadillo, OpenMP)"]
        har["hit_and_run.cpp<br/>hit_and_run_cpp · loop_hit_and_run_cpp · check_feasible_dual_cpp"]
        tum["tum_check.cpp<br/>check_tum_cpp"]
    end

    crr --> prim
    bench --> crr
    bench --> diag
    diag --> methods
    sim --> cons
    methods --> prim
    prim --> har
    cons --> tum
```

## The three layers

### Layer 2: orchestration & analysis
User-facing entry points that compose the primitives.

| File | Role | Exported | Internal |
|---|---|---|---|
| `crr.R` | main entry (covariate matrix **or one-sided formula + `data`**); dispatch to a sampler; warmup-phase adaptive block controller | `crr` | `fit_mhwg`, `fit_unconstrained`, `.zeta_block_controller` |
| `methods.R` | S3 methods on `crr_fit` | `print`, `summary`, `coef`, `plot`, `predict`, `as.mcmc`, `as_draws` | `kept_draws`, `plot_*` helpers |
| `diagnostics.R` | MCMC + regression diagnostics | `crr_diagnostics`, `crr_ess`, `crr_rhat`, `crr_ppc`, `fitted`, `residuals` | `split_rhat`, `ppc_table` |
| `benchmark.R` | compare methods on shared data | `crr_benchmark` | — |
| `simulate.R` | synthetic data + random constraint systems | `simulate_crr`, `random_constraints` | `ilp_argmax` (needs lpSolve), `require_lpsolve` |

### Layer 1: sampling primitives & constructors
Exported building blocks; reusable for custom mean structures.

| File | Role | Exported | Calls C++ |
|---|---|---|---|
| `primitives.R` | one Gibbs / MH block each | `init_dual`, `sample_dual`, `draw_utility`, `dual_feasible`, `sample_utility`, `coef_precompute`, `update_coef` | `sample_dual` → `loop_hit_and_run_cpp`; `dual_feasible` → `check_feasible_dual_cpp` |
| `constraints.R` | feasible-set definition, redundant-row **de-duplication** (`dedup`), TUM check | `crr_constraints`, `is_tum`, `is_feasible` | `crr_constraints`, `is_tum` → `check_tum_cpp` |
| `priors.R`, `control.R` | validated parameter objects; `crr_control` tunes the **adaptive `zeta_block`** and the two hit-and-run sweep counts (`n_iter_hit_and_run` for the dual refresh, `n_iter_hit_and_run_mh` for the MH certificate) | `crr_prior`, `crr_control` | — |

### Layer 0: native kernels
Rcpp / RcppArmadillo, OpenMP-parallel, no R state.

| File | Role | Entry points |
|---|---|---|
| `hit_and_run.cpp` | dual-polytope hit-and-run; per-chain RNG streams (thread-count invariant) | `hit_and_run_cpp`, `loop_hit_and_run_cpp` (OpenMP over observations), `check_feasible_dual_cpp` |
| `tum_check.cpp` | exhaustive total-unimodularity test | `check_tum_cpp` |
| `R/RcppExports.R`, `src/RcppExports.cpp` | generated R ↔ C++ bridge | — |

## The sampler, traced

What `crr(method = "mhwg")` runs each iteration:

1. **Dual refresh** — `U | ζ, y` — hit-and-run over the dual-certificate polytope, one chain per observation.
2. **Utility MH** — `ζ | U, β, y` — propose a coordinate block, form the envelope, accept per row by feasibility.
3. **Conjugate draw** — `β | ζ` — Gaussian update via a cached Cholesky factor.

```
crr()  →  fit_mhwg()                                   [Layer 2]  seed, allocate draws, loop n_iter
  1 ·  sample_dual()      → loop_hit_and_run_cpp        [L1→L0]   draw U
  2 ·  sample_utility()                                 [L1]      orchestrates the MH step:
         → draw_utility()                               [L1]      propose zeta (for a block of coordinates)
         → sample_dual()  → loop_hit_and_run_cpp        [L1→L0]   draw U* 
         → dual_feasible() → check_feasible_dual_cpp    [L1→L0]   decide if accept or reject
  3 ·  update_coef()                                    [L1]      uses coef_precompute() factor (V, L)
```

During warmup the step-2 block size is retuned by a dual-averaging controller (`.zeta_block_controller`) toward `zeta_target_accept`, then frozen for sampling. `fit_unconstrained` (the
Albert–Chib probit baseline) shares `coef_precompute`, `draw_utility`, and
`update_coef` but skips the dual machinery entirely.

## Around the edges

**External dependencies**

| Package | Used for | Requirement |
|---|---|---|
| Rcpp / RcppArmadillo | Armadillo kernels, `useDynLib` | Imports / LinkingTo |
| coda | `as.mcmc`, `effectiveSize` (ESS) | Imports |
| truncnorm | `rtruncnorm` in `draw_utility` | Imports |
| stats, graphics, grDevices, utils | `chol`, `solve`, `rnorm`, `rexp`, `acf`, plotting | Imports |
| lpSolve | ILP argmax — `simulate_crr` & `predict(type = "response")` | Suggests |
| posterior | `as_draws` converter | Suggests |

**S3 classes:**

 `crr_fit` (the model), `crr_constraints` (feasible set),
`crr_prior`, `crr_control`, `crr_diagnostics`, `crr_ppc` (predictive checks),
`crr_benchmark`.

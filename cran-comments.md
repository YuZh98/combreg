## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* local macOS, R 4.4.3
* win-builder, R-devel (2026-07-17 r90265 ucrt) — 0 errors, 0 warnings, 1 note

## Notes

This is a new release, so a "New submission" NOTE is expected.

On win-builder the incoming-feasibility NOTE also flags possibly-misspelled
words in DESCRIPTION: these are technical terms ("combinatorial",
"Combinatorial", "unimodularity", "Rhat"), an author surname ("Zheng"), and the
Latin "et al." — all spelled correctly.

The method reference is given as a DOI, <doi:10.48550/arXiv.2504.11630>.

## Compiled code

The package contains C++ (Rcpp / RcppArmadillo) with optional OpenMP
parallelism. OpenMP is enabled by a `configure` probe that falls back to a
serial build where it is unavailable; Windows uses `SHLIB_OPENMP_CXXFLAGS` via
`src/Makevars.win`. No more than two threads are used in any example, test, or
vignette.

## Downstream dependencies

None; this is a new package.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* local macOS, R 4.4.3
* win-builder (R-devel) — pending
* macOS builder, R-hub — pending

## Notes

The only NOTE on CRAN's machines is the expected "New submission".

Local checks additionally emit environment-specific messages that do not arise
on CRAN's check infrastructure: missing `qpdf` and `checkbashisms` helper tools,
and an outdated system HTML `tidy` that rejects HTML5 tags in the rendered
manual. The `configure` script is POSIX `sh` with no bashisms.

## Compiled code

The package contains C++ (Rcpp / RcppArmadillo) with optional OpenMP
parallelism. OpenMP is enabled by a `configure` probe that falls back to a
serial build where it is unavailable; Windows uses `SHLIB_OPENMP_CXXFLAGS` via
`src/Makevars.win`. No more than two threads are used in any example, test, or
vignette.

## Downstream dependencies

None; this is a new package.

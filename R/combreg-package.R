#' @keywords internal
"_PACKAGE"

#' @useDynLib combreg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rnorm rexp acf coef predict fitted residuals qnorm
#' @importFrom truncnorm rtruncnorm
NULL

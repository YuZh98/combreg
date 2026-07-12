#' Prior specification for crr models
#'
#' Currently a single family: independent Gaussian priors
#' \eqn{\beta_{kj} \sim N(0, sd^2)} on all regression coefficients.
#'
#' @param family Prior family; only `"gaussian"` is implemented.
#' @param sd Prior standard deviation (scalar, positive).
#'
#' @return An object of class `crr_prior`.
#'
#' @examples
#' crr_prior()               # N(0, 1)
#' crr_prior(sd = 10)        # weakly informative
#'
#' @export
crr_prior <- function(family = "gaussian", sd = 1) {
  family <- match.arg(family, "gaussian")
  if (!is.numeric(sd) || length(sd) != 1 || sd <= 0) {
    stop("sd must be a positive scalar", call. = FALSE)
  }
  structure(
    list(family = family, sd = sd),
    class = "crr_prior"
  )
}

#' @export
print.crr_prior <- function(x, ...) {
  cat("<crr_prior>: beta_kj ~ N(0, ", format(x$sd^2), ")\n", sep = "")
  invisible(x)
}

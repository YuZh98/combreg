#' Sampler control parameters
#'
#' Tuning parameters for the MH-within-Gibbs sampler in [crr()].
#'
#' @param n_iter_hit_and_run Inner hit-and-run steps per dual update. The
#'   default (100) matches the simulation study of the paper; smaller values
#'   are faster but mix the dual variable less well.
#' @param rho Rate of the exponential dual kernel (ignored for
#'   `kernel = "half_gaussian"`).
#' @param zeta_block Maximum number of response coordinates updated per
#'   MH step for the latent utilities (the paper uses `min(d, 100)`).
#' @param max_dir_tries Maximum attempts to draw a valid hit-and-run
#'   direction.
#' @param bound_truncation Truncation applied to infinite line-segment
#'   bounds inside hit-and-run.
#' @param n_threads OpenMP threads for the per-observation dual updates.
#'   Results do not depend on this value (each observation has its own RNG
#'   stream).
#'
#' @return An object of class `crr_control`.
#'
#' @export
crr_control <- function(n_iter_hit_and_run = 100,
                        rho = 1,
                        zeta_block = 100,
                        max_dir_tries = 1000,
                        bound_truncation = 1e5,
                        n_threads = 1) {
  stopifnot(
    n_iter_hit_and_run >= 2,
    rho > 0,
    zeta_block >= 1,
    max_dir_tries >= 1,
    bound_truncation > 0,
    n_threads >= 1
  )
  structure(
    list(
      n_iter_hit_and_run = as.integer(n_iter_hit_and_run),
      rho = rho,
      zeta_block = as.integer(zeta_block),
      max_dir_tries = as.integer(max_dir_tries),
      bound_truncation = bound_truncation,
      n_threads = as.integer(n_threads)
    ),
    class = "crr_control"
  )
}

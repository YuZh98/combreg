#' Sampler control parameters
#'
#' Tuning parameters for the MH-within-Gibbs sampler in [crr()].
#'
#' @param n_iter_hit_and_run Inner hit-and-run steps per dual update. The
#'   default (100) matches the simulation study of the paper; smaller values
#'   are faster but mix the dual variable less well.
#' @param rho Rate of the exponential dual kernel (ignored for
#'   `kernel = "half_gaussian"`).
#' @param zeta_block Number of response coordinates updated per MH step for the
#'   latent utilities. Either a positive integer for a fixed block size (the
#'   paper uses `min(d, 100)`), or `"adaptive"` (the default) to tune the block
#'   size automatically during warmup to hit `zeta_target_accept`. The optimal
#'   block size depends on the constraint geometry, not just `d`, so the
#'   adaptive rule is more robust across problems; pass an integer to reproduce
#'   a fixed-block run.
#' @param zeta_target_accept Target MH acceptance rate for the adaptive block
#'   controller (ignored when `zeta_block` is an integer). The default `0.6`
#'   maximizes mixing efficiency across a wide range of problems.
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
                        zeta_block = "adaptive",
                        zeta_target_accept = 0.6,
                        max_dir_tries = 1000,
                        bound_truncation = 1e5,
                        n_threads = 1) {
  adaptive <- identical(zeta_block, "adaptive")
  if (!adaptive && !(is.numeric(zeta_block) && length(zeta_block) == 1L &&
                     is.finite(zeta_block) && zeta_block >= 1)) {
    stop("zeta_block must be a positive integer or \"adaptive\"", call. = FALSE)
  }
  stopifnot(
    n_iter_hit_and_run >= 2,
    rho > 0,
    zeta_target_accept > 0, zeta_target_accept < 1,
    max_dir_tries >= 1,
    bound_truncation > 0,
    n_threads >= 1
  )
  structure(
    list(
      n_iter_hit_and_run = as.integer(n_iter_hit_and_run),
      rho = rho,
      zeta_block = if (adaptive) "adaptive" else as.integer(zeta_block),
      zeta_target_accept = zeta_target_accept,
      max_dir_tries = as.integer(max_dir_tries),
      bound_truncation = bound_truncation,
      n_threads = as.integer(n_threads)
    ),
    class = "crr_control"
  )
}

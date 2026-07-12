#' Constraint system for combinatorial responses
#'
#' Defines the feasible set \eqn{\{y \in \{0,1\}^d : Ay \le b\}} of a
#' combinatorial response. The matrix `A` must be integer-valued and, for the
#' exactness of the dual-certificate augmentation used by [crr()], totally
#' unimodular.
#'
#' @param A Integer constraint matrix (`m` x `d`).
#' @param b Integer right-hand side vector of length `m`.
#' @param direction Either a single string or a length-`m` character vector
#'   with entries `"<="` or `">="`. Internally all constraints are stored in
#'   `"<="` form.
#' @param check_tum Verify total unimodularity of `A` (exhaustive check,
#'   exponential in `min(m, d)`). Set to `FALSE` for large matrices whose TUM
#'   structure is known, e.g. network/incidence matrices.
#'
#' @return An object of class `crr_constraints` with elements `A`, `b`
#'   (in `<=` form), `m`, `d`.
#'
#' @examples
#' A <- rbind(c(1, 1, 0), c(0, 1, 1))
#' con <- crr_constraints(A, b = c(1, 1))
#'
#' @export
crr_constraints <- function(A, b, direction = "<=", check_tum = TRUE) {
  A <- as.matrix(A)
  if (!all(A == round(A))) {
    stop("A must be integer-valued", call. = FALSE)
  }
  m <- nrow(A)
  d <- ncol(A)

  b <- as.numeric(b)
  if (length(b) != m) {
    stop("length(b) must equal nrow(A)", call. = FALSE)
  }
  if (!all(b == round(b))) {
    stop("b must be integer-valued", call. = FALSE)
  }

  if (length(direction) == 1) direction <- rep(direction, m)
  if (length(direction) != m || !all(direction %in% c("<=", ">="))) {
    stop("direction must be \"<=\" or \">=\", length 1 or nrow(A)", call. = FALSE)
  }

  flip <- direction == ">="
  A[flip, ] <- -A[flip, , drop = FALSE]
  b[flip] <- -b[flip]

  if (any(rowSums(abs(A)) == 0)) {
    stop("A contains all-zero rows", call. = FALSE)
  }

  if (check_tum) {
    res <- check_tum_cpp(matrix(as.integer(A), m, d))
    if (!res$isTUM) {
      stop(
        "A is not totally unimodular (submatrix rows ",
        paste(res$rows, collapse = ","), ", cols ",
        paste(res$cols, collapse = ","), " has determinant ",
        res$determinant, "). The dual-certificate augmentation requires a ",
        "totally unimodular constraint matrix. If TUM is known by ",
        "construction, use check_tum = FALSE.",
        call. = FALSE
      )
    }
  }

  structure(
    list(A = A, b = b, m = m, d = d),
    class = "crr_constraints"
  )
}

#' @export
print.crr_constraints <- function(x, ...) {
  cat("<crr_constraints>: ", x$m, " constraints on {0,1}^", x$d,
      " (stored as A y <= b)\n", sep = "")
  invisible(x)
}

#' Check total unimodularity
#'
#' Exhaustively checks whether every square submatrix of `A` has determinant
#' in \{-1, 0, 1\}. Cost grows exponentially in `min(nrow(A), ncol(A))`; use
#' only for small-to-moderate matrices.
#'
#' @param A Integer matrix.
#'
#' @return `TRUE` or `FALSE`, with attribute `"witness"` (a list with `rows`,
#'   `cols`, `determinant`) when `FALSE`.
#'
#' @examples
#' is_tum(rbind(c(1, 1, 0), c(0, 1, 1)))       # TRUE
#' is_tum(rbind(c(1, 1, 0), c(1, 0, 1), c(0, 1, 1)))  # FALSE
#'
#' @export
is_tum <- function(A) {
  A <- as.matrix(A)
  if (!all(A == round(A))) stop("A must be integer-valued", call. = FALSE)
  res <- check_tum_cpp(matrix(as.integer(A), nrow(A), ncol(A)))
  out <- res$isTUM
  if (!out) {
    attr(out, "witness") <- list(
      rows = res$rows, cols = res$cols, determinant = res$determinant
    )
  }
  out
}

#' Check feasibility of responses
#'
#' @param constraints A [crr_constraints] object.
#' @param Y Response matrix (`n` x `d`), rows are responses.
#'
#' @return Logical vector of length `n`: does each row satisfy `A y <= b` and
#'   `y` in \{0,1\}?
#'
#' @export
is_feasible <- function(constraints, Y) {
  stopifnot(inherits(constraints, "crr_constraints"))
  Y <- as.matrix(Y)
  if (ncol(Y) != constraints$d) {
    stop("ncol(Y) must equal constraints$d", call. = FALSE)
  }
  binary <- apply(Y == 0 | Y == 1, 1, all)
  lhs <- Y %*% t(constraints$A)
  satisfied <- apply(lhs <= rep(constraints$b, each = nrow(Y)) + 1e-9, 1, all)
  binary & satisfied
}

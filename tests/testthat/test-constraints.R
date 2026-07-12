test_that("crr_constraints validates and normalizes inputs", {
  A <- rbind(c(1, 1, 0), c(0, 1, 1))
  con <- crr_constraints(A, b = c(1, 1))
  expect_s3_class(con, "crr_constraints")
  expect_equal(con$m, 2)
  expect_equal(con$d, 3)

  con_ge <- crr_constraints(A, b = c(1, 1), direction = ">=")
  expect_equal(con_ge$A, -A, ignore_attr = TRUE)
  expect_equal(con_ge$b, c(-1, -1))

  expect_error(crr_constraints(A, b = 1), "length\\(b\\)")
  expect_error(crr_constraints(A + 0.5, b = c(1, 1)), "integer-valued")
  expect_error(crr_constraints(rbind(A, 0), b = c(1, 1, 1)), "all-zero")
})

test_that("is_tum detects TUM and non-TUM matrices", {
  # Interval matrix: TUM.
  expect_true(is_tum(rbind(c(1, 1, 0), c(0, 1, 1))))
  # Odd-cycle incidence matrix: det = 2, not TUM.
  bad <- rbind(c(1, 1, 0), c(1, 0, 1), c(0, 1, 1))
  res <- is_tum(bad)
  expect_false(res)
  expect_equal(abs(attr(res, "witness")$determinant), 2)
  expect_error(crr_constraints(bad, b = rep(1, 3)), "not totally unimodular")
})

test_that("is_feasible checks binary support and A y <= b", {
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  Y <- rbind(
    c(1, 0, 1),  # feasible
    c(1, 1, 0),  # violates first constraint
    c(0, 2, 0)   # not binary
  )
  expect_equal(is_feasible(con, Y), c(TRUE, FALSE, FALSE))
})

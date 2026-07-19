test_that("exact duplicate rows are removed by default", {
  A <- rbind(c(1, 1, 0), c(1, 1, 0), c(0, 1, 1))   # row 1 == row 2
  con <- crr_constraints(A, b = c(1, 1, 1), check_tum = FALSE)
  expect_equal(con$m, 2L)
})

test_that("a repeated A row keeps the tighter (smaller) b", {
  A <- rbind(c(1, 1, 0), c(1, 1, 0))
  con <- crr_constraints(A, b = c(2, 1), check_tum = FALSE)
  expect_equal(con$m, 1L)
  expect_equal(con$b, 1)
})

test_that("dedup leaves the feasible set unchanged", {
  A <- rbind(c(1, 1, 0), c(1, 1, 0), c(0, 1, 1))
  b <- c(2, 1, 1)
  raw <- crr_constraints(A, b, check_tum = FALSE, dedup = FALSE)
  ded <- crr_constraints(A, b, check_tum = FALSE, dedup = TRUE)
  Y <- as.matrix(expand.grid(rep(list(0:1), 3)))   # all 8 binary vectors
  expect_identical(is_feasible(raw, Y), is_feasible(ded, Y))
})

test_that("dedup preserves total unimodularity (check passes)", {
  A <- rbind(c(1, 1, 0), c(1, 1, 0), c(0, 1, 1))
  ded <- crr_constraints(A, b = c(1, 1, 1), check_tum = TRUE)  # errors if not TU
  expect_true(is_tum(ded$A))
})

test_that("dedup = FALSE keeps every row", {
  A <- rbind(c(1, 1, 0), c(1, 1, 0))
  con <- crr_constraints(A, b = c(1, 1), check_tum = FALSE, dedup = FALSE)
  expect_equal(con$m, 2L)
})

test_that("random_constraints preserves the requested m (faithful, no dedup)", {
  for (s in 1:30) {
    set.seed(s)
    con <- random_constraints(20, 15)  # ~27% of seeds yield a duplicate row
    expect_equal(con$m, 15L)
  }
})

testthat::test_that("deterministic_schedule: handles small sample sizes across many arms", {
  set.seed(1)
  sch <- deterministic_schedule(
    sample_size = 2,
    arms = c("A", "B", "C"),
    allocation = c(1, 1, 1),
    enrollment = list(end_time = 1, rate = 3),
    dropout = list(end_time = 1, rate = 0)
  )

  testthat::expect_true(is.data.frame(sch))
  testthat::expect_equal(sum(sch$enroll), 2L)
  testthat::expect_true(all(sch$enroll >= 0L))
  testthat::expect_true(all(sch$drop >= 0L))
})

testthat::test_that("deterministic_schedule: keeps per-arm enrollment at target", {
  arms <- c("A", "B", "C")
  allocation <- c(3, 2, 1)
  sample_size <- 12

  sch <- deterministic_schedule(
    sample_size = sample_size,
    arms = arms,
    allocation = allocation,
    enrollment = list(end_time = c(2, 4), rate = c(4, 3)),
    dropout = list(end_time = c(2, 4), rate = c(0, 1))
  )

  expected <- as.integer(allocation / sum(allocation) * sample_size)
  names(expected) <- arms
  observed <- vapply(arms, function(a) sum(sch$enroll[sch$arm == a]), integer(1))
  testthat::expect_equal(observed, expected)
  testthat::expect_equal(sum(observed), sample_size)
})

testthat::test_that("stochastic_schedule: enrollment totals match sample_size", {
  set.seed(99)
  sch <- stochastic_schedule(
    sample_size = 25,
    arms = c("A", "B"),
    allocation = c(1, 2),
    enrollment = function(n) rexp(n, 1),
    dropout = function(n) rexp(n, 0.5)
  )

  testthat::expect_equal(sum(sch$enroll), 25L)
  testthat::expect_equal(sum(sch$drop), 25L)
})

testthat::test_that("deterministic_schedule: correction row before first dropout boundary", {
  # findInterval(next_t, end_time) is 0 when the correction lands before the
  # first period boundary; rate[0] used to crash df_add with a length-0 column.
  sch <- deterministic_schedule(
    sample_size = 7,
    arms = c("x", "y"),
    allocation = c(1, 3),
    enrollment = list(end_time = 3, rate = 10),
    dropout = list(end_time = 3, rate = 2)
  )
  testthat::expect_equal(sum(sch$enroll), 7L)
  testthat::expect_s3_class(sch, "data.frame")
})

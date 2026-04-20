# Timer: initialize

testthat::test_that("Timer initialize: creates instance with correct defaults", {
  t <- Timer$new(name = "test_timer")

  testthat::expect_equal(t$name, "test_timer")
  testthat::expect_true(is.list(t$timelist))
  testthat::expect_equal(length(t$timelist), 0L)
  testthat::expect_r6_class(t, "Timer")
})

testthat::test_that("Timer initialize: accepts custom timelist", {
  tp <- list(list(time = 1, arm = "A", dropper = 1L, enroller = 5L))

  t <- Timer$new(name = "t", timelist = tp)

  testthat::expect_equal(length(t$timelist), 1L)
})

testthat::test_that("Timer initialize: errors on non-character name", {
  testthat::expect_error(Timer$new(name = 123))
  testthat::expect_error(Timer$new(name = TRUE))
  testthat::expect_error(Timer$new(name = NULL))
})

# Timer: add_timepoint

testthat::test_that("add_timepoint: appends a single timepoint", {
  t <- Timer$new(name = "t")
  t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 10L)

  testthat::expect_equal(length(t$timelist), 1L)
  testthat::expect_equal(t$timelist[[1]]$time, 1)
  testthat::expect_equal(t$timelist[[1]]$arm, "A")
  testthat::expect_equal(t$timelist[[1]]$dropper, 2L)
  testthat::expect_equal(t$timelist[[1]]$enroller, 10L)
})

testthat::test_that("add_timepoint: appends multiple timepoints preserving order", {
  t <- Timer$new(name = "t")
  t$add_timepoint(time = 1, arm = "A", dropper = 1L, enroller = 5L)
  t$add_timepoint(time = 2, arm = "B", dropper = 2L, enroller = 8L)
  t$add_timepoint(time = 3, arm = "A", dropper = 0L, enroller = 3L)

  testthat::expect_equal(length(t$timelist), 3L)
  testthat::expect_equal(t$timelist[[2]]$arm, "B")
  testthat::expect_equal(t$timelist[[3]]$time, 3)
})

testthat::test_that("add_timepoint: errors when dropper is not integer", {
  t <- Timer$new(name = "t")
  testthat::expect_error(
    t$add_timepoint(time = 1, arm = "A", dropper = 2, enroller = 10L)
  )
})

testthat::test_that("add_timepoint: errors when enroller is not integer", {
  t <- Timer$new(name = "t")
  testthat::expect_error(
    t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 10)
  )
})

# Timer: get_end_timepoint

testthat::test_that("get_end_timepoint: returns max time across all arms", {
  t <- Timer$new(name = "t")
  t$add_timepoint(time = 1, arm = "A", dropper = 1L, enroller = 5L)
  t$add_timepoint(time = 5, arm = "B", dropper = 1L, enroller = 5L)
  t$add_timepoint(time = 3, arm = "A", dropper = 0L, enroller = 2L)

  testthat::expect_equal(t$get_end_timepoint(), 5)
})

testthat::test_that("get_end_timepoint: works with a single timepoint", {
  t <- Timer$new(name = "t")
  t$add_timepoint(time = 3.5, arm = "A", dropper = 7L, enroller = 22L)

  testthat::expect_equal(t$get_end_timepoint(), 3.5)
})

testthat::test_that("get_end_timepoint: errors on empty timelist", {
  t <- Timer$new(name = "t")
  testthat::expect_error(t$get_end_timepoint())
})

# Timer: get_n_arms

testthat::test_that("get_n_arms: counts unique arms", {
  t <- Timer$new(name = "t")
  t$add_timepoint(time = 1, arm = "A", dropper = 1L, enroller = 5L)
  t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 5L)
  t$add_timepoint(time = 1, arm = "B", dropper = 1L, enroller = 5L)
  t$add_timepoint(time = 1, arm = "C", dropper = 1L, enroller = 5L)

  testthat::expect_equal(t$get_n_arms(), 3L)
})

# Timer: get_unique_times

testthat::test_that("get_unique_times: returns sorted unique time values", {
  t <- Timer$new(name = "t")
  t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 5L)
  t$add_timepoint(time = 1, arm = "B", dropper = 1L, enroller = 5L)
  t$add_timepoint(time = 2, arm = "B", dropper = 1L, enroller = 5L)

  result <- t$get_unique_times()
  testthat::expect_equal(sort(result), c(1, 2))
})

testthat::test_that("get_unique_times: returns empty for empty timelist", {
  t <- Timer$new(name = "t")
  result <- t$get_unique_times()
  testthat::expect_equal(length(result), 0L)
})

# Timer: get_timepoint

testthat::test_that("get_timepoint: returns matching timepoint", {
  t <- Timer$new(name = "t")
  t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 10L)
  t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 12L)
  t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 8L)

  tp <- t$get_timepoint("A", 1)
  testthat::expect_true(is.list(tp))
  testthat::expect_equal(tp$time, 1)
  testthat::expect_equal(tp$arm, "A")
  testthat::expect_equal(tp$dropper, 2L)
  testthat::expect_equal(tp$enroller, 10L)
})

testthat::test_that("get_timepoint: returns NULL when no match", {
  t <- Timer$new(name = "t")
  t$add_timepoint(time = 1, arm = "A", dropper = 1L, enroller = 5L)

  testthat::expect_null(t$get_timepoint("B", 1))
  testthat::expect_null(t$get_timepoint("A", 9))
})

testthat::test_that("get_timepoint: errors on multiple matches", {
  t <- Timer$new(name = "t")
  t$add_timepoint(time = 1, arm = "A", dropper = 1L, enroller = 5L)
  t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 6L)

  testthat::expect_error(t$get_timepoint("A", 1), "Multiple timepoints")
})

testthat::test_that("get_timepoint: errors when arm is missing", {
  t <- Timer$new(name = "t")
  testthat::expect_error(t$get_timepoint(i = 1), "`arm` is required")
})

testthat::test_that("get_timepoint: errors when i is missing", {
  t <- Timer$new(name = "t")
  testthat::expect_error(t$get_timepoint(arm = "A"), "`i` is required")
})


# Timer: initialize

testthat::test_that("Timer initialize: creates instance with correct defaults", {
  t <- Timer$new(name = "test_timer")

  testthat::expect_equal(t$name, "test_timer")
  testthat::expect_true(is.data.frame(t$timelist))
  testthat::expect_equal(nrow(t$timelist), 0L)
  testthat::expect_r6_class(t, "Timer")
})

testthat::test_that("Timer initialize: accepts custom timelist", {
  tp <- data.frame(
    time = 1,
    arm = "A",
    drop = 1L,
    enroll = 5L,
    stringsAsFactors = FALSE
  )

  t <- Timer$new(name = "t", timelist = tp)

  testthat::expect_true(is.data.frame(t$timelist))
  testthat::expect_equal(nrow(t$timelist), 1L)
  testthat::expect_equal(t$timelist$arm, "A")
})

testthat::test_that("Timer initialize: errors on non-character name", {
  testthat::expect_error(Timer$new(name = 123))
  testthat::expect_error(Timer$new(name = TRUE))
  testthat::expect_error(Timer$new(name = NULL))
})

# Timer: add_schedule

testthat::test_that("add_schedule: appends a single timepoint", {
  t <- Timer$new(name = "t")
  t$add_schedule(data.frame(time = 1, arm = "A", drop = 2L, enroll = 10L))

  testthat::expect_equal(nrow(t$timelist), 1L)
  testthat::expect_equal(t$timelist$time[[1]], 1)
  testthat::expect_equal(t$timelist$arm[[1]], "A")
  testthat::expect_equal(t$timelist$drop[[1]], 2L)
  testthat::expect_equal(t$timelist$enroll[[1]], 10L)
})

testthat::test_that("add_schedule: appends multiple timepoints preserving order", {
  t <- Timer$new(name = "t")
  t$add_schedule(data.frame(time = 1, arm = "A", drop = 1L, enroll = 5L))
  t$add_schedule(data.frame(time = 2, arm = "B", drop = 2L, enroll = 8L))
  t$add_schedule(data.frame(time = 3, arm = "A", drop = 0L, enroll = 3L))

  testthat::expect_equal(nrow(t$timelist), 3L)
  testthat::expect_equal(t$timelist$arm[[2]], "B")
  testthat::expect_equal(t$timelist$time[[3]], 3)
})

testthat::test_that("add_schedule: errors when drop is not integer", {
  t <- Timer$new(name = "t")
  testthat::expect_error(
    t$add_schedule(data.frame(time = 1, arm = "A", drop = 2, enroll = 10L)),
    "`drop` must be an integer vector"
  )
})

testthat::test_that("add_schedule: errors when enroll is not integer", {
  t <- Timer$new(name = "t")
  testthat::expect_error(
    t$add_schedule(data.frame(time = 1, arm = "A", drop = 2L, enroll = 10)),
    "`enroll` must be an integer vector"
  )
})

testthat::test_that("add_schedule: errors when schedule is not a data.frame", {
  t <- Timer$new(name = "t")
  testthat::expect_error(t$add_schedule("not_a_df"), "must be a data.frame")
  testthat::expect_error(t$add_schedule(1:5), "must be a data.frame")
  testthat::expect_error(t$add_schedule(), "must be a data.frame")
})

testthat::test_that("add_schedule: errors when required columns are missing", {
  t <- Timer$new(name = "t")
  testthat::expect_error(
    t$add_schedule(data.frame(time = 1, arm = "A")),
    "Missing required columns in schedule: enroll, drop"
  )
})

testthat::test_that("add_schedule: adds a whole multi-row schedule and returns timer invisibly", {
  t <- Timer$new(name = "t")
  df <- data.frame(time = c(1, 2, 3), arm = "A", enroll = c(5L, 3L, 2L), drop = c(0L, 1L, 1L))
  result <- t$add_schedule(df)

  testthat::expect_true(inherits(result, "Timer"))
  testthat::expect_equal(nrow(t$timelist), 3L)
  testthat::expect_equal(t$timelist$enroll, c(5L, 3L, 2L))
})

testthat::test_that("add_schedule: normalizes tibble input to a plain data.frame", {
  t <- Timer$new(name = "t")
  t$add_schedule(dplyr::group_by(
    data.frame(time = 1, arm = "A", enroll = 3L, drop = 0L), arm
  ))
  testthat::expect_identical(class(t$timelist), "data.frame")
})

# Timer: get_end_timepoint

testthat::test_that("get_end_timepoint: returns max time across all arms", {
  t <- Timer$new(name = "t")
  t$add_schedule(data.frame(time = 1, arm = "A", drop = 1L, enroll = 5L))
  t$add_schedule(data.frame(time = 5, arm = "B", drop = 1L, enroll = 5L))
  t$add_schedule(data.frame(time = 3, arm = "A", drop = 0L, enroll = 2L))

  testthat::expect_equal(t$get_end_timepoint(), 5)
})

testthat::test_that("get_end_timepoint: works with a single timepoint", {
  t <- Timer$new(name = "t")
  t$add_schedule(data.frame(time = 3.5, arm = "A", drop = 7L, enroll = 22L))

  testthat::expect_equal(t$get_end_timepoint(), 3.5)
})

testthat::test_that("get_end_timepoint: errors on empty timelist", {
  t <- Timer$new(name = "t")
  testthat::expect_error(t$get_end_timepoint(), "`timelist` is empty")
})

# Timer: get_n_arms

testthat::test_that("get_n_arms: counts unique arms", {
  t <- Timer$new(name = "t")
  t$add_schedule(data.frame(time = 1, arm = "A", drop = 1L, enroll = 5L))
  t$add_schedule(data.frame(time = 2, arm = "A", drop = 1L, enroll = 5L))
  t$add_schedule(data.frame(time = 1, arm = "B", drop = 1L, enroll = 5L))
  t$add_schedule(data.frame(time = 1, arm = "C", drop = 1L, enroll = 5L))

  testthat::expect_equal(t$get_n_arms(), 3L)
})

# Timer: get_unique_times

testthat::test_that("get_unique_times: returns sorted unique time values", {
  t <- Timer$new(name = "t")
  t$add_schedule(data.frame(time = 2, arm = "A", drop = 1L, enroll = 5L))
  t$add_schedule(data.frame(time = 1, arm = "B", drop = 1L, enroll = 5L))
  t$add_schedule(data.frame(time = 2, arm = "B", drop = 1L, enroll = 5L))

  result <- t$get_unique_times()
  testthat::expect_equal(sort(result), c(1, 2))
})

testthat::test_that("get_unique_times: returns empty for empty timelist", {
  t <- Timer$new(name = "t")
  result <- t$get_unique_times()
  testthat::expect_equal(length(result), 0L)
})


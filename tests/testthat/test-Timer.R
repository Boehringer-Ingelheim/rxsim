# Timer: initialize

testthat::test_that("Timer initialize: creates instance with correct defaults", {
  t <- Timer$new(name = "test_timer")

  testthat::expect_equal(t$name, "test_timer")
  testthat::expect_true(is.list(t$timelist))
  testthat::expect_equal(length(t$timelist), 0L)
  testthat::expect_true(is.list(t$conditions))
  testthat::expect_equal(length(t$conditions), 0L)
  testthat::expect_r6_class(t, "Timer")
})

testthat::test_that("Timer initialize: accepts custom timelist and conditions", {
  tp <- list(list(time = 1, arm = "A", dropper = 1L, enroller = 5L))
  cond <- list(list(where = NULL, analysis = NULL, name = "cond1"))

  t <- Timer$new(name = "t", timelist = tp, conditions = cond)

  testthat::expect_equal(length(t$timelist), 1L)
  testthat::expect_equal(length(t$conditions), 1L)
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

# Timer: add_condition

testthat::test_that("add_condition: stores condition with quosures", {
  t <- Timer$new(name = "t")
  t$add_condition(status == "active", name = "cond1")

  testthat::expect_equal(length(t$conditions), 1L)
  testthat::expect_equal(t$conditions[[1]]$name, "cond1")
  testthat::expect_true(is.list(t$conditions[[1]]$where))
})

testthat::test_that("add_condition: stores analysis function", {
  t <- Timer$new(name = "t")
  fn <- function(dat, current_time) nrow(dat)

  t$add_condition(arm == "A", analysis = fn, name = "with_fn")

  testthat::expect_true(is.function(t$conditions[[1]]$analysis))
  testthat::expect_null(t$conditions[[1]]$analysis(NULL, 1))
})

testthat::test_that("add_condition: stores multiple where predicates", {
  t <- Timer$new(name = "t")
  t$add_condition(arm == "A", status == "active", name = "multi_pred")

  testthat::expect_equal(length(t$conditions[[1]]$where), 2L)
})

# Timer: check_conditions

testthat::test_that("check_conditions: applies filter and returns filtered data with warning", {
  t <- Timer$new(name = "t")
  t$add_condition(arm == "A", name = "armA_only")

  df <- data.frame(
    id = 1:6,
    arm = c("A", "A", "B", "B", "A", "B"),
    stringsAsFactors = FALSE
  )

  testthat::expect_warning(
    res <- t$check_conditions(locked_data = df, current_time = 1),
    "returning filtered data"
  )
  testthat::expect_true(is.data.frame(res$armA_only))
  testthat::expect_equal(nrow(res$armA_only), 3L)
  testthat::expect_true(all(res$armA_only$arm == "A"))
})

testthat::test_that("check_conditions: runs analysis function when provided", {
  t <- Timer$new(name = "t")
  my_analysis <- function(dat, current_time) {
    data.frame(n = nrow(dat), ct = current_time)
  }

  t$add_condition(
    status == "active",
    analysis = my_analysis,
    name = "active_count"
  )

  df <- data.frame(
    id = 1:5,
    status = c("active", "inactive", "active", "active", "inactive"),
    stringsAsFactors = FALSE
  )

  res <- t$check_conditions(locked_data = df, current_time = 3)

  testthat::expect_equal(res$active_count$n, 3L)
  testthat::expect_equal(res$active_count$ct, 3)
})

testthat::test_that("check_conditions: skips conditions with empty filtered data", {
  t <- Timer$new(name = "t")
  t$add_condition(arm == "Z", name = "empty_cond")

  df <- data.frame(
    id = 1:3,
    arm = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )

  res <- t$check_conditions(locked_data = df, current_time = 1)
  testthat::expect_equal(length(res), 0L)
})

testthat::test_that("check_conditions: uses integer index as key when name is NULL", {
  t <- Timer$new(name = "t")
  t$add_condition(id > 0)

  df <- data.frame(id = 1:3, stringsAsFactors = FALSE)

  testthat::expect_warning(
    res <- t$check_conditions(locked_data = df, current_time = 1),
    "returning filtered data"
  )

  testthat::expect_equal(length(res), 1L)
  testthat::expect_equal(nrow(res[[1]]), 3L)
})

testthat::test_that("check_conditions: handles multiple conditions", {
  t <- Timer$new(name = "t")

  analysis_a <- function(dat, ct) data.frame(n_a = nrow(dat))
  analysis_b <- function(dat, ct) data.frame(n_b = nrow(dat))

  t$add_condition(arm == "A", analysis = analysis_a, name = "cond_a")
  t$add_condition(arm == "B", analysis = analysis_b, name = "cond_b")

  df <- data.frame(
    id = 1:4,
    arm = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  res <- t$check_conditions(locked_data = df, current_time = 1)

  testthat::expect_equal(res$cond_a$n_a, 2L)
  testthat::expect_equal(res$cond_b$n_b, 2L)
})

testthat::test_that("check_conditions: condition with no where returns full data", {
  t <- Timer$new(name = "t")
  t$add_condition(analysis = function(dat, ct) nrow(dat), name = "all_data")

  df <- data.frame(id = 1:5, stringsAsFactors = FALSE)
  res <- t$check_conditions(locked_data = df, current_time = 1)

  testthat::expect_equal(res$all_data, 5L)
})

testthat::test_that("check_conditions: condition with multiple where predicates filters correctly", {
  t <- Timer$new(name = "t")
  t$add_condition(
    arm == "A", status == "active",
    analysis = function(dat, ct) nrow(dat),
    name = "multi_filter"
  )

  df <- data.frame(
    id = 1:6,
    arm = c("A", "A", "B", "A", "B", "A"),
    status = c("active", "inactive", "active", "active", "active", "active"),
    stringsAsFactors = FALSE
  )

  res <- t$check_conditions(locked_data = df, current_time = 1)
  testthat::expect_equal(res$multi_filter, 3L)
})

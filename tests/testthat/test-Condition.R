# Condition: initialize

testthat::test_that("Condition initialize: creates instance with correct defaults", {
  cond <- Condition$new(where = NULL)

  testthat::expect_null(cond$where)
  testthat::expect_null(cond$analysis)
  testthat::expect_null(cond$name)
  testthat::expect_equal(cond$cooldown, 0)
  testthat::expect_equal(cond$max_triggers, 1L)
  testthat::expect_equal(cond$trigger_count, 0L)
  testthat::expect_true(is.na(cond$last_trigger_time))
  testthat::expect_r6_class(cond, "Condition")
})

testthat::test_that("Condition initialize: stores all arguments", {
  fn <- function(dat, ct) nrow(dat)
  cond <- Condition$new(
    where        = rlang::quos(arm == "A"),
    analysis     = fn,
    name         = "my_cond",
    cooldown     = 2,
    max_triggers = 3L
  )

  testthat::expect_equal(length(cond$where), 1L)
  testthat::expect_true(is.function(cond$analysis))
  testthat::expect_equal(cond$name, "my_cond")
  testthat::expect_equal(cond$cooldown, 2)
  testthat::expect_equal(cond$max_triggers, 3L)
})

testthat::test_that("Condition initialize: errors on invalid cooldown", {
  testthat::expect_error(Condition$new(cooldown = -1), "`cooldown`")
  testthat::expect_error(Condition$new(cooldown = NA_real_), "`cooldown`")
  testthat::expect_error(Condition$new(cooldown = c(1, 2)), "`cooldown`")
})

testthat::test_that("Condition initialize: errors on invalid max_triggers", {
  testthat::expect_error(Condition$new(max_triggers = -1L), "`max_triggers`")
  testthat::expect_error(Condition$new(max_triggers = NA_integer_), "`max_triggers`")
  testthat::expect_error(Condition$new(max_triggers = c(1L, 2L)), "`max_triggers`")
})

# Condition: check_conditions — filtering and analysis

testthat::test_that("check_conditions: applies filter and returns filtered data with warning", {
  cond <- Condition$new(
    where = rlang::quos(arm == "A"),
    name  = "armA_only"
  )

  df <- data.frame(
    id  = 1:6,
    arm = c("A", "A", "B", "B", "A", "B"),
    stringsAsFactors = FALSE
  )

  testthat::expect_warning(
    res <- cond$check_conditions(locked_data = df, current_time = 1),
    "returning filtered data"
  )
  testthat::expect_true(is.data.frame(res$armA_only))
  testthat::expect_equal(nrow(res$armA_only), 3L)
  testthat::expect_true(all(res$armA_only$arm == "A"))
})

testthat::test_that("check_conditions: runs analysis function when provided", {
  my_analysis <- function(dat, current_time) {
    data.frame(n = nrow(dat), ct = current_time)
  }
  cond <- Condition$new(
    where    = rlang::quos(status == "active"),
    analysis = my_analysis,
    name     = "active_count",
    max_triggers = 10L
  )

  df <- data.frame(
    id     = 1:5,
    status = c("active", "inactive", "active", "active", "inactive"),
    stringsAsFactors = FALSE
  )

  res <- cond$check_conditions(locked_data = df, current_time = 3)
  testthat::expect_equal(res$active_count$n, 3L)
  testthat::expect_equal(res$active_count$ct, 3)
})

testthat::test_that("check_conditions: skips when filter result is empty", {
  cond <- Condition$new(
    where = rlang::quos(arm == "Z"),
    name  = "empty_cond"
  )

  df <- data.frame(
    id  = 1:3,
    arm = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )

  res <- cond$check_conditions(locked_data = df, current_time = 1)
  testthat::expect_equal(length(res), 0L)
})

testthat::test_that("check_conditions: uses integer 1L as key when name is NULL", {
  cond <- Condition$new(where = rlang::quos(id > 0))

  df <- data.frame(id = 1:3, stringsAsFactors = FALSE)

  testthat::expect_warning(
    res <- cond$check_conditions(locked_data = df, current_time = 1),
    "returning filtered data"
  )
  testthat::expect_equal(length(res), 1L)
  testthat::expect_equal(nrow(res[[1]]), 3L)
})

testthat::test_that("check_conditions: condition with no where returns full data", {
  cond <- Condition$new(
    analysis = function(dat, ct) nrow(dat),
    name     = "all_data"
  )

  df <- data.frame(id = 1:5, stringsAsFactors = FALSE)
  res <- cond$check_conditions(locked_data = df, current_time = 1)
  testthat::expect_equal(res$all_data, 5L)
})

testthat::test_that("check_conditions: multiple where predicates are ANDed", {
  cond <- Condition$new(
    where    = rlang::quos(arm == "A", status == "active"),
    analysis = function(dat, ct) nrow(dat),
    name     = "multi_filter",
    max_triggers = 10L
  )

  df <- data.frame(
    id     = 1:6,
    arm    = c("A", "A", "B", "A", "B", "A"),
    status = c("active", "inactive", "active", "active", "active", "active"),
    stringsAsFactors = FALSE
  )

  res <- cond$check_conditions(locked_data = df, current_time = 1)
  testthat::expect_equal(res$multi_filter, 3L)
})

testthat::test_that("check_conditions: errors when locked_data is not a data.frame", {
  cond <- Condition$new()
  testthat::expect_error(cond$check_conditions(locked_data = list(a = 1), current_time = 1))
  testthat::expect_error(cond$check_conditions(locked_data = "bad", current_time = 1))
})

# Condition: trigger guards

testthat::test_that("check_conditions: max_triggers=1 prevents second fire", {
  cond <- Condition$new(
    analysis     = function(dat, ct) data.frame(n = nrow(dat)),
    name         = "once_only",
    max_triggers = 1L
  )
  df <- data.frame(id = 1:3, stringsAsFactors = FALSE)

  res1 <- cond$check_conditions(locked_data = df, current_time = 1)
  testthat::expect_equal(length(res1), 1L)

  # Second call: should not fire
  res2 <- cond$check_conditions(locked_data = df, current_time = 2)
  testthat::expect_equal(length(res2), 0L)
})

testthat::test_that("check_conditions: max_triggers=Inf allows unlimited fires", {
  cond <- Condition$new(
    analysis     = function(dat, ct) data.frame(n = nrow(dat)),
    name         = "unlimited",
    max_triggers = Inf
  )
  df <- data.frame(id = 1:3, stringsAsFactors = FALSE)

  for (t in 1:5) {
    res <- cond$check_conditions(locked_data = df, current_time = t)
    testthat::expect_equal(length(res), 1L)
  }
  testthat::expect_equal(cond$trigger_count, 5L)
})

testthat::test_that("check_conditions: cooldown prevents re-trigger within period", {
  cond <- Condition$new(
    analysis     = function(dat, ct) data.frame(n = nrow(dat)),
    name         = "cooled",
    cooldown     = 5,
    max_triggers = 10L
  )
  df <- data.frame(id = 1:3, stringsAsFactors = FALSE)

  # Fires at t=1
  res1 <- cond$check_conditions(locked_data = df, current_time = 1)
  testthat::expect_equal(length(res1), 1L)

  # Does not fire at t=4 (only 3 units elapsed, need 5)
  res2 <- cond$check_conditions(locked_data = df, current_time = 4)
  testthat::expect_equal(length(res2), 0L)

  # Fires again at t=6 (5 units elapsed since t=1)
  res3 <- cond$check_conditions(locked_data = df, current_time = 6)
  testthat::expect_equal(length(res3), 1L)
})

testthat::test_that("check_conditions: state persists across calls", {
  cond <- Condition$new(
    analysis     = function(dat, ct) data.frame(n = nrow(dat)),
    name         = "stateful",
    max_triggers = 3L
  )
  df <- data.frame(id = 1:3, stringsAsFactors = FALSE)

  cond$check_conditions(locked_data = df, current_time = 1)
  cond$check_conditions(locked_data = df, current_time = 2)
  cond$check_conditions(locked_data = df, current_time = 3)

  testthat::expect_equal(cond$trigger_count, 3L)
  testthat::expect_equal(cond$last_trigger_time, 3)

  # 4th call should not fire
  res4 <- cond$check_conditions(locked_data = df, current_time = 4)
  testthat::expect_equal(length(res4), 0L)
  testthat::expect_equal(cond$trigger_count, 3L)
})

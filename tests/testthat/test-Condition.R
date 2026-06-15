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

# Condition: check_conditions  -  filtering and analysis

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

# --- analysis_args tests ---

testthat::test_that("analysis_args: extra named arg is passed to analysis function", {
  threshold <- 99  # wrong value  -  should be overridden by analysis_args
  cond <- Condition$new(
    analysis      = function(df, current_time, threshold) current_time > threshold,
    analysis_args = list(threshold = 5),
    name          = "threshold_test"
  )
  df  <- data.frame(id = 1:3)
  res <- cond$check_conditions(locked_data = df, current_time = 10)
  testthat::expect_true(res$threshold_test)
})

testthat::test_that("analysis_args: multiple extra args are all injected", {
  cond <- Condition$new(
    analysis = function(df, current_time, alpha, beta) {
      list(a = alpha, b = beta, n = nrow(df), t = current_time)
    },
    analysis_args = list(alpha = 0.05, beta = 0.8),
    name          = "multi_args"
  )
  df  <- data.frame(id = 1:4)
  res <- cond$check_conditions(locked_data = df, current_time = 7)$multi_args
  testthat::expect_equal(res$a, 0.05)
  testthat::expect_equal(res$b, 0.8)
  testthat::expect_equal(res$n, 4L)
  testthat::expect_equal(res$t, 7)
})

testthat::test_that("analysis_args: NULL analysis_args leaves call unchanged", {
  cond <- Condition$new(
    analysis      = function(df, current_time) nrow(df) + current_time,
    analysis_args = NULL,
    name          = "no_extra"
  )
  df  <- data.frame(id = 1:3)
  res <- cond$check_conditions(locked_data = df, current_time = 2)
  testthat::expect_equal(res$no_extra, 5)
})

testthat::test_that("analysis_args: scenario isolation  -  two conditions with different args", {
  make_cond <- function(cov_val) {
    Condition$new(
      analysis      = function(df, current_time, cov_matrix) cov_matrix,
      analysis_args = list(cov_matrix = cov_val),
      name          = "cov_test"
    )
  }
  cond1 <- make_cond(0.2)
  cond2 <- make_cond(0.9)
  df <- data.frame(id = 1:2)

  res1 <- cond1$check_conditions(locked_data = df, current_time = 1)$cov_test
  res2 <- cond2$check_conditions(locked_data = df, current_time = 1)$cov_test
  testthat::expect_equal(res1, 0.2)
  testthat::expect_equal(res2, 0.9)
})

testthat::test_that("analysis_args: works end-to-end through Condition$new", {
  trig <- value_trigger("id", ">=", 1)
  cond <- Condition$new(
    where         = trig,
    analysis      = function(df, current_time, scale) nrow(df) * scale,
    analysis_args = list(scale = 10),
    name          = "scaled"
  )
  df  <- data.frame(id = 1:3)
  res <- cond$check_conditions(locked_data = df, current_time = 1)
  testthat::expect_equal(res$scaled, 30)
})

testthat::test_that("analysis_args: works end-to-end through replicate_trial", {
  set.seed(1)
  ag <- list(final = list(
    trigger       = enroll_trigger(1.0, 6L),
    analysis      = function(df, current_time, multiplier) nrow(df) * multiplier,
    analysis_args = list(multiplier = 3L)
  ))
  pop_gen <- list(
    A = function(n) data.frame(id = seq_len(n), value = rnorm(n), readout_time = 0),
    B = function(n) data.frame(id = seq_len(n), value = rnorm(n), readout_time = 0)
  )
  trials <- replicate_trial(
    "t", 6L, c("A", "B"), c(1, 1),
    function(n) rexp(n, 1), function(n) rep(Inf, n),
    ag, pop_gen, 1L
  )
  testthat::expect_r6_class(trials[[1L]]$conditions[[1L]], "Condition")
  testthat::expect_equal(trials[[1L]]$conditions[[1L]]$analysis_args, list(multiplier = 3L))
})

testthat::test_that("collect_results: combines analyses with different columns", {
  sample_size   <- 100L
  arms          <- c("treatment", "placebo")
  allocation    <- c(1, 1)
  enrollment_fn <- function(n) rexp(n, 1)
  dropout_fn    <- function(n) rexp(n, 0.05)

  trt <- Population$new(
    name = "treatment",
    data = as_population_data(rnorm(sample_size / 2))
  )
  pbo <- Population$new(
    name = "placebo",
    data = as_population_data(rnorm(sample_size / 2, mean = 2))
  )

  t <- Timer$new(name = "timer")
  add_timepoints(t, stochastic_schedule(sample_size, arms, allocation, enrollment_fn, dropout_fn))

  final <- condition_enrollment_fraction(1.0, sample_size, analysis = function(df, time) {
    enrolled <- subset(df, !is.na(enroll_time))
    data.frame(
      n = nrow(enrolled),
      p_value = t.test(data ~ arm, data = enrolled)$p.value
    )
  })

  interim <- condition_enrollment_fraction(0.5, sample_size, analysis = function(df, time) {
    enrolled <- subset(df, !is.na(enroll_time))
    data.frame(
      n = nrow(enrolled),
      mean = mean(enrolled$data)
    )
  })

  trial <- Trial$new(name = "simple", timer = t, population = list(trt, pbo), conditions = list(final, interim))
  trial$run()

  res <- collect_results(list(trial))
  testthat::expect_equal(nrow(res), 2L)
  testthat::expect_setequal(names(res), c("replicate", "timepoint", "analysis", "n", "p_value", "mean"))
  testthat::expect_true(any(is.na(res$p_value)))
  testthat::expect_true(any(is.na(res$mean)))
})

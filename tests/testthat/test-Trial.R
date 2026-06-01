make_pop <- function(name, n_subj = 4, n_read = 2) {
  df <- data.frame(
    id = rep(seq_len(n_subj), each = n_read),
    data = seq_len(n_subj * n_read) * 0.1,
    readout_time = rep(seq_len(n_read) - 1L, n_subj)
  )
  Population$new(name = name, data = df)
}

make_trial <- function(name = "trial",
                       seed = NULL,
                       n_subj = 4,
                       n_read = 1,
                       enroller = n_subj,
                       dropper = 0L,
                       cal_time = 1,
                       analysis = function(df, ct) df) {
  pop <- make_pop("A", n_subj, n_read)
  timer <- Timer$new("t")
  timer$add_timepoint(
    time = cal_time, arm = "A",
    enroller = as.integer(enroller),
    dropper  = as.integer(dropper)
  )
  cal_cond <- trigger_by_calendar(cal_time, analysis = analysis)
  Trial$new(name = name, seed = seed, timer = timer, population = list(pop), conditions = list(cal_cond))
}

### initialize() ###

test_that("Trial initialize: timer must be a Timer instance if supplied", {
  testthat::expect_error(
    Trial$new(name = "x", timer = list(foo = 1), population = list()),
    "`timer` must be a Timer instance."
  )
})

test_that("Trial initialize: errors when populations have different n_readouts", {
  popA <- Population$new("A", data = data.frame(
    id = 1:3, data = 1:3, readout_time = 0
  ))
  popB <- Population$new("B", data = data.frame(
    id = rep(1:3, each = 2), data = rep(1:3, each = 2), readout_time = rep(c(0, 1), 3)
  ))

  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 3L, dropper = 0L)
  timer$add_timepoint(time = 1, arm = "B", enroller = 3L, dropper = 0L)

  testthat::expect_error(
    Trial$new(name = "bad_trial", timer = timer, population = list(popA, popB)),
    "All populations must have the same n_readouts"
  )
})

test_that("Trial initialize: errors when neither timer nor population has enrollment", {
  timer <- Timer$new("t")
  pop <- make_pop("A", 3, 1)

  testthat::expect_error(
    Trial$new(name = "no_data", timer = timer, population = list(pop)),
    "Neither Timer nor Population has enrollment data"
  )
})

test_that("Trial initialize: auto-builds timer from pre-enrolled population", {
  pop <- make_pop("A", 3, 1)
  pop$enrolled <- c(1, 1, 2)

  trial <- Trial$new(name = "auto_timer", timer = NULL, population = list(pop))

  testthat::expect_false(is.null(trial$timer))
  testthat::expect_true(length(trial$timer$timelist) > 0)
})

### run() ###

test_that("Trial run: errors when population list is empty", {
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 3L, dropper = 0L)
  trial <- Trial$new("empty_pop", timer = timer, population = list())
  testthat::expect_error(
    trial$run(),
    "Timer and population list must be set before running run()"
  )
})

test_that("Trial run: locked_data and results remain empty when no conditions fire", {
  pop <- make_pop("A", 4, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 4L, dropper = 0L)
  trial <- Trial$new("no_cond", timer = timer, population = list(pop))
  trial$run()

  testthat::expect_equal(trial$locked_data, list())
  testthat::expect_equal(trial$results, list())
})

test_that("Trial run: measurement_time equals readout_time + enroll_time", {
  trial <- make_trial("meas_time", seed = 42, n_subj = 4, n_read = 2)
  trial$run()

  snap <- trial$locked_data[["time_1"]]
  testthat::expect_equal(snap$measurement_time, snap$readout_time + snap$enroll_time)
})

test_that("Trial run: only enrolled subjects appear in snapshot", {
  trial <- make_trial("partial_enroll", seed = 1, n_subj = 6, n_read = 1, enroller = 3L)
  trial$run()

  snap <- trial$locked_data[["time_1"]]
  testthat::expect_equal(nrow(snap), 3L)
})

test_that("Trial run: snapshot grows cumulatively across timepoints", {
  pop <- make_pop("A", 6, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 3L, dropper = 0L)
  timer$add_timepoint(time = 2, arm = "A", enroller = 3L, dropper = 0L)
  cal_cond_1 <- trigger_by_calendar(1, analysis = function(df, ct) df)
  cal_cond_2 <-trigger_by_calendar(2, analysis = function(df, ct) df)

  trial <- Trial$new("cumulative", seed = 1, timer = timer, population = list(pop), conditions = list(cal_cond_1, cal_cond_2))
  trial$run()

  testthat::expect_equal(nrow(trial$locked_data[["time_1"]]), 3L)
  testthat::expect_equal(nrow(trial$locked_data[["time_2"]]), 6L)
})

test_that("Trial run: dropped subjects remain in snapshot with non-NA drop_time", {
  pop <- make_pop("A", 4, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 4L, dropper = 0L)
  timer$add_timepoint(time = 2, arm = "A", enroller = 0L, dropper = 2L)
  cal_cond_1 <- trigger_by_calendar(2, analysis = function(df, ct) df)

  trial <- Trial$new("dropout", seed = 1, timer = timer, population = list(pop), conditions = list(cal_cond_1))
  trial$run()

  snap <- trial$locked_data[["time_2"]]
  testthat::expect_equal(nrow(snap), 4L)
  testthat::expect_equal(sum(!is.na(snap$drop_time)), 2L)
})

test_that("Trial run: analysis function result is stored in results", {
  analysis <- function(df, ct) list(n_enrolled = nrow(df), cal_time = ct)
  trial <- make_trial("analysis_result", seed = 1, analysis = analysis)
  trial$run()

  testthat::expect_true("time_1" %in% names(trial$results))
  res <- trial$results[["time_1"]]
  testthat::expect_true(is.list(res))
})

test_that("Trial run: subject_id globally unique across arms with same n_readouts", {
  # Arm A: 3 subjects, 2 readouts each -> 6 rows
  dataA <- data.frame(
    id           = rep(1:3, each = 2),
    data         = c(1, 1, 2, 2, 3, 3),
    readout_time = rep(c(0, 1), 3)
  )
  popA <- Population$new("A", data = dataA)

  # Arm B: 3 subjects, 2 readouts each -> 6 rows
  dataB <- data.frame(
    id           = rep(1:3, each = 2),
    data         = c(4, 4, 5, 5, 6, 6),
    readout_time = rep(c(0, 1), 3)
  )
  popB <- Population$new("B", data = dataB)

  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 3L, dropper = 0L)
  timer$add_timepoint(time = 1, arm = "B", enroller = 3L, dropper = 0L)

  cal_cond_1 <- trigger_by_calendar(1, analysis = function(df, current_time) df)

  trial <- Trial$new(
    name       = "two_arm_repeated",
    seed       = 1,
    timer      = timer,
    population = list(popA, popB),
    conditions = list(cal_cond_1)
  )

  trial$run()

  snap <- trial$locked_data[["time_1"]]

  # 6 subjects x 2 readouts = 12 rows
  testthat::expect_equal(nrow(snap), 12L)

  # 6 globally unique subject IDs
  testthat::expect_equal(snap$subject_id, rep(1:6, each = 2))

  # Each subject ID appears exactly 2 times (each, not cycled with times)
  testthat::expect_true(all(table(snap$subject_id) == 2L))
})

test_that("Trial run: same seed produces identical enrolled and dropped assignments", {
  make_seeded_trial <- function() {
    pop <- make_pop("A", 6, 1)
    timer <- Timer$new("t")
    timer$add_timepoint(time = 1, arm = "A", enroller = 4L, dropper = 2L)
    cal_cond_1 <- trigger_by_calendar(1, analysis = function(df, ct) df)
    Trial$new("repro", seed = 7654, timer = timer, population = list(pop), conditions = list(cal_cond_1))
  }

  t1 <- make_seeded_trial()
  t1$run()
  t2 <- make_seeded_trial()
  t2$run()

  snap1 <- t1$locked_data[["time_1"]]
  snap2 <- t2$locked_data[["time_1"]]

  testthat::expect_equal(snap1$enroll_time, snap2$enroll_time)
  testthat::expect_equal(snap1$drop_time, snap2$drop_time)
})

test_that("Trial run: returns invisible self for method chaining", {
  trial <- make_trial("chain", seed = 1)
  ret <- trial$run()
  testthat::expect_identical(ret, trial)
})

test_that("Trial run: second run() call is idempotent", {
  trial <- make_trial("idem", seed = 1, n_subj = 4, dropper = 2L)
  trial$run()
  snap_first <- trial$locked_data[["time_1"]]
  results_first <- trial$results[["time_1"]]

  trial$run()
  snap_second <- trial$locked_data[["time_1"]]
  results_second <- trial$results[["time_1"]]

  testthat::expect_equal(snap_first, snap_second)
  testthat::expect_equal(results_first, results_second)
})

test_that("Trial run: enrollment capped when enroller exceeds available subjects", {
  trial <- make_trial("cap", seed = 1, n_subj = 3, enroller = 10L)
  trial$run()

  snap <- trial$locked_data[["time_1"]]
  testthat::expect_equal(nrow(snap), 3L)
  testthat::expect_true(all(!is.na(snap$enroll_time)))
})

test_that("Trial run: arm column preserved in multi-arm snapshot", {
  popA <- make_pop("A", 3, 1)
  popB <- make_pop("B", 3, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 3L, dropper = 0L)
  timer$add_timepoint(time = 1, arm = "B", enroller = 3L, dropper = 0L)
  cal_cond_1 <- trigger_by_calendar(1, analysis = function(df, ct) df)

  trial <- Trial$new("arms", seed = 1, timer = timer, population = list(popA, popB), conditions = list(cal_cond_1))
  trial$run()

  snap <- trial$locked_data[["time_1"]]
  testthat::expect_true("arm" %in% names(snap))
  testthat::expect_setequal(snap$arm, rep(c("A", "B"), each = 3))
})

test_that("Trial run: drop_time >= enroll_time for all dropped subjects", {
  trial <- make_trial("inv", seed = 1, n_subj = 6, dropper = 0L)

  trial$timer$add_timepoint(time = 2, arm = "A", enroller = 0L, dropper = 3L)
  cal_cond_1 <- trigger_by_calendar(2, analysis = function(df, ct) df)
  trial$conditions <- append(trial$conditions, cal_cond_1)
  trial$run()

  snap <- trial$locked_data[["time_2"]]
  dropped <- snap[!is.na(snap$drop_time), ]
  testthat::expect_true(nrow(dropped) > 0)
  testthat::expect_true(all(dropped$drop_time >= dropped$enroll_time))
})

test_that("Trial run: timepoints processed in sorted order regardless of insertion order", {
  pop <- make_pop("A", 6, 1)

  timer_rev <- Timer$new("t_rev")
  timer_rev$add_timepoint(time = 3, arm = "A", enroller = 3L, dropper = 0L)
  timer_rev$add_timepoint(time = 1, arm = "A", enroller = 3L, dropper = 0L)
  cal_cond_1 <- trigger_by_calendar(1, analysis = function(df, ct) df)
  cal_cond_2 <- trigger_by_calendar(3, analysis = function(df, ct) df)

  trial <- Trial$new("sort_check", seed = 123, timer = timer_rev, population = list(pop), conditions = list(cal_cond_1, cal_cond_2))
  trial$run()

  testthat::expect_equal(nrow(trial$locked_data[["time_1"]]), 3L)
  testthat::expect_equal(nrow(trial$locked_data[["time_3"]]), 6L)
})

test_that("Trial run: duplicate time/arm timepoint rows are aggregated", {
  pop <- make_pop("A", 6, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 2L, dropper = 0L)
  timer$add_timepoint(time = 1, arm = "A", enroller = 2L, dropper = 0L)
  cal_cond_1 <- trigger_by_calendar(1, analysis = function(df, ct) df)

  trial <- Trial$new("agg", seed = 1, timer = timer, population = list(pop), conditions = list(cal_cond_1))
  trial$run()

  snap <- trial$locked_data[["time_1"]]
  testthat::expect_equal(nrow(snap), 4L)
})

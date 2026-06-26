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
                       enroll = n_subj,
                       drop = 0L,
                       cal_time = 1,
                       analysis = function(df, ct) df) {
  pop <- make_pop("A", n_subj, n_read)
  timer <- Timer$new("t")
  timer$add_timepoint(
    time = cal_time, arm = "A",
    enroll = as.integer(enroll),
    drop   = as.integer(drop)
  )
  cal_cond <- condition_calendar_time(cal_time, analysis = analysis)
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
  timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
  timer$add_timepoint(time = 1, arm = "B", enroll = 3L, drop = 0L)

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
  testthat::expect_true(nrow(trial$timer$timelist) > 0L)
})

### run() ###

test_that("Trial run: errors when population list is empty", {
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
  trial <- Trial$new("empty_pop", timer = timer, population = list())
  testthat::expect_error(
    trial$run(),
    "Timer and population list must be set before running run()"
  )
})

test_that("Trial run: locked_data and results remain empty when no conditions fire", {
  pop <- make_pop("A", 4, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 4L, drop = 0L)
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

test_that("Trial run: empty snapshot timepoints are skipped and later snapshots still store", {
  timer <- Timer$new("t")
  timer$add_timepoint(time = 0, arm = "A", enroll = 0L, drop = 0L)
  timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
  pop <- Population$new("A", as_population_data(rnorm(5)))
  cal_cond_0 <- condition_calendar_time(0, analysis = function(df, ct) df)
  cal_cond_1 <- condition_calendar_time(1, analysis = function(df, ct) df)
  trial <- Trial$new(
    name = "empty_snapshot_guard",
    timer = timer,
    population = list(pop),
    conditions = list(cal_cond_0, cal_cond_1)
  )

  testthat::expect_no_error(trial$run())
  testthat::expect_false("time_0" %in% names(trial$locked_data))
  testthat::expect_false("time_0" %in% names(trial$results))
  testthat::expect_true("time_1" %in% names(trial$locked_data))
  testthat::expect_true("time_1" %in% names(trial$results))
  testthat::expect_equal(nrow(trial$locked_data[["time_1"]]), 3L)
})

test_that("Trial run: only enrolled subjects appear in snapshot", {
  trial <- make_trial("partial_enroll", seed = 1, n_subj = 6, n_read = 1, enroll = 3L)
  trial$run()

  snap <- trial$locked_data[["time_1"]]
  testthat::expect_equal(nrow(snap), 3L)
})

test_that("Trial run: snapshot grows cumulatively across timepoints", {
  pop <- make_pop("A", 6, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
  timer$add_timepoint(time = 2, arm = "A", enroll = 3L, drop = 0L)
  cal_cond_1 <- condition_calendar_time(1, analysis = function(df, ct) df)
  cal_cond_2 <- condition_calendar_time(2, analysis = function(df, ct) df)

  trial <- Trial$new("cumulative", seed = 1, timer = timer, population = list(pop), conditions = list(cal_cond_1, cal_cond_2))
  trial$run()

  testthat::expect_equal(nrow(trial$locked_data[["time_1"]]), 3L)
  testthat::expect_equal(nrow(trial$locked_data[["time_2"]]), 6L)
})

test_that("Trial run: dropped subjects remain in snapshot with non-NA drop_time", {
  pop <- make_pop("A", 4, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 4L, drop = 0L)
  timer$add_timepoint(time = 2, arm = "A", enroll = 0L, drop = 2L)
  cal_cond_1 <- condition_calendar_time(2, analysis = function(df, ct) df)

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
  timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
  timer$add_timepoint(time = 1, arm = "B", enroll = 3L, drop = 0L)

  cal_cond_1 <- condition_calendar_time(1, analysis = function(df, current_time) df)

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
    timer$add_timepoint(time = 1, arm = "A", enroll = 4L, drop = 2L)
    cal_cond_1 <- condition_calendar_time(1, analysis = function(df, ct) df)
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
  trial <- make_trial("idem", seed = 1, n_subj = 4, drop = 2L)
  trial$run()
  snap_first <- trial$locked_data[["time_1"]]
  results_first <- trial$results[["time_1"]]

  trial$run()
  snap_second <- trial$locked_data[["time_1"]]
  results_second <- trial$results[["time_1"]]

  testthat::expect_equal(snap_first, snap_second)
  testthat::expect_equal(results_first, results_second)
})

test_that("Trial run: enrollment capped when enroll exceeds available subjects", {
  trial <- make_trial("cap", seed = 1, n_subj = 3, enroll = 10L)
  trial$run()

  snap <- trial$locked_data[["time_1"]]
  testthat::expect_equal(nrow(snap), 3L)
  testthat::expect_true(all(!is.na(snap$enroll_time)))
})

test_that("Trial run: arm column preserved in multi-arm snapshot", {
  popA <- make_pop("A", 3, 1)
  popB <- make_pop("B", 3, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
  timer$add_timepoint(time = 1, arm = "B", enroll = 3L, drop = 0L)
  cal_cond_1 <- condition_calendar_time(1, analysis = function(df, ct) df)

  trial <- Trial$new("arms", seed = 1, timer = timer, population = list(popA, popB), conditions = list(cal_cond_1))
  trial$run()

  snap <- trial$locked_data[["time_1"]]
  testthat::expect_true("arm" %in% names(snap))
  testthat::expect_setequal(snap$arm, rep(c("A", "B"), each = 3))
})

test_that("Trial run: drop_time >= enroll_time for all dropped subjects", {
  trial <- make_trial("inv", seed = 1, n_subj = 6, drop = 0L)

  trial$timer$add_timepoint(time = 2, arm = "A", enroll = 0L, drop = 3L)
  cal_cond_1 <- condition_calendar_time(2, analysis = function(df, ct) df)
  trial$conditions <- append(trial$conditions, cal_cond_1)
  trial$run()

  snap <- trial$locked_data[["time_2"]]
  dropped <- snap[!is.na(snap$drop_time), ]
  testthat::expect_true(nrow(dropped) > 0)
  testthat::expect_true(all(dropped$drop_time >= dropped$enroll_time))
})

test_that("Trial run: output is invariant to timepoint insertion order", {
  # The engine sorts timepoints internally, so scrambled insertion must yield
  # the same snapshots as ascending insertion. Compare the two directly: a
  # snapshot keyed by insertion index instead of time would diverge here.
  build <- function(times) {
    pop <- make_pop("A", 6, 1)
    timer <- Timer$new("t")
    for (tm in times) timer$add_timepoint(time = tm, arm = "A", enroll = 3L, drop = 0L)
    cal_cond_1 <- condition_calendar_time(1, analysis = function(df, ct) df)
    cal_cond_2 <- condition_calendar_time(3, analysis = function(df, ct) df)
    trial <- Trial$new("order", seed = 123, timer = timer,
                       population = list(pop), conditions = list(cal_cond_1, cal_cond_2))
    trial$run()
    trial$locked_data
  }
  scrambled <- build(c(3, 1))
  ascending <- build(c(1, 3))

  testthat::expect_equal(scrambled[["time_1"]], ascending[["time_1"]])
  testthat::expect_equal(scrambled[["time_3"]], ascending[["time_3"]])
})

test_that("Trial run: duplicate time/arm timepoint rows are aggregated", {
  pop <- make_pop("A", 6, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 2L, drop = 0L)
  timer$add_timepoint(time = 1, arm = "A", enroll = 2L, drop = 0L)
  cal_cond_1 <- condition_calendar_time(1, analysis = function(df, ct) df)

  trial <- Trial$new("agg", seed = 1, timer = timer, population = list(pop), conditions = list(cal_cond_1))
  trial$run()

  snap <- trial$locked_data[["time_1"]]
  testthat::expect_equal(nrow(snap), 4L)
})

### adaptive = FALSE (fixed fast path) ###

test_that("fixed path: adaptive flag defaults to FALSE", {
  trial <- make_trial("default_adaptive")
  testthat::expect_false(trial$adaptive)
})

test_that("fixed path: enroll_times sorted ascending (deterministic)", {
  pop <- make_pop("A", 6, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
  timer$add_timepoint(time = 2, arm = "A", enroll = 3L, drop = 0L)
  cal_cond <- condition_calendar_time(2, analysis = function(df, ct) df)
  trial <- Trial$new("fixed_enroll_order", timer = timer, population = list(pop),
                     conditions = list(cal_cond), adaptive = FALSE)
  trial$run()
  snap <- trial$locked_data[["time_2"]]
  testthat::expect_equal(sum(snap$enroll_time == 1), 3L)
  testthat::expect_equal(sum(snap$enroll_time == 2), 3L)
})

test_that("fixed path: drop_time >= enroll_time for all dropped subjects", {
  pop <- make_pop("A", 4, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 4L, drop = 0L)
  timer$add_timepoint(time = 2, arm = "A", enroll = 0L, drop = 2L)
  cal_cond <- condition_calendar_time(2, analysis = function(df, ct) df)
  trial <- Trial$new("fixed_drop_order", timer = timer, population = list(pop),
                     conditions = list(cal_cond), adaptive = FALSE)
  trial$run()
  snap <- trial$locked_data[["time_2"]]
  dropped <- snap[!is.na(snap$drop_time), ]
  testthat::expect_equal(nrow(dropped), 2L)
  testthat::expect_true(all(dropped$drop_time >= dropped$enroll_time))
})

test_that("fixed path: drop_time masked to NA for future drops in earlier snapshots", {
  pop <- make_pop("A", 4, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 4L, drop = 0L)
  timer$add_timepoint(time = 2, arm = "A", enroll = 0L, drop = 2L)
  cal_cond_1 <- condition_calendar_time(1, analysis = function(df, ct) df)
  cal_cond_2 <- condition_calendar_time(2, analysis = function(df, ct) df)
  trial <- Trial$new("fixed_mask", timer = timer, population = list(pop),
                     conditions = list(cal_cond_1, cal_cond_2), adaptive = FALSE)
  trial$run()
  snap_t1 <- trial$locked_data[["time_1"]]
  snap_t2 <- trial$locked_data[["time_2"]]
  testthat::expect_true(all(is.na(snap_t1$drop_time)))
  testthat::expect_equal(sum(!is.na(snap_t2$drop_time)), 2L)
})

test_that("fixed path: prefix snapshots grow cumulatively across timepoints", {
  pop <- make_pop("A", 6, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
  timer$add_timepoint(time = 2, arm = "A", enroll = 3L, drop = 0L)
  cal_cond_1 <- condition_calendar_time(1, analysis = function(df, ct) df)
  cal_cond_2 <- condition_calendar_time(2, analysis = function(df, ct) df)
  trial <- Trial$new("fixed_cumulative", timer = timer, population = list(pop),
                     conditions = list(cal_cond_1, cal_cond_2), adaptive = FALSE)
  trial$run()
  testthat::expect_equal(nrow(trial$locked_data[["time_1"]]), 3L)
  testthat::expect_equal(nrow(trial$locked_data[["time_2"]]), 6L)
})

test_that("fixed path: NULL drops produce all-NA drop_time", {
  pop <- make_pop("A", 4, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 4L, drop = NULL)
  cal_cond <- condition_calendar_time(1, analysis = function(df, ct) df)
  trial <- Trial$new("fixed_null_drop", timer = timer, population = list(pop),
                     conditions = list(cal_cond), adaptive = FALSE)
  trial$run()
  snap <- trial$locked_data[["time_1"]]
  testthat::expect_true(all(is.na(snap$drop_time)))
})

test_that("fixed path: parity with adaptive — same nrow and enrolled count", {
  make_parity_trial <- function(mode) {
    pop <- make_pop("A", 6, 1)
    timer <- Timer$new("t")
    timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
    timer$add_timepoint(time = 2, arm = "A", enroll = 3L, drop = 1L)
    cal_cond <- condition_calendar_time(2, analysis = function(df, ct) df)
    Trial$new("parity", seed = 42, timer = timer, population = list(pop),
              conditions = list(cal_cond), adaptive = mode)
  }
  t_fixed    <- make_parity_trial(FALSE)
  t_adaptive <- make_parity_trial(TRUE)
  t_fixed$run()
  t_adaptive$run()
  snap_f <- t_fixed$locked_data[["time_2"]]
  snap_a <- t_adaptive$locked_data[["time_2"]]
  testthat::expect_equal(nrow(snap_f), nrow(snap_a))
  testthat::expect_equal(sum(!is.na(snap_f$enroll_time)), sum(!is.na(snap_a$enroll_time)))
  testthat::expect_true("time_2" %in% names(t_fixed$results))
  testthat::expect_true("time_2" %in% names(t_adaptive$results))
})

test_that("fixed path: clone_trial carries adaptive flag", {
  trial <- make_trial("src")
  clones <- clone_trial(trial, n = 2)
  testthat::expect_false(clones[[1]]$adaptive)
  testthat::expect_false(clones[[2]]$adaptive)
})

test_that("fixed path: interim + final both fire (mid-gap skip is correct)", {
  pop <- make_pop("A", 10, 1)
  timer <- Timer$new("t")
  for (k in 1:10) timer$add_timepoint(time = k, arm = "A", enroll = 1L, drop = 0L)
  interim <- Condition$new(where = count_trigger("enroll_time", ">=", 5L),
                           analysis = function(df, ct) nrow(df), name = "interim")
  final   <- Condition$new(where = count_trigger("enroll_time", ">=", 10L),
                           analysis = function(df, ct) nrow(df), name = "final")
  trial <- Trial$new("skip_midgap", timer = timer, population = list(pop),
                     conditions = list(interim, final), adaptive = FALSE)
  trial$run()
  testthat::expect_equal(trial$results[["time_5"]][["interim"]], 5L)
  testthat::expect_equal(trial$results[["time_10"]][["final"]], 10L)
})

test_that("fixed path: unknown (non-monotone) trigger still fires via fallback", {
  pop <- make_pop("A", 6, 1)
  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroll = 3L, drop = 0L)
  timer$add_timepoint(time = 2, arm = "A", enroll = 3L, drop = 0L)
  # value_trigger on a data column with `>` is not reasoned about (fire_time=-Inf)
  cond <- Condition$new(where = value_trigger("subject_id", "==", 1),
                        analysis = function(df, ct) nrow(df), name = "any")
  trial <- Trial$new("skip_fallback", timer = timer, population = list(pop),
                     conditions = list(cond), adaptive = FALSE)
  trial$run()
  testthat::expect_true(any(grepl("^time_", names(trial$results))))
})

test_that(".trigger_fire_time: count beyond n never fires (Inf)", {
  testthat::expect_equal(
    rxsim:::.trigger_fire_time(count_trigger("enroll_time", ">=", 99L), subj_enroll = 1:10),
    Inf
  )
  testthat::expect_equal(
    rxsim:::.trigger_fire_time(calendar_trigger(7), subj_enroll = 1:10), 7
  )
  testthat::expect_equal(
    rxsim:::.trigger_fire_time(NULL, subj_enroll = 1:10), -Inf
  )
})

test_that(".trigger_fire_time: combinators and strict ops", {
  e <- 1:10  # subj_enroll: k-th subject enrolls at time k

  # & = max of children (both must hold): calendar>=3, count>=5 -> max(3, 5)
  testthat::expect_equal(
    rxsim:::.trigger_fire_time(
      calendar_trigger(3) & count_trigger("enroll_time", ">=", 5L), e), 5)

  # | = min of children (either fires): calendar>=8, count>=2 -> min(8, 2)
  testthat::expect_equal(
    rxsim:::.trigger_fire_time(
      calendar_trigger(8) | count_trigger("enroll_time", ">=", 2L), e), 2)

  # strict count `>` k uses the (k+1)-th enrollment: >5 -> subj_enroll[6]
  testthat::expect_equal(
    rxsim:::.trigger_fire_time(count_trigger("enroll_time", ">", 5L), e), 6)

  # strict calendar `>` rhs is treated conservatively as rhs (evaluate from rhs)
  testthat::expect_equal(
    rxsim:::.trigger_fire_time(value_trigger("time", ">", 4), e), 4)
})

test_that("fixed path: leading timepoints are actually skipped (not just correct)", {
  SpyCond <- R6::R6Class("SpyCond", inherit = Condition, public = list(
    calls = 0L,
    check_conditions = function(locked_data, current_time) {
      self$calls <- self$calls + 1L
      super$check_conditions(locked_data, current_time)
    }))

  pop <- make_pop("A", 10, 1)
  timer <- Timer$new("t")
  for (k in 1:10) timer$add_timepoint(time = k, arm = "A", enroll = 1L, drop = 0L)
  # 1 enrollment/timepoint => 5th subject enrolls at time 5 => fire_time = 5
  # max_triggers high so it keeps evaluating after the first fire, proving the
  # skip is leading-only (times 1-4 never reached, 5..10 all evaluated).
  spy <- SpyCond$new(where = count_trigger("enroll_time", ">=", 5L),
                     analysis = function(df, ct) nrow(df), name = "spy",
                     max_triggers = 100L)
  trial <- Trial$new("skip_spy", timer = timer, population = list(pop),
                     conditions = list(spy), adaptive = FALSE)
  trial$run()

  # Evaluated only at times 5..10 (6 calls); == 10 would mean no skip happened.
  testthat::expect_equal(spy$calls, 6L)
  testthat::expect_equal(spy$trigger_count, 6L)
})

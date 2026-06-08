# ── value_trigger ────────────────────────────────────────────────────────────

test_that("value_trigger: returns rxsim_trigger with correct fields", {
  t <- value_trigger("time", ">=", 52)

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$type, "value")
  expect_equal(t$col, "time")
  expect_equal(t$op, ">=")
  expect_equal(t$rhs, 52)
})

test_that("value_trigger: accepts vector rhs", {
  t <- value_trigger("time", "%in%", c(26, 52))

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$rhs, c(26, 52))
})

test_that("value_trigger: errors on invalid col", {
  expect_error(value_trigger(123,   ">=", 52),  "`col`")
  expect_error(value_trigger(NA_character_, ">=", 52), "`col`")
  expect_error(value_trigger(c("a", "b"), ">=", 52), "`col`")
})

test_that("value_trigger: errors on invalid op", {
  expect_error(value_trigger("time", "INVALID", 52), "`op`")
  expect_error(value_trigger("time", "=>",      52), "`op`")
})

test_that("value_trigger: errors on non-atomic rhs", {
  expect_error(value_trigger("time", ">=", list(52)), "`rhs`")
})

# ── count_trigger ────────────────────────────────────────────────────────────

test_that("count_trigger: returns rxsim_trigger with correct fields", {
  t <- count_trigger("enroll_time", ">=", 100)

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$type, "count")
  expect_equal(t$col, "enroll_time")
  expect_equal(t$op, ">=")
  expect_equal(t$rhs, 100)
})

test_that("count_trigger: errors on non-numeric rhs", {
  expect_error(count_trigger("enroll_time", ">=", "100"), "`rhs`")
  expect_error(count_trigger("enroll_time", ">=", list(100)), "`rhs`")
})

test_that("count_trigger: errors on invalid col and op", {
  expect_error(count_trigger(123, ">=", 100), "`col`")
  expect_error(count_trigger("enroll_time", "INVALID", 100), "`op`")
})

# ── enroll_trigger ───────────────────────────────────────────────────────────

test_that("enroll_trigger: wraps count_trigger with enroll_time", {
  t <- enroll_trigger(0.5, 200)

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$type, "count")
  expect_equal(t$col, "enroll_time")
  expect_equal(t$op, ">=")
  expect_equal(t$rhs, 100)  # 0.5 * 200
})

test_that("enroll_trigger: errors on fraction out of range", {
  expect_error(enroll_trigger(0,    100), "`fraction`")
  expect_error(enroll_trigger(1.1,  100), "`fraction`")
  expect_error(enroll_trigger(-0.1, 100), "`fraction`")
})

test_that("enroll_trigger: errors on invalid sample_size", {
  expect_error(enroll_trigger(0.5, "100"), "`sample_size`")
  expect_error(enroll_trigger(0.5, c(50, 50)), "`sample_size`")
  expect_error(enroll_trigger(0.5, NA_real_), "`sample_size`")
})

test_that("enroll_trigger: errors when arguments missing", {
  expect_error(enroll_trigger(0.5), "`fraction` and `sample_size`")
})

# ── calendar_trigger ─────────────────────────────────────────────────────────

test_that("calendar_trigger: wraps value_trigger with time column", {
  t <- calendar_trigger(52)

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$type, "value")
  expect_equal(t$col, "time")
  expect_equal(t$op, "%in%")
  expect_equal(t$rhs, 52)
})

test_that("calendar_trigger: accepts vector of times", {
  t <- calendar_trigger(c(26, 52))
  expect_equal(t$rhs, c(26, 52))
})

test_that("calendar_trigger: errors on non-numeric", {
  expect_error(calendar_trigger("52"), "`cal_time`")
  expect_error(calendar_trigger(),     "`cal_time`")
})

# ── &.rxsim_trigger / |.rxsim_trigger ────────────────────────────────────────

test_that("& produces rxsim_trigger with AND combinator", {
  t <- enroll_trigger(0.5, 200) & calendar_trigger(52)

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$combinator, "&")
  expect_length(t$predicates, 2L)
})

test_that("| produces rxsim_trigger with OR combinator", {
  t <- enroll_trigger(0.5, 200) | calendar_trigger(26)

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$combinator, "|")
  expect_s3_class(t$left,  "rxsim_trigger")
  expect_s3_class(t$right, "rxsim_trigger")
})

test_that("& and | error when operand is not rxsim_trigger", {
  t <- enroll_trigger(0.5, 200)
  expect_error(t & "not a trigger", "rxsim_trigger")
  expect_error(t | 42,              "rxsim_trigger")
})

test_that("triggers compose nested: (A & B) | C", {
  t <- (enroll_trigger(0.5, 200) & calendar_trigger(52)) | value_trigger("time", ">=", 100)

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$combinator, "|")
  expect_equal(t$left$combinator, "&")
})

# ── Condition$new with rxsim_trigger ─────────────────────────────────────────

test_that("Condition$new: rxsim_trigger value_trigger produces correct quosure", {
  trig <- value_trigger("time", ">=", 52)
  cond <- Condition$new(where = trig, analysis = function(df, t) nrow(df), name = "test")

  expect_r6_class(cond, "Condition")
  expect_equal(cond$name, "test")
  expect_length(cond$where, 1L)
  expect_true(is.function(cond$analysis))
})

test_that("Condition$new: rxsim_trigger count_trigger produces correct quosure", {
  trig <- count_trigger("enroll_time", ">=", 50)
  cond <- Condition$new(where = trig, name = "enrolled_50")

  expect_r6_class(cond, "Condition")
  expect_length(cond$where, 1L)
})

test_that("Condition$new: AND trigger produces two quosures (dplyr ANDs them)", {
  trig <- enroll_trigger(0.5, 200) & calendar_trigger(52)
  cond <- Condition$new(where = trig, name = "and_test")

  expect_r6_class(cond, "Condition")
  expect_length(cond$where, 2L)  # dplyr::filter(..., pred1, pred2) ANDs them
})

test_that("Condition$new: OR trigger produces one quosure with | expression", {
  trig <- enroll_trigger(0.5, 200) | calendar_trigger(26)
  cond <- Condition$new(where = trig, name = "or_test")

  expect_r6_class(cond, "Condition")
  expect_length(cond$where, 1L)
  expect_equal(rlang::call_name(rlang::get_expr(cond$where[[1L]])), "|")
})

test_that("Condition$new: cooldown and max_triggers are passed through", {
  trig <- calendar_trigger(52)
  cond <- Condition$new(where = trig, cooldown = 5, max_triggers = 3L)

  expect_equal(cond$cooldown, 5)
  expect_equal(cond$max_triggers, 3L)
})

# ── Condition$new: quosure correctness (evaluated against snapshot) ───────────

test_that("Condition$new value_trigger quosure evaluates correctly in filter", {
  trig <- value_trigger("time", ">=", 3)
  cond <- Condition$new(where = trig)
  df   <- data.frame(time = 1:5)

  result <- dplyr::filter(df, !!!cond$where)
  expect_equal(result$time, 3:5)
})

test_that("Condition$new count_trigger quosure evaluates correctly in filter", {
  trig <- count_trigger("enroll_time", ">=", 2)
  cond <- Condition$new(where = trig)
  df   <- data.frame(enroll_time = c(NA, 1, 2, 3), x = 1:4)

  result <- dplyr::filter(df, !!!cond$where)
  expect_equal(nrow(result), 4L)
})

test_that("Condition$new AND quosures: both conditions must hold", {
  trig <- count_trigger("enroll_time", ">=", 1) & value_trigger("time", ">=", 3)
  cond <- Condition$new(where = trig)
  df   <- data.frame(time = 1:5, enroll_time = c(1, NA, 1, 1, 1))

  result <- dplyr::filter(df, !!!cond$where)
  expect_true(all(result$time >= 3))
})

test_that("Condition$new OR quosure: either condition fires", {
  trig <- value_trigger("time", ">=", 4) | value_trigger("time", "<=", 2)
  cond <- Condition$new(where = trig)
  df   <- data.frame(time = 1:5)

  result <- dplyr::filter(df, !!!cond$where)
  expect_equal(result$time, c(1L, 2L, 4L, 5L))
})

# ── trigger_by_calendar ───────────────────────────────────────────────────────

test_that("trigger_by_calendar: returns Condition firing at correct time", {
  cond <- trigger_by_calendar(cal_time = 3, analysis = function(df, t) nrow(df))

  df <- data.frame(time = 1:5, enroll_time = rep(1, 5))
  expect_length(cond$check_conditions(df, 3L), 1L)
  expect_length(cond$check_conditions(df, 4L), 0L)
})

test_that("trigger_by_calendar: default name contains cal_time", {
  cond <- trigger_by_calendar(52)
  expect_match(cond$name, "52")
})

test_that("trigger_by_calendar: custom name is used", {
  cond <- trigger_by_calendar(52, name = "final")
  expect_equal(cond$name, "final")
})

test_that("trigger_by_calendar: analysis function is stored", {
  fn   <- function(df, t) nrow(df)
  cond <- trigger_by_calendar(52, analysis = fn)
  expect_identical(cond$analysis, fn)
})

# ── trigger_by_fraction ───────────────────────────────────────────────────────

test_that("trigger_by_fraction: returns Condition firing at correct enrollment", {
  cond <- trigger_by_fraction(fraction = 0.5, sample_size = 4, analysis = function(df, t) nrow(df))

  df_below <- data.frame(enroll_time = c(1, NA, NA, NA))   # 1/4 enrolled
  df_above <- data.frame(enroll_time = c(1,  2, NA, NA))   # 2/4 enrolled

  expect_length(cond$check_conditions(df_below, 1), 0L)
  expect_length(cond$check_conditions(df_above, 2), 1L)
})

test_that("trigger_by_fraction: default name contains fraction", {
  cond <- trigger_by_fraction(0.5, 100)
  expect_match(cond$name, "0.5")
})

test_that("trigger_by_fraction: errors on missing arguments", {
  expect_error(trigger_by_fraction(0.5),       "`fraction` and `sample_size`")
  expect_error(trigger_by_fraction(sample_size = 100), "`fraction` and `sample_size`")
})

# ── replicate_trial: security gate ───────────────────────────────────────────

make_pop_gen <- function() {
  list(A = function(n) data.frame(id = 1:n, value = rnorm(n), readout_time = 0))
}

test_that("replicate_trial: rejects raw quosures as trigger", {
  ag <- list(final = list(
    trigger  = rlang::quos(sum(!is.na(enroll_time)) >= 10),
    analysis = function(df, t) nrow(df)
  ))
  expect_error(
    replicate_trial("t", 10L, "A", 1, function(n) rexp(n, 1), function(n) rep(Inf, n), ag, make_pop_gen(), 1L),
    "rxsim_trigger"
  )
})

test_that("replicate_trial: accepts rxsim_trigger and builds Condition", {
  set.seed(1)
  ag <- list(final = list(
    trigger  = enroll_trigger(1.0, 10L),
    analysis = function(df, t) nrow(df)
  ))
  trials <- replicate_trial("t", 10L, "A", 1, function(n) rexp(n, 1), function(n) rep(Inf, n), ag, make_pop_gen(), 1L)

  expect_length(trials, 1L)
  expect_r6_class(trials[[1L]], "Trial")
  expect_length(trials[[1L]]$conditions, 1L)
  expect_r6_class(trials[[1L]]$conditions[[1L]], "Condition")
})

# ── notna_trigger ─────────────────────────────────────────────────────────────

test_that("notna_trigger: returns rxsim_trigger with correct fields", {
  t <- notna_trigger("enroll_time")

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$type, "notna")
  expect_equal(t$col, "enroll_time")
})

test_that("notna_trigger: errors on invalid col", {
  expect_error(notna_trigger(123),            "`col`")
  expect_error(notna_trigger(NA_character_),  "`col`")
  expect_error(notna_trigger(c("a", "b")),    "`col`")
})

test_that("notna_trigger quosure evaluates correctly in filter", {
  trig <- notna_trigger("enroll_time")
  cond <- Condition$new(where = trig)
  df   <- data.frame(enroll_time = c(NA, 1, NA, 3))

  result <- dplyr::filter(df, !!!cond$where)
  expect_equal(result$enroll_time, c(1, 3))
})

# ── col_trigger ───────────────────────────────────────────────────────────────

test_that("col_trigger: returns rxsim_trigger with correct fields", {
  t <- col_trigger("enroll_time", "<=", "time")

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$type, "col_compare")
  expect_equal(t$col, "enroll_time")
  expect_equal(t$op, "<=")
  expect_equal(t$ref_col, "time")
})

test_that("col_trigger: errors on invalid arguments", {
  expect_error(col_trigger(123,   "<=", "time"),         "`col`")
  expect_error(col_trigger("a",   "INVALID", "time"),    "`op`")
  expect_error(col_trigger("a",   "<=", NA_character_),  "`ref_col`")
  expect_error(col_trigger("a",   "<=", c("x", "y")),   "`ref_col`")
})

test_that("col_trigger quosure evaluates correctly in filter", {
  trig <- col_trigger("enroll_time", "<=", "time")
  cond <- Condition$new(where = trig)
  df   <- data.frame(enroll_time = c(1, 5, 10), time = rep(6, 3))

  result <- dplyr::filter(df, !!!cond$where)
  expect_equal(result$enroll_time, c(1, 5))
})

# ── timed_count_trigger ───────────────────────────────────────────────────────

test_that("timed_count_trigger: returns rxsim_trigger with correct fields", {
  t <- timed_count_trigger("pfs_event_time", "time", ">=", 100)

  expect_s3_class(t, "rxsim_trigger")
  expect_equal(t$type, "timed_count")
  expect_equal(t$col, "pfs_event_time")
  expect_equal(t$time_col, "time")
  expect_equal(t$op, ">=")
  expect_equal(t$threshold, 100)
})

test_that("timed_count_trigger: errors on invalid arguments", {
  expect_error(timed_count_trigger(123, "time", ">=", 100),               "`col`")
  expect_error(timed_count_trigger("ev", c("a","b"), ">=", 100),          "`time_col`")
  expect_error(timed_count_trigger("ev", "time", "INVALID", 100),         "`op`")
  expect_error(timed_count_trigger("ev", "time", ">=", "100"),            "`threshold`")
  expect_error(timed_count_trigger("ev", "time", ">=", c(100, 200)),      "`threshold`")
})

test_that("timed_count_trigger quosure evaluates correctly: counts events <= time", {
  trig <- timed_count_trigger("pfs_event_time", "time", ">=", 3)
  cond <- Condition$new(where = trig)

  # 2 events at or before time=8: should NOT fire
  df_below <- data.frame(
    pfs_event_time = c(3, 7, NA, NA),
    time           = rep(8, 4)
  )
  result_below <- dplyr::filter(df_below, !!!cond$where)
  expect_equal(nrow(result_below), 0L)

  # 3 events at or before time=8: should fire (all rows pass)
  df_above <- data.frame(
    pfs_event_time = c(3, 7, NA, 5),
    time           = rep(8, 4)
  )
  result_above <- dplyr::filter(df_above, !!!cond$where)
  expect_equal(nrow(result_above), 4L)
})

test_that("timed_count_trigger quosure: events after current time do not count", {
  trig <- timed_count_trigger("event_time", "time", ">=", 2)
  cond <- Condition$new(where = trig)

  # event at t=15 should NOT count when time=10
  df <- data.frame(event_time = c(5, 15, NA), time = rep(10, 3))
  result <- dplyr::filter(df, !!!cond$where)
  expect_equal(nrow(result), 0L)  # only 1 event at t<=10, threshold=2
})

# ── trigger_by_events ─────────────────────────────────────────────────────────

test_that("trigger_by_events: returns Condition with correct default name", {
  cond <- trigger_by_events("pfs_event_time", 100)
  expect_r6_class(cond, "Condition")
  expect_equal(cond$name, "events_100")
})

test_that("trigger_by_events: custom name is used", {
  cond <- trigger_by_events("pfs_event_time", 100, name = "pfs_ia")
  expect_equal(cond$name, "pfs_ia")
})

test_that("trigger_by_events: errors on missing required arguments", {
  expect_error(trigger_by_events(n_events = 100),                "`event_col` and `n_events`")
  expect_error(trigger_by_events(event_col = "pfs_event_time"),  "`event_col` and `n_events`")
})

test_that("trigger_by_events: errors on invalid arguments", {
  expect_error(trigger_by_events(123, 100),            "`event_col`")
  expect_error(trigger_by_events("ev", "100"),          "`n_events`")
  expect_error(trigger_by_events("ev", 100, op = "!"), "`op`")
})

test_that("trigger_by_events: produces 3 quosures (ANDed)", {
  cond <- trigger_by_events("pfs_event_time", 50)
  expect_length(cond$where, 3L)
})

test_that("trigger_by_events: fires when event threshold is reached", {
  cond <- trigger_by_events(
    event_col = "pfs_event_time",
    n_events  = 3,
    analysis  = function(df, t) data.frame(n = nrow(df), fired_at = t)
  )

  # 2 events at time <= 10: should NOT fire
  df_below <- data.frame(
    time           = rep(10, 4),
    enroll_time    = c(1, 2, 3, 4),
    pfs_event_time = c(5, 8, NA, NA)
  )
  expect_length(cond$check_conditions(df_below, 10), 0L)

  # 3 events at time <= 10: should fire
  df_above <- data.frame(
    time           = rep(10, 4),
    enroll_time    = c(1, 2, 3, 4),
    pfs_event_time = c(5, 8, 9, NA)
  )
  res <- cond$check_conditions(df_above, 10)
  expect_length(res, 1L)
  expect_equal(res[["events_3"]]$fired_at, 10)
})

test_that("trigger_by_events: unenrolled subjects are excluded from filtered result", {
  cond <- trigger_by_events(
    event_col = "pfs_event_time",
    n_events  = 2,
    analysis  = function(df, t) data.frame(n_enrolled = nrow(df))
  )

  # 3 events but 1 subject unenrolled (NA enroll_time); threshold=2 met
  df <- data.frame(
    time           = rep(10, 4),
    enroll_time    = c(1, 2, NA, 4),
    pfs_event_time = c(5, 8, 6, NA)
  )
  res <- cond$check_conditions(df, 10)
  expect_length(res, 1L)
  # analysis df contains only enrolled subjects enrolled before time
  expect_equal(res[["events_2"]]$n_enrolled, 3L)
})

test_that("trigger_by_events: events after current time do not count towards threshold", {
  cond <- trigger_by_events(
    event_col = "pfs_event_time",
    n_events  = 2,
    analysis  = function(df, t) data.frame(n = nrow(df))
  )

  # 1 event at t=5, 1 event at t=15 (future) — only 1 counts
  df <- data.frame(
    time           = rep(10, 3),
    enroll_time    = c(1, 2, 3),
    pfs_event_time = c(5, 15, NA)
  )
  expect_length(cond$check_conditions(df, 10), 0L)
})

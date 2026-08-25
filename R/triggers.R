.trigger_ops <- c(">=", "<=", ">", "<", "==", "!=", "%in%")

# Earliest calendar time at which a trigger spec *could* fire, given the sorted
# per-subject enrollment times. Used by the fixed path (Trial$run, adaptive =
# FALSE) to skip leading timepoints where no condition can fire.
#
# Returns:
#   -Inf  : cannot reason about this trigger -> evaluate at every timepoint (safe)
#   Inf   : threshold can never be reached -> the condition never fires
#   t     : earliest time the trigger's gate-1 (non-empty filter) can hold
#
# Monotone helper triggers only (`>=`/`>` on calendar time or counts, and their
# `&`/`|` combinations). Anything else falls back to -Inf.
.trigger_fire_time <- function(spec, subj_enroll) {
  if (is.null(spec)) return(-Inf)

  if (!is.null(spec$combinator)) {
    if (identical(spec$combinator, "&")) {
      parts <- if (!is.null(spec$predicates)) spec$predicates else list(spec$left, spec$right)
      return(max(vapply(parts, .trigger_fire_time, numeric(1), subj_enroll)))
    }
    return(min(
      .trigger_fire_time(spec$left, subj_enroll),
      .trigger_fire_time(spec$right, subj_enroll)
    ))
  }

  op <- spec$op
  rhs <- spec$rhs

  # calendar time: prefix$time == i is constant, so fires at first i >= rhs
  if (identical(spec$type, "value") && identical(spec$col, "time") && op %in% c(">=", ">")) {
    return(rhs)
  }

  # count threshold: fires when the rhs-th subject has enrolled
  if (identical(spec$type, "count") && op %in% c(">=", ">")) {
    k <- if (identical(op, ">=")) ceiling(rhs) else floor(rhs) + 1
    if (k < 1) return(-Inf)
    if (k > length(subj_enroll)) return(Inf)
    return(subj_enroll[k])
  }

  -Inf
}

.check_col_op <- function(col, op) {
  if (!is.character(col) || length(col) != 1L || is.na(col)) {
    stop("`col` must be a single character string.")
  }
# Validation shared by the trigger constructors. `arg` names the parameter so
# the message still points at `ref_col`/`time_col`, not a generic "col".
.check_string <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single character string.", arg))
  }
}

.check_op <- function(op) {
  if (!is.character(op) || length(op) != 1L || is.na(op) || !op %in% .trigger_ops) {
    stop("`op` must be one of: >=, <=, >, <, ==, !=, %in%.")
  }
}

#' @name trigger_primitives
#' @title Build Trial Triggers
#' @description Create trigger specifications that can be passed to
#'   [`Condition`] or composed with `&` and `|`.
#'
#' @param col `character` Column name referenced by the trigger.
#' @param op `character` Comparison operator. Must be one of
#'   `c(">=", "<=", ">", "<", "==", "!=", "%in%")`.
#' @param rhs Right-hand side value. Must be atomic for `value_trigger()` and
#'   numeric for `count_trigger()`.
#'
#' @return A `trigger` object.
#'
#' @seealso [Condition], `enroll_trigger()`, `calendar_trigger()`.
#'
#' @export
#'
#' @examples
#' t1 <- value_trigger("time", ">=", 52)
#' t2 <- count_trigger("enroll_time", ">=", 100)
#' t3 <- enroll_trigger(0.5, 200) & calendar_trigger(52)
value_trigger <- function(col, op, rhs) {
  .check_string(col, "col")
  .check_op(op)
  if (!is.atomic(rhs)) {
    stop("`rhs` must be atomic for `value_trigger()`.")
  }

  structure(list(type = "value", col = col, op = op, rhs = rhs), class = "trigger")
}

#' @rdname trigger_primitives
#' @export
count_trigger <- function(col, op, rhs) {
  .check_string(col, "col")
  .check_op(op)
  if (!is.atomic(rhs) || !is.numeric(rhs)) {
    stop("`rhs` must be numeric for `count_trigger()`.")
  }

  structure(list(type = "count", col = col, op = op, rhs = rhs), class = "trigger")
}

#' @rdname trigger_primitives
#' @export
notna_trigger <- function(col) {
  .check_string(col, "col")
  structure(list(type = "notna", col = col), class = "trigger")
}

#' @rdname trigger_primitives
#' @param ref_col `character` Column name on the right-hand side of the
#'   comparison.
#' @export
col_trigger <- function(col, op, ref_col) {
  .check_string(col, "col")
  .check_op(op)
  .check_string(ref_col, "ref_col")
  structure(list(type = "col_compare", col = col, op = op, ref_col = ref_col), class = "trigger")
}

#' @rdname trigger_primitives
#' @param time_col `character` Column name of the current-time column (e.g.
#'   `"time"`).
#' @param threshold `numeric` Event count threshold used with `op`.
#' @export
timed_count_trigger <- function(col, time_col, op, threshold) {
  .check_string(col, "col")
  .check_string(time_col, "time_col")
  .check_op(op)
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) stop("`threshold` must be a single numeric value.")
  structure(list(type = "timed_count", col = col, time_col = time_col, op = op, threshold = threshold), class = "trigger")
}

#' @rdname trigger_primitives
#' @param fraction `numeric` Sample fraction (0 < fraction <= 1).
#' @param sample_size `numeric` Target sample size.
#' @export
enroll_trigger <- function(fraction, sample_size) {
  if (missing(fraction) || missing(sample_size)) stop("`fraction` and `sample_size` are required.")
  if (!is.numeric(sample_size) || length(sample_size) != 1L || is.na(sample_size)) {
    stop("`sample_size` must be a single numeric value.")
  }
  if (!is.numeric(fraction) || length(fraction) != 1L || is.na(fraction) || fraction <= 0 || fraction > 1) {
    stop("`fraction` must be a single number in (0, 1].")
  }

  # count_trigger's aggregate predicate decides *whether* to fire but, on its
  # own, doesn't drop unenrolled rows (dplyr::filter recycles a length-1
  # logical across all rows). AND it with notna_trigger so the analysis only
  # ever sees enrolled subjects.
  count_trigger("enroll_time", ">=", fraction * sample_size) & notna_trigger("enroll_time")
}

#' @rdname trigger_primitives
#' @param cal_time `numeric` Calendar time at or after which to trigger.
#'   Uses `>=` so it works with both stochastic (continuous) and deterministic
#'   (integer) schedules.
#' @export
calendar_trigger <- function(cal_time) {
  if (missing(cal_time)) stop("`cal_time` is required.")
  if (!is.numeric(cal_time) || length(cal_time) != 1L) {
    stop("`cal_time` must be a single numeric value.")
  }

  value_trigger("time", ">=", cal_time)
}

#' @rdname trigger_primitives
#' @param e1,e2 `trigger` objects to combine.
#' @export
`&.trigger` <- function(e1, e2) {
  stopifnot(inherits(e1, "trigger"), inherits(e2, "trigger"))
  structure(list(predicates = list(e1, e2), combinator = "&"), class = "trigger")
}

#' @rdname trigger_primitives
#' @param e1,e2 `trigger` objects to combine.
#' @export
`|.trigger` <- function(e1, e2) {
  stopifnot(inherits(e1, "trigger"), inherits(e2, "trigger"))
  structure(list(left = e1, right = e2, combinator = "|"), class = "trigger")
}

#' Build a Condition that Fires at a Calendar Time
#'
#' Convenience constructor: wraps `Condition$new()` with a
#' `calendar_trigger()`. The returned `Condition` should be passed to
#' `Trial$new(conditions = list(...))`.
#'
#' @param cal_time `numeric` Calendar time at or after which to fire.
#' @param analysis `function` or `NULL` Called as
#'   `analysis(df, current_time, ...)`. If `NULL`, the filtered snapshot
#'   is returned as-is with a warning.
#' @param name `character` or `NULL` Result key. Defaults to
#'   `"cal_time_<cal_time>"`.
#'
#' @return A [`Condition`] object.
#'
#' @seealso [Condition], [condition_enrollment_fraction()], [Trial],
#'   `calendar_trigger()`.
#'
#' @export
#'
#' @examples
#' cond <- condition_calendar_time(
#'   cal_time = 12,
#'   analysis = function(df, current_time) {
#'     data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
#'   }
#' )
condition_calendar_time <- function(cal_time, analysis = NULL, name = NULL) {
  if (missing(cal_time)) stop("`cal_time` is required.")
  stopifnot(is.numeric(cal_time))
  if (is.null(name)) name <- paste0("cal_time_", paste(cal_time, collapse = "_"))

  Condition$new(where = calendar_trigger(cal_time), analysis = analysis, name = name)
}

#' Build a Condition that Fires at an Enrollment Fraction
#'
#' Convenience constructor: wraps `Condition$new()` with an
#' `enroll_trigger()`. The returned `Condition` should be passed to
#' `Trial$new(conditions = list(...))`.
#'
#' @param fraction `numeric` Sample fraction (0 < fraction <= 1).
#' @param sample_size `integer` Target sample size.
#' @param analysis `function` or `NULL` Called as
#'   `analysis(df, current_time, ...)`. If `NULL`, the filtered snapshot
#'   is returned as-is with a warning.
#' @param name `character` or `NULL` Result key. Defaults to
#'   `"frac_<fraction>"`.
#'
#' @return A [`Condition`] object.
#'
#' @seealso [Condition], [condition_calendar_time()], [Trial],
#'   `enroll_trigger()`.
#'
#' @export
#'
#' @examples
#' cond <- condition_enrollment_fraction(
#'   fraction    = 0.5,
#'   sample_size = 100,
#'   analysis    = function(df, current_time) {
#'     data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
#'   }
#' )
condition_enrollment_fraction <- function(fraction, sample_size, analysis = NULL, name = NULL) {
  if (missing(fraction) || missing(sample_size)) stop("`fraction` and `sample_size` are required.")
  stopifnot(is.numeric(sample_size) && length(sample_size) == 1L)
  stopifnot(fraction > 0 && fraction <= 1)
  if (is.null(name)) name <- paste0("frac_", fraction)

  Condition$new(where = enroll_trigger(fraction, sample_size), analysis = analysis, name = name)
}

#' Trigger Analysis When Enough TTE Events Have Occurred
#'
#' Builds a [`Condition`] that fires once a target number of time-to-event
#' (TTE) endpoint events have been observed by the current trial time.
#' Enrollment guards are always included: only subjects who are enrolled
#' before the current time are eligible.
#'
#' @param event_col `character` Column holding each subject's event calendar
#'   time. `NA` indicates no event has occurred yet.
#' @param n_events `numeric` Number of events required to fire.
#' @param time_col `character` Column holding the current trial time.
#'   Defaults to `"time"` (the column added by `Trial$run()`).
#' @param enroll_col `character` Column holding each subject's enrollment
#'   time. Defaults to `"enroll_time"`.
#' @param analysis `function` or `NULL` Optional analysis function called as
#'   `analysis(df, current_time, ...)`. If `NULL`, the filtered snapshot is
#'   returned as-is with a warning.
#' @param name `character` or `NULL` Result key. Defaults to
#'   `"events_<n_events>"`.
#' @param op `character` Comparison operator for the event-count threshold.
#'   Must be one of `c(">=", "<=", ">", "<", "==", "!=", "%in%")`.
#'   Defaults to `">="`.
#'
#' @return A [`Condition`] object.
#'
#' @seealso [Condition], [condition_calendar_time()], [condition_enrollment_fraction()],
#'   [timed_count_trigger()], [notna_trigger()], [col_trigger()].
#'
#' @export
#'
#' @examples
#' cond <- trigger_by_events(
#'   event_col = "pfs_event_time",
#'   n_events  = 100,
#'   analysis  = function(df, current_time) {
#'     data.frame(n_events = sum(!is.na(df$pfs_event_time)), fired_at = current_time)
#'   }
#' )
trigger_by_events <- function(
    event_col,
    n_events,
    time_col   = "time",
    enroll_col = "enroll_time",
    analysis   = NULL,
    name       = NULL,
    op         = ">="
) {
  if (missing(event_col) || missing(n_events)) stop("`event_col` and `n_events` are required.")
  .check_string(event_col, "event_col")
  if (!is.numeric(n_events) || length(n_events) != 1L || is.na(n_events)) stop("`n_events` must be a single numeric value.")
  .check_op(op)
  if (is.null(name)) name <- paste0("events_", n_events)

  trig <- timed_count_trigger(event_col, time_col, op, n_events) &
    notna_trigger(enroll_col) &
    col_trigger(enroll_col, "<=", time_col)

  Condition$new(where = trig, analysis = analysis, name = name)
}

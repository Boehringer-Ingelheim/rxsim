.trigger_ops <- c(">=", "<=", ">", "<", "==", "!=", "%in%")

#' @name trigger_primitives
#' @title Build Safe Trial Triggers
#' @description Create inert trigger specifications that can be safely
#'   passed to [`Condition`] or composed with `&` and `|`.
#'
#' @param col `character` Column name referenced by the trigger.
#' @param op `character` Comparison operator. Must be one of
#'   `c(">=", "<=", ">", "<", "==", "!=", "%in%")`.
#' @param rhs Right-hand side value. Must be atomic for `value_trigger()` and
#'   numeric for `count_trigger()`.
#'
#' @return An `rxsim_trigger` object.
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
  if (!is.character(col) || length(col) != 1L || is.na(col)) stop("`col` must be a single character string.")
  if (!is.character(op) || length(op) != 1L || is.na(op) || !op %in% .trigger_ops) stop("`op` must be one of: >=, <=, >, <, ==, !=, %in%.")
  if (!is.atomic(rhs)) {
    stop("`rhs` must be atomic for `value_trigger()`.")
  }

  structure(list(type = "value", col = col, op = op, rhs = rhs), class = "rxsim_trigger")
}

#' @rdname trigger_primitives
#' @export
count_trigger <- function(col, op, rhs) {
  if (!is.character(col) || length(col) != 1L || is.na(col)) stop("`col` must be a single character string.")
  if (!is.character(op) || length(op) != 1L || is.na(op) || !op %in% .trigger_ops) stop("`op` must be one of: >=, <=, >, <, ==, !=, %in%.")
  if (!is.atomic(rhs) || !is.numeric(rhs)) {
    stop("`rhs` must be numeric for `count_trigger()`.")
  }

  structure(list(type = "count", col = col, op = op, rhs = rhs), class = "rxsim_trigger")
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

  count_trigger("enroll_time", ">=", fraction * sample_size)
}

#' @rdname trigger_primitives
#' @param cal_time `numeric` Calendar time(s) at which to trigger.
#' @export
calendar_trigger <- function(cal_time) {
  if (missing(cal_time)) stop("`cal_time` is required.")
  if (!is.numeric(cal_time)) {
    stop("`cal_time` must be numeric.")
  }

  value_trigger("time", "%in%", cal_time)
}

#' @rdname trigger_primitives
#' @param e1,e2 `rxsim_trigger` objects to combine.
#' @export
`&.rxsim_trigger` <- function(e1, e2) {
  stopifnot(inherits(e1, "rxsim_trigger"), inherits(e2, "rxsim_trigger"))
  structure(list(predicates = list(e1, e2), combinator = "&"), class = "rxsim_trigger")
}

#' @rdname trigger_primitives
#' @param e1,e2 `rxsim_trigger` objects to combine.
#' @export
`|.rxsim_trigger` <- function(e1, e2) {
  stopifnot(inherits(e1, "rxsim_trigger"), inherits(e2, "rxsim_trigger"))
  structure(list(left = e1, right = e2, combinator = "|"), class = "rxsim_trigger")
}

#' Trigger Analysis at a Calendar Time
#'
#' Builds a [`Condition`] that fires when the trial clock reaches a specified
#' calendar time. The returned `Condition` should be passed to
#' `Trial$new(conditions = list(...))`.
#'
#' @param cal_time `numeric` Calendar time(s) at which to trigger.
#' @param analysis `function` or `NULL` Optional analysis function called as
#'   `analysis(df, current_time, ...)`. If `NULL`, the filtered snapshot
#'   is returned as-is with a warning.
#' @param name `character` or `NULL` Result key. Defaults to
#'   `"cal_time_<cal_time>"`.
#'
#' @return A [`Condition`] object.
#'
#' @seealso [Condition], [trigger_by_fraction()], [Trial],
#'   `calendar_trigger()`.
#'
#' @export
#'
#' @examples
#' cond <- trigger_by_calendar(
#'   cal_time = 12,
#'   analysis = function(df, current_time) {
#'     data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
#'   }
#' )
trigger_by_calendar <- function(cal_time, analysis = NULL, name = NULL) {
  if (missing(cal_time)) stop("`cal_time` is required.")
  stopifnot(is.numeric(cal_time))
  if (is.null(name)) name <- paste0("cal_time_", paste(cal_time, collapse = "_"))

  Condition$new(where = calendar_trigger(cal_time), analysis = analysis, name = name)
}

#' Trigger Analysis at a Sample Fraction
#'
#' Builds a [`Condition`] that fires when a given fraction of the target sample
#' has been enrolled. The returned `Condition` should be passed to
#' `Trial$new(conditions = list(...))`.
#'
#' @param fraction `numeric` Sample fraction (0 < fraction <= 1).
#' @param sample_size `integer` Target sample size.
#' @param analysis `function` or `NULL` Optional analysis function called as
#'   `analysis(df, current_time, ...)`. If `NULL`, the filtered snapshot
#'   is returned as-is with a warning.
#' @param name `character` or `NULL` Result key. Defaults to
#'   `"frac_<fraction>"`.
#'
#' @return A [`Condition`] object.
#'
#' @seealso [Condition], [trigger_by_calendar()], [Trial],
#'   `enroll_trigger()`.
#'
#' @export
#'
#' @examples
#' cond <- trigger_by_fraction(
#'   fraction    = 0.5,
#'   sample_size = 100,
#'   analysis    = function(df, current_time) {
#'     data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
#'   }
#' )
trigger_by_fraction <- function(fraction, sample_size, analysis = NULL, name = NULL) {
  if (missing(fraction) || missing(sample_size)) stop("`fraction` and `sample_size` are required.")
  stopifnot(is.numeric(sample_size) && length(sample_size) == 1L)
  stopifnot(fraction > 0 && fraction <= 1)
  if (is.null(name)) name <- paste0("frac_", fraction)

  Condition$new(where = enroll_trigger(fraction, sample_size), analysis = analysis, name = name)
}

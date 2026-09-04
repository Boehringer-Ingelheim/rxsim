# Build one nested call for the whole trigger tree.
.trigger_expr <- function(trigger) {
  if (is.null(trigger$type)) {
    op    <- if (identical(trigger$combinator, "&")) "&" else "|"
    preds <- if (op == "&") trigger$predicates else list(trigger$left, trigger$right)
    return(Reduce(function(a, b) call(op, a, b), lapply(preds, .trigger_expr)))
  }

  if (identical(trigger$type, "value")) {
    call(trigger$op, call("[[", quote(.data), trigger$col), trigger$rhs)
  } else if (identical(trigger$type, "notna")) {
    call("!", call("is.na", call("[[", quote(.data), trigger$col)))
  } else if (identical(trigger$type, "col_compare")) {
    call(trigger$op, call("[[", quote(.data), trigger$col), call("[[", quote(.data), trigger$ref_col))
  } else if (identical(trigger$type, "timed_count")) {
    col_expr  <- call("[[", quote(.data), trigger$col)
    time_expr <- call("[[", quote(.data), trigger$time_col)
    call(trigger$op,
         call("sum", call("&", call("!", call("is.na", col_expr)), call("<=", col_expr, time_expr))),
         trigger$threshold)
  } else {
    call(trigger$op,
         call("sum", call("!", call("is.na", call("[[", quote(.data), trigger$col)))),
         trigger$rhs)
  }
}

.trigger_to_quosures <- function(trigger) {
  list(rlang::new_quosure(.trigger_expr(trigger), env = rlang::current_env()))
}

#' Condition: Stateful trigger and analysis unit
#'
#' @description
#' A `Condition` encapsulates a single trigger rule that is evaluated against
#' a data snapshot at each simulated timepoint. It combines three concerns:
#'
#' \enumerate{
#'   \item **Filtering**  -  a `dplyr::filter()` expression selects the rows
#'     relevant to this condition (e.g. "only enrolled subjects in arm A").
#'   \item **Analysis**  -  an optional function transforms the filtered snapshot
#'     into a result (e.g. a t-test, a subject count, a Go/No-Go decision).
#'   \item **Trigger bookkeeping**  -  the condition fires only when the
#'     filtered data is non-empty, the cooldown period has elapsed since the
#'     last trigger, and the maximum trigger count has not been reached.
#' }
#'
#' `Condition` objects are stored in `trial$conditions` and evaluated by
#' [`Trial`]`$run()` at each timepoint.
#'
#' @details
#' **Three-gate logic.** A trigger fires only when all three gates pass:
#' \enumerate{
#'   \item The filtered snapshot contains at least one row.
#'   \item `trigger_count < max_triggers`.
#'   \item `current_time - last_trigger_time >= cooldown` (or the condition
#'     has never fired before).
#' }
#' If any gate fails, `check_conditions()` returns an empty list and state
#' is not updated.
#'
#' On a successful trigger, the condition calls
#' `analysis(df, current_time, ...)` and stores the result under
#' `name` (or `1L` when no name is set). Any values in `analysis_args` are
#' appended as additional named arguments. If no analysis function is
#' provided, the filtered data frame is returned as-is with a warning.
#'
#' @seealso
#' \itemize{
#'   \item [`Timer`] for managing trial timepoints
#'   \item [`Trial`] for running the simulation and iterating over conditions
#'   \item [value_trigger()], [count_trigger()], [enroll_trigger()],
#'     [calendar_trigger()] for building trigger specifications
#'   \item [`condition_calendar_time()`] and [`condition_enrollment_fraction()`]
#'     for convenient `Condition` constructors
#' }
#'
#' @examples
#' # Build a snapshot data frame (enroll_time = NA means not yet enrolled)
#' snapshot <- data.frame(
#'   arm         = c("A", "A", "A", "B"),
#'   status      = c("active", "active", "active", "active"),
#'   enroll_time = c(1, 2, 3, NA_real_),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Analysis function: count enrolled subjects and record fire time
#' count_fn <- function(df, current_time) {
#'   data.frame(n_active = nrow(df), fired_at = current_time)
#' }
#'
#' # Condition fires once when 3+ of 4 subjects are enrolled (max_triggers = 1)
#' cond <- Condition$new(
#'   where        = enroll_trigger(fraction = 0.75, sample_size = 4),
#'   analysis     = count_fn,
#'   name         = "interim_A",
#'   cooldown     = 0,
#'   max_triggers = 1L
#' )
#'
#' # First call: fires and returns analysis result
#' res <- cond$check_conditions(snapshot, current_time = 5)
#' res[["interim_A"]]  # data.frame(n_active = 4, fired_at = 5)
#'
#' # Second call: does not fire (max_triggers already reached)
#' res2 <- cond$check_conditions(snapshot, current_time = 6)
#' length(res2)  # 0
#'
#' @importFrom dplyr filter
#' @export
Condition <- R6::R6Class(
  classname = "Condition",

  public = list(
    # --- fields ---
    #' @field where `list` of quosures (`rlang::quos()`) used as `dplyr::filter()`
    #'   predicates, or a `trigger` object (converted automatically). `NULL` or
    #'   empty list passes the full snapshot.
    where = NULL,

    #' @field trigger_spec The original `trigger` object passed to `where`
    #'   (before quosure conversion), or `NULL`. Used by the fixed fast path to
    #'   compute the earliest possible firing time and skip non-firing
    #'   timepoints. `NULL` means "evaluate at every timepoint" (safe fallback).
    trigger_spec = NULL,

    #' @field analysis `function` or `NULL`. Called as
    #'   `analysis(df, current_time, ...)` on a successful trigger, where `...`
    #'   are any values from `analysis_args`.
    analysis = NULL,

    #' @field analysis_args `list` or `NULL`. Named list of extra values injected
    #'   into the analysis function call as additional named arguments.
    analysis_args = NULL,

    #' @field name `character` or `NULL`. Key labelling the result in the output
    #'   list. Falls back to `1L` when `NULL`.
    name = NULL,

    #' @field cooldown `numeric`. Minimum time units between consecutive
    #'   triggers. Default `0`.
    cooldown = 0,

    #' @field max_triggers `integer` or `Inf`. Maximum number of times this
    #'   condition may fire. Default `1L`.
    max_triggers = 1L,

    #' @field trigger_count `integer`. Number of successful triggers so far.
    #'   Initialised to `0L`.
    trigger_count = 0L,

    #' @field last_trigger_time `numeric`. Calendar time of the most recent
    #'   successful trigger. `NA_real_` until first trigger.
    last_trigger_time = NA_real_,

    # --- constructor ---
    #' @description
    #' Create a new `Condition` instance.
    #'
    #' @param where `trigger` object (from `value_trigger()`, `count_trigger()`,
    #'   `enroll_trigger()`, or `calendar_trigger()`), a `list` of quosures
    #'   from `rlang::quos()`, or `NULL` to use the full snapshot.
    #' @param analysis `function` or `NULL`. Called as
    #'   `analysis(df, current_time, ...)` on a successful trigger, where `...`
    #'   are the values from `analysis_args`.
    #' @param analysis_args `list` or `NULL`. Named list of extra arguments
    #'   passed to the analysis function after `df` and `current_time`.
    #' @param name `character` or `NULL`. Result key. Defaults to `1L`.
    #' @param cooldown `numeric`. Minimum time between triggers. Default `0`.
    #' @param max_triggers `integer`. Maximum trigger count. Default `1L`.
    #'   Use `Inf` for unlimited.
    #'
    #' @return A new `Condition` instance.
    initialize = function(
      where         = NULL,
      analysis      = NULL,
      analysis_args = NULL,
      name          = NULL,
      cooldown      = 0,
      max_triggers  = 1L
    ) {
      if (inherits(where, "trigger")) {
        self$trigger_spec <- where
        where <- .trigger_to_quosures(where)
      }
      self$where         <- where
      self$analysis      <- analysis
      self$analysis_args <- analysis_args
      self$name          <- name

      cooldown <- as.numeric(cooldown)
      if (length(cooldown) != 1L || cooldown < 0 || is.na(cooldown)) {
        stop("`cooldown` must be a single non-negative number.")
      }

      if (length(max_triggers) == 1L && is.infinite(max_triggers) && max_triggers > 0) {
        # Inf means unlimited  -  keep as-is
      } else {
        max_triggers <- as.integer(max_triggers)
        if (length(max_triggers) != 1L || is.na(max_triggers) || max_triggers < 0L) {
          stop("`max_triggers` must be a non-negative integer (use Inf for unlimited).")
        }
      }

      self$cooldown     <- cooldown
      self$max_triggers <- max_triggers
    },

    # --- methods ---

    #' @description
    #' Evaluate this condition against a data snapshot.
    #'
    #' Applies the three-gate logic: non-empty filter result, cooldown
    #' elapsed, and trigger count below `max_triggers`. Returns the analysis
    #' result (or filtered data) on a successful trigger, or an empty list
    #' otherwise.
    #'
    #' @param locked_data `data.frame` The trial snapshot at the current time.
    #' @param current_time `numeric` Calendar time of the current timepoint.
    #'
    #' @return Named `list` with one entry (the analysis result) on success,
    #'   or an empty `list` if the condition did not fire.
    check_conditions = function(locked_data, current_time) {
      stopifnot(is.data.frame(locked_data))

      results <- list()

      key <- if (!is.null(self$name) && nzchar(self$name)) self$name else 1L

      # Filter snapshot (dplyr semantics: NA in predicates drops rows)
      df_i <- if (!is.null(self$where) && length(self$where) > 0) {
        dplyr::filter(locked_data, !!!self$where)
      } else {
        locked_data
      }

      # Gate 1: non-empty match
      if (nrow(df_i) == 0L) return(results)

      # Gate 2: hard cap on number of triggers
      if (is.finite(self$max_triggers) && self$trigger_count >= self$max_triggers) {
        return(results)
      }

      # Gate 3: cooldown
      if (is.finite(self$last_trigger_time)) {
        if ((current_time - self$last_trigger_time) < self$cooldown) {
          return(results)
        }
      }

      if (is.function(self$analysis)) {
        results[[key]] <- do.call(self$analysis, c(list(df_i, current_time), self$analysis_args))
      } else {
        results[[key]] <- df_i
        warning(
          sprintf(
            " returning filtered data as is because condition '%s' has no applicable analysis \n",
            key
          ),
          call. = FALSE
        )
      }

      # Update trigger state after a successful trigger
      self$trigger_count     <- self$trigger_count + 1L
      self$last_trigger_time <- current_time

      results
    }

  ) # end public
) # end class

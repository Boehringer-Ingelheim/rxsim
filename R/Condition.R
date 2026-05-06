#' Condition: Stateful trigger and analysis unit
#'
#' @description
#' A `Condition` encapsulates a single trigger rule that is evaluated against
#' a data snapshot at each simulated timepoint. It combines three concerns:
#'
#' \enumerate{
#'   \item **Filtering** — a `dplyr::filter()` expression selects the rows
#'     relevant to this condition (e.g. "only enrolled subjects in arm A").
#'   \item **Analysis** — an optional function transforms the filtered snapshot
#'     into a result (e.g. a t-test, a subject count, a Go/No-Go decision).
#'   \item **Trigger bookkeeping** — the condition fires only when the
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
#'   \item `current_time - last_trigger_time >= cooldown` (or the condition
#'     has never fired before).
#'   \item `trigger_count < max_triggers`.
#' }
#' If any gate fails, `check_conditions()` returns an empty list and state
#' is not updated.
#'
#' On a successful trigger, the condition calls
#' `analysis(filtered_data, current_time)` and stores the result under
#' `name` (or `1L` when no name is set). If no analysis function is
#' provided, the filtered data frame is returned as-is with a warning.
#'
#' @section Fields:
#' \describe{
#'   \item{\code{where}}{`list` of quosures (from `rlang::quos()`) used as
#'     `dplyr::filter()` predicates. Pass `NULL` or an empty list to skip
#'     filtering and pass the full snapshot to the analysis.}
#'   \item{\code{analysis}}{`function` or `NULL`. Called as
#'     `analysis(filtered_data, current_time)` on a successful trigger.
#'     Should return a `data.frame` or named list. If `NULL`, the filtered
#'     data frame is returned with a warning.}
#'   \item{\code{name}}{`character` or `NULL`. Key used to label the result
#'     in the returned list. Falls back to `1L` when `NULL`.}
#'   \item{\code{cooldown}}{`numeric`. Minimum time units that must elapse
#'     between consecutive triggers. Default `0` (no cooldown).}
#'   \item{\code{max_triggers}}{`integer`. Maximum number of times this
#'     condition may fire. Use `Inf` for unlimited. Default `1L`.}
#'   \item{\code{trigger_count}}{`integer`. Number of successful triggers so
#'     far. Initialised to `0L`.}
#'   \item{\code{last_trigger_time}}{`numeric`. Calendar time of the most
#'     recent successful trigger. Initialised to `NA_real_`.}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$new(where, analysis, name, cooldown, max_triggers)}}{
#'     Construct a new `Condition`. All arguments except `where` are
#'     optional. `cooldown` must be a single non-negative number;
#'     `max_triggers` must be a single non-negative integer or `Inf`.}
#'   \item{\code{$check_conditions(locked_data, current_time)}}{
#'     Evaluate the condition against `locked_data` at `current_time`.
#'     Returns a named `list` containing the analysis result (or filtered
#'     data frame) if the condition fires, or an empty `list` otherwise.
#'     On a successful trigger, `trigger_count` is incremented and
#'     `last_trigger_time` is updated.}
#' }
#'
#' @seealso
#' \itemize{
#'   \item [`Timer`] for managing trial timepoints
#'   \item [`Trial`] for running the simulation and iterating over conditions
#'   \item [`trigger_by_calendar()`] and [`trigger_by_fraction()`] for
#'     convenient `Condition` constructors
#'   \item [`dplyr::filter()`] for predicate syntax
#' }
#'
#' @examples
#' # Build a snapshot data frame
#' snapshot <- data.frame(
#'   arm    = c("A", "A", "A", "B"),
#'   status = c("active", "active", "active", "active"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Analysis function: count active subjects per arm
#' count_fn <- function(dat, current_time) {
#'   data.frame(n_active = nrow(dat), fired_at = current_time)
#' }
#'
#' # Condition fires once when arm A has active subjects (max_triggers = 1)
#' cond <- Condition$new(
#'   where        = rlang::quos(arm == "A", status == "active"),
#'   analysis     = count_fn,
#'   name         = "interim_A",
#'   cooldown     = 0,
#'   max_triggers = 1L
#' )
#'
#' # First call: fires and returns analysis result
#' res <- cond$check_conditions(snapshot, current_time = 5)
#' res[["interim_A"]]  # data.frame(n_active = 3, fired_at = 5)
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
    #'   predicates. `NULL` or empty list passes the full snapshot.
    where = NULL,

    #' @field analysis `function` or `NULL`. Called as
    #'   `analysis(filtered_data, current_time)` on a successful trigger.
    analysis = NULL,

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
    #' @param where `list` of quosures (from `rlang::quos()`) used as filter
    #'   predicates. Pass `NULL` or omit to use the full snapshot.
    #' @param analysis `function` or `NULL`. Called as
    #'   `analysis(filtered_data, current_time)` on a successful trigger.
    #' @param name `character` or `NULL`. Result key. Defaults to `1L`.
    #' @param cooldown `numeric`. Minimum time between triggers. Default `0`.
    #' @param max_triggers `integer`. Maximum trigger count. Default `1L`.
    #'   Use `Inf` for unlimited.
    #'
    #' @return A new `Condition` instance.
    initialize = function(
      where        = NULL,
      analysis     = NULL,
      name         = NULL,
      cooldown     = 0,
      max_triggers = 1L
    ) {
      self$where    <- where
      self$analysis <- analysis
      self$name     <- name

      cooldown <- as.numeric(cooldown)
      if (length(cooldown) != 1L || cooldown < 0 || is.na(cooldown)) {
        stop("`cooldown` must be a single non-negative number.")
      }

      if (length(max_triggers) == 1L && is.infinite(max_triggers) && max_triggers > 0) {
        # Inf means unlimited — keep as-is
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
        results[[key]] <- self$analysis(df_i, current_time)
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

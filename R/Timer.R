

#' Timer: Track timed events across arms
#'
#' @description
#' A class to collect and query _timepoints_ - time-based enrollment and
#' dropout events - across trial arms.
#'
#' Use `add_schedule()` to register events and
#' `get_end_timepoint()` / `get_n_arms()` / `get_unique_times()` for
#' summary queries. The full event table is the public `timelist` field.
#'
#' @details
#' Trigger conditions (filtering + analysis) are now managed by the separate
#' [`Condition`] class. `Condition` objects are stored in `trial$conditions`
#' and evaluated by [`Trial`]`$run()` at each timepoint.
#'
#' Helper functions [`condition_calendar_time()`] and [`condition_enrollment_fraction()`]
#' provide convenient shortcuts for building `Condition` objects; both return
#' a [`Condition`] that you pass to `Trial$new(conditions = list(...))`.
#'
#' @seealso [`Trial`] to coordinate simulations with populations,
#'   [`Condition`] for trigger/analysis logic,
#'   [stochastic_schedule()] / [deterministic_schedule()] to build a schedule.
#'
#' @examples
#' # Basic construction
#' t <- Timer$new(name = "Timer")
#'
#' # Add timepoints
#' t$add_schedule(data.frame(
#'   time   = c(1, 2, 1),
#'   arm    = c("A", "A", "B"),
#'   drop   = c(2L, 1L, 0L),
#'   enroll = c(10L, 12L, 8L)
#' ))
#'
#' # Query
#' t$get_end_timepoint() # max time => 2
#' t$get_n_arms()        # unique arms => 2
#' t$get_unique_times()  # unique times => c(1, 2)
#' t$timelist            # the full event table
#'
#' @importFrom rlang enquos
#' @importFrom dplyr filter
#' @export
Timer <- R6::R6Class(
  classname = "Timer",
  public = list(
    #' @field name `character` Unique identifier for the `Timer` instance.
    name = NULL,

    #' @field timelist `data.frame` A data.frame of timepoints with columns:
    #' - `time` `numeric` Calendar time
    #' - `arm` `character` Unique identifier of the arm
    #' - `drop` `integer` # of subjects dropped at `time`
    #' - `enroll` `integer` # of subjects enrolled at `time`
    timelist = NULL,

    #' @description
    #' Create a new `Timer` instance.
    #'
    #' @param name `character` Unique identifier.
    #' @param timelist `data.frame` Optional data.frame of timepoints with columns
    #'   `time`, `arm`, `drop`, `enroll`. If `NULL`, an empty frame is created.
    #'
    #' @return A new `Timer` instance.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    initialize = function(name, timelist = NULL) {
      stopifnot(is.character(name))
      self$name <- name
      if (is.null(timelist)) {
        self$timelist <- data.frame(
          time   = numeric(0),
          arm    = character(0),
          drop   = integer(0),
          enroll = integer(0),
          stringsAsFactors = FALSE
        )
      } else {
        stopifnot(is.data.frame(timelist))
        self$timelist <- timelist
      }
    },

    #' @description
    #' Add a schedule of timepoints to the timer.
    #'
    #' @param schedule `data.frame` with columns `time` (numeric), `arm`
    #'   (character), `enroll` (integer), `drop` (integer). One row per event;
    #'   a single event is a one-row data frame. Typically the output of
    #'   [stochastic_schedule()] or [deterministic_schedule()].
    #'
    #'   `enroll` and `drop` are subject counts and must be integer
    #'   (`3L`, not `3`) - fractional counts are silently truncated
    #'   downstream, so they are rejected here.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #'
    #' # single event
    #' t$add_schedule(data.frame(time = 1, arm = "A", drop = 1L, enroll = 3L))
    #'
    #' # whole schedule (data.frame() recycles the constant columns)
    #' t$add_schedule(data.frame(time = 2:3, arm = "A", enroll = 2L, drop = 0L))
    add_schedule = function(schedule) {
      if (missing(schedule) || !is.data.frame(schedule)) {
        stop("`schedule` must be a data.frame with columns: time, arm, enroll, drop")
      }
      required_cols <- c("time", "arm", "enroll", "drop")
      missing_cols <- setdiff(required_cols, names(schedule))
      if (length(missing_cols) > 0L) {
        stop(sprintf("Missing required columns in schedule: %s", paste(missing_cols, collapse = ", ")))
      }
      for (col in c("enroll", "drop")) {
        if (!is.integer(schedule[[col]])) {
          stop(sprintf(
            "`%s` must be an integer vector (e.g. 3L), not %s. Subject counts cannot be fractional.",
            col, class(schedule[[col]])[1]
          ))
        }
      }

      # ponytail: rebuild as a plain data.frame so tibble/grouped_df input
      # doesn't leak its class into timelist.
      self$timelist <- rbind(
        self$timelist,
        data.frame(
          time   = schedule$time,
          arm    = schedule$arm,
          drop   = schedule$drop,
          enroll = schedule$enroll,
          stringsAsFactors = FALSE
        )
      )
      invisible(self)
    },

    #' @description
    #' Determine the last timepoint for a given instance of `Timer` class.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #' t$add_schedule(data.frame(time = 3.14, arm = "A", drop = 7L, enroll = 22L))
    #' t$get_end_timepoint()
    get_end_timepoint = function() {
      if (nrow(self$timelist) == 0L) stop("`timelist` is empty.")
      max(self$timelist$time)
    },

    #' @description
    #' Get number of unique arms.
    #'
    #' @return `integer` Number of unique arms.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #' t$add_schedule(data.frame(time = 3.14, arm = "A", drop = 7L, enroll = 22L))
    #' t$add_schedule(data.frame(time = 3.28, arm = "B", drop = 6L, enroll = 23L))
    #' t$get_n_arms()
    get_n_arms = function() length(unique(self$timelist$arm)),

    #' @description
    #' Get unique timepoints.
    #'
    #' @return `numeric` vector of unique times.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #' t$add_schedule(data.frame(time = 3.14, arm = "A", drop = 7L, enroll = 22L))
    #' t$add_schedule(data.frame(time = 3.28, arm = "B", drop = 6L, enroll = 23L))
    #' t$get_unique_times()
    get_unique_times = function() unique(self$timelist$time)
  )
)

#' Timer: Track timed events across arms
#'
#' @description
#' A class to collect and query _timepoints_ - time-based enrollment and
#' dropout events - across trial arms.
#'
#' Use `add_timepoint()` to register events, `get_timepoint()` for lookup,
#' `get_end_timepoint()` / `get_n_arms()` / `get_unique_times()` for
#' summary queries.
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
#'   [`add_timepoints()`] to attach multiple timepoints.
#'
#' @examples
#' # Basic construction
#' t <- Timer$new(name = "Timer")
#'
#' # Add timepoints
#' t$add_timepoint(time = 1, arm = "A", drop = 2L, enroll = 10L)
#' t$add_timepoint(time = 2, arm = "A", drop = 1L, enroll = 12L)
#' t$add_timepoint(time = 1, arm = "B", drop = 0L, enroll = 8L)
#'
#' # Query
#' t$get_end_timepoint() # max time => 2
#' t$get_n_arms()        # unique arms => 2
#' t$get_unique_times()  # unique times => c(1, 2)
#' t$get_timepoint("A", 1) # returns a single-row data.frame
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
    #' Add a timepoint to the timer.
    #'
    #' @param time `numeric` Calendar time.
    #' @param arm `character` Arm identifier.
    #' @param drop `integer` Count of subjects to drop.
    #' @param enroll `integer` Count of subjects to enroll.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #' t$add_timepoint(time = 1, arm = "A", drop = 1L, enroll = 3L)
    add_timepoint = function(time, arm, drop, enroll) {
      stopifnot(is.integer(drop), is.integer(enroll))
      self$timelist <- rbind(
        self$timelist,
        data.frame(
          time = time,
          arm = arm,
          drop = drop,
          enroll = enroll,
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
    #' t$add_timepoint(time = 3.14, arm = "A", drop = 7L, enroll = 22L)
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
    #' t$add_timepoint(time = 3.14, arm = "A", drop = 7L, enroll = 22L)
    #' t$add_timepoint(time = 3.28, arm = "B", drop = 6L, enroll = 23L)
    #' t$get_n_arms()
    get_n_arms = function() length(unique(self$timelist$arm)),

    #' @description
    #' Get unique timepoints.
    #'
    #' @return `numeric` vector of unique times.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #' t$add_timepoint(time = 3.14, arm = "A", drop = 7L, enroll = 22L)
    #' t$add_timepoint(time = 3.28, arm = "B", drop = 6L, enroll = 23L)
    #' t$get_unique_times()
    get_unique_times = function() unique(self$timelist$time),

    #' @description
    #' Get a timepoint by arm and time value.
    #'
    #' @param arm `character` Arm identifier.
    #' @param i `numeric` Time value to look up.
    #'
    #' @return Single-row `data.frame` or `NULL` if not found.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #' t$add_timepoint(time = 3.14, arm = "A", drop = 7L, enroll = 22L)
    #' t$add_timepoint(time = 3.28, arm = "B", drop = 6L, enroll = 23L)
    #' t$get_timepoint("A", 3.14)
    get_timepoint = function(arm, i) {
      if (missing(arm)) stop("`arm` is required.")
      if (missing(i)) stop("`i` is required.")

      idx <- which(self$timelist$time == i & self$timelist$arm == arm)

      if (length(idx) == 0L) return(NULL)
      if (length(idx) > 1L) {
        stop(sprintf("Multiple timepoints found for arm = %s and time = %s.", arm, as.character(i)))
      }

      self$timelist[idx, , drop = FALSE]
    }
  )
)

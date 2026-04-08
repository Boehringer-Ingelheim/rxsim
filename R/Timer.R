#' Timer: Track timed events and apply condition-triggered analyses
#'
#' @description
#' A class to collect and query _timepoints_, time-based events, across arms.
#' Timer class also supports conditions that filter data using [dplyr::filter()]
#' and apply custom analyses.
#'
#' Use `add_timepoint()` to append timepoints, `get_timepoint()` for a lookup,
#' and `check_conditions()` to filter a data frame based on a trigger condition
#' and return either analysis results or the filtered data.
#'
#' @details
#' Helper functions [trigger_by_calendar()] and [trigger_by_fraction()] provide
#' convenient shortcuts for common trigger patterns.
#'
#' @seealso [Trial] to coordinate simulations with populations, [add_timepoints()]
#'   to attach multiple timepoints, [dplyr::filter()] for condition syntax.
#'
#' @examples
#' # Basic construction
#' t <- Timer$new(name = "Timer")
#'
#' # Add timepoints
#' t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 10L)
#' t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 12L)
#' t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 8L)
#'
#' # Query
#' t$get_end_timepoint() # max time => 2
#' t$get_n_arms() # unique arms => 2
#' t$get_unique_times() # unique times => c(1, 2)
#' t$get_timepoint("A", 1) # returns a single timepoint
#'
#' # Add conditions using trigger helpers or dplyr style
#' # Suppose you have a data.frame:
#' df <- data.frame(
#'   id = 1:6,
#'   arm = c("A", "A", "B", "B", "A", "B"),
#'   status = c("active", "inactive", "active", "active", "inactive", "active"),
#'   visit = c(1, 2, 1, 3, 3, 2)
#' )
#'
#' # Analysis function: count rows at/after a given visit, per arm
#' my_analysis <- function(dat, current_time) {
#'   out <- aggregate(id ~ arm, dat, length)
#'   out$current_time <- current_time
#'   out
#' }
#'

#'
#' @importFrom rlang enquos
#' @importFrom dplyr filter
#' @export
Timer <- R6::R6Class(
  classname = "Timer",
  public = list(
    # --- fields ---
    #' @field name `character` Unique identifier for the `Timer` instance.
    name = NULL,

    #' @field timelist `list` A list of timepoints. Each timepoint is a list with keys:
    #' - `time` `numeric` Calendar time
    #' - `arm` `character` Unique identifier of the arm
    #' - `dropper` `integer` # of subjects dropper at `time`
    #' - `enroller` `integer` # of subjects enrolled at `time`
    timelist = NULL,



    # --- constructor ---
    #' @description
    #' Create a new `Timer` instance.
    #'
    #' @param name `character` Unique identifier.
    #' @param timelist `list` Optional list of timepoints.
    #' @param conditions `list` Optional list of condition entries.
    #'
    #' @return A new `Timer` instance.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    initialize = function(
      name,
      timelist = NULL#,
      #conditions = NULL
    ) {
      stopifnot(is.character(name))
      self$name <- name
      self$timelist <- if (is.null(timelist)) list() else timelist
     # self$conditions <- if (is.null(conditions)) list() else conditions
    },

    # --- methods ---
    #' @description
    #' Add a timepoint to a timer.
    #'
    #' @param time `numeric` Calendar time.
    #' @param arm `character` Arm identifier.
    #' @param dropper `integer` Count of subjects to drop.
    #' @param enroller `integer` Count of subjects to enroll.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #' t$add_timepoint(
    #'   time = 1,
    #'   arm = "A",
    #'   dropper = 1L,
    #'   enroller = 3L
    #' )
    add_timepoint = function(time, arm, dropper, enroller) {
      stopifnot(is.integer(dropper), is.integer(enroller))
      tp <- list(time = time, arm = arm, dropper = dropper, enroller = enroller)
      self$timelist <- append(self$timelist, list(tp))
      invisible(self)
    },

    #' @description
    #' Determine the last timepoint for a given instance of `Timer` class.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #' t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
    #' t$get_end_timepoint()
    get_end_timepoint = function() {
      max(sapply(self$timelist, function(x) {
        x$time
      }))
    },

    #' @description
    #' Get number of unique arms.
    #'
    #' @return `integer` Number of unique arms.
    #'
    #' @examples
     #' t <- Timer$new(name = "Timer")
    #' t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
    #' t$add_timepoint(time = 3.28, arm = "B", dropper = 6L, enroller = 23L)
    #' t$get_n_arms()
    get_n_arms = function() length(unique(sapply(self$timelist, function(x) x$arm))),

    #' @description
    #' Get unique timepoints.
    #'
    #' @return `numeric` vector of unique times.
    #'
    #' @examples
     #' t <- Timer$new(name = "Timer")
    #' t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
    #' t$add_timepoint(time = 3.28, arm = "B", dropper = 6L, enroller = 23L)
    #' t$get_unique_times()
    get_unique_times = function() unique(sapply(self$timelist, function(x) x$time)),

    #' @description
    #' Get a timepoint by arm and index.
    #'
    #' @param arm `character` Arm identifier.
    #' @param i `integer` Timepoint index.
    #'
     #' @return `list` timepoint or `NULL` if not found.
    #'
    #' @examples
    #' t <- Timer$new(name = "Timer")
    #' t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
    #' t$add_timepoint(time = 3.28, arm = "B", dropper = 6L, enroller = 23L)
    #'
    #' t$get_timepoint("A", 1)
    get_timepoint = function(arm, i) {
      # Basic validation
      if (missing(arm)) stop("`arm` is required.")
      if (missing(i)) stop("`i` is required.")

      # Extract columns from list with type safety
      times <- vapply(self$timelist, function(x) x$time, FUN.VALUE = numeric(1))
      arms <- vapply(self$timelist, function(x) x$arm, FUN.VALUE = character(1))

      # Find matching indices
      idx <- which(times == i & arms == arm)

      # Handle match outcomes
      if (length(idx) == 0L) {
        return(NULL)
      }
      if (length(idx) > 1L) {
        stop(sprintf("Multiple timepoints found for arm = %s and time = %s.", arm, as.character(i)))
      }

      self$timelist[[idx]]
    }

  ) # end public
) # end class

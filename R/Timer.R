# --- Fix run_readers so it passes current_time into cond$func ---
Timer <- R6::R6Class(
  classname = "Timer",
  public = list(
    # fields
    name = NULL,
    timelist = NULL,
    conditions = NULL,

    # constructor
    initialize = function(
    name,
    timelist = NULL,
    conditions = NULL
    ) {
      stopifnot(is.character(name))
      self$name <- name
      self$timelist <- if (is.null(timelist)) list() else timelist
      self$conditions <- if (is.null(conditions)) list() else conditions
    },

    # methods
    add_timepoint = function(dropper, enroller) {
      stopifnot(is.numeric(dropper), is.numeric(enroller))
      tp <- list(dropper = dropper, enroller = enroller)
      self$timelist <- append(self$timelist, list(tp))
    },

    add_condition = function(
    time = NULL,
    n_events = NULL,
    threshold = NULL,
    func = NULL
    ) {
      cond <- list(time = time, n_events = n_events, threshold = threshold, func = func)
      self$conditions <- append(self$conditions, list(cond))
    },

    get_n_timepoints = function() length(self$timelist),

    get_timepoint = function(i) {
      if (i <= 0 || i > length(self$timelist)) stop("Index out of range")
      self$timelist[[i]]
    },

    check_conditions = function(
    locked_data,
    current_time,
    n_events = NULL
    ) {
      results <- list()
      for (i in seq_along(self$conditions)) {
        cond <- self$conditions[[i]]
        time_ok <- is.null(cond$time) || cond$time == current_time
        events_ok <- is.null(cond$n_events) || (!is.null(n_events) && n_events >= cond$n_events)
        threshold_ok <- is.null(cond$threshold) || (is.numeric(cond$threshold) && n_events >= cond$threshold)

        if (time_ok && events_ok && threshold_ok) {
          if (is.function(cond$func)) {
            # Pass both locked_data and current_time into the function
            results[[paste0("", i)]] <- cond$func(locked_data, current_time)
          } else {
            results[[paste0("", i)]] <- locked_data
          }
        }
      }
      results
    }

  ) # end public
) # end class

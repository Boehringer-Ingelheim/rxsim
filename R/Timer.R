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
    n_enrolled = NULL,
    threshold = NULL,
    analysis = NULL # this is the action function
    ) {
      cond <- list(time = time, n_enrolled = n_enrolled, threshold = threshold, analysis = analysis)
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
    n_enrolled = NULL
    ) {
      results <- list()
      for (i in seq_along(self$conditions)) {
        cond <- self$conditions[[i]]
        time_ok <- is.null(cond$time) || cond$time == current_time
        enrolled_ok <- is.null(cond$n_enrolled) || (!is.null(n_enrolled) && n_enrolled >= cond$n_enrolled)
        threshold_ok <- is.null(cond$threshold) || (is.numeric(cond$threshold) && n_enrolled >= cond$threshold)

        if (time_ok && enrolled_ok && threshold_ok) {
          if (is.function(cond$analysis)) {
            # Pass both locked_data and current_time into the function
            results[[paste0("",i)]] <- cond$analysis(locked_data, current_time)
          } else {
            results[[paste0("",i)]] <- locked_data
          }
        }
      }
      results
    }

  ) # end public
) # end class

Trial <- R6::R6Class(
  classname = "Trial",
  public = list(
    # fields
    name = NULL,
    description = NULL,
    silent = FALSE,
    seed = NULL,
    timer = NULL,        # a Timers object
    population = NULL,   # now a list of Population objects
    locked_data = NULL,
    results = NULL,

    # constructor
    initialize = function(
    name,
    description = name,
    silent = FALSE,
    seed = NULL,
    timer = NULL,
    population = list(),   # default empty list
    locked_data = list(),
    results = list(),
    ...
    ) {
      stopifnot(is.character(name))
      self$name <- name
      self$description <- description
      self$silent <- silent
      self$seed <- seed
      if (!is.null(seed)) set.seed(seed)

      self$timer <- timer
      stopifnot(is.list(population))   # enforce list of Population objects
      self$population <- population
      self$locked_data <- locked_data
      self$results <- results
    },

    # methods
    # run method: loop through timelist and apply conditions
    run = function() {
      if (is.null(self$timer) || length(self$population) == 0) {
        stop("Timer and population list must be set before running run()")
      }

      n_timepoints <- self$timer$get_n_timepoints()

      for (i in seq_len(n_timepoints)) {
        tp <- self$timer$get_timepoint(i)

        # apply enrollment/dropout to each Population object in the list
        for (p in self$population) {
          if (!is.null(tp$enroller)) {
            p$set_enrolled(tp$enroller, time = i)
          }
          if (!is.null(tp$dropper)) {
            p$set_dropped(tp$dropper, time = i)
          }
        }

        # collect snapshots from all populations
        locked_snapshot <- lapply(self$population, function(p) p$get_data())
        n_enrolled <- sum(sapply(self$population, function(p) sum(!is.na(p$enrolled)))) # this right here! is enrollment

        # check all conditions
        results <- self$timer$check_conditions(
          locked_data = locked_snapshot,
          current_time = i,
          n_enrolled = n_enrolled
        )

        if (length(results) > 0) {
          self$locked_data[[paste0("time_", i)]] <- locked_snapshot
          self$results[[paste0("time_", i)]] <- results
        }
      }
    }

  ) # end public
) # end class

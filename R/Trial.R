Trial <- R6::R6Class(
  classname = "Trial",
  public = list(
    # fields
    name = NULL,
    seed = NULL,
    timer = NULL,        # a Timers object
    population = NULL,   # list of Population objects
    locked_data = NULL,  # list of snapshots per time
    results = NULL,      # list of results per time

    # constructor
    initialize = function(
    name,
    seed = NULL,
    timer = NULL,
    population = list(),   # default empty list
    locked_data = list(),
    results = list(),
    ...
    ) {
      stopifnot(is.character(name))
      self$name <- name
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

      if( self$timer$get_n_arms() != length(self$population))
      {
        stop("Need timers for the same amount of arms run()")

      }

        for (i in  sort(self$timer$get_unique_times())){



        # apply enrollment/dropout to each Population object in the list
        for (p in self$population) {
          #add an error statement if the time/population dne
            tp <- self$timer$get_timepoint(p$name,i)
          if(length(tp)>1){
          if (!is.null(tp$enroller)) {
            p$set_enrolled(tp$enroller, time = i)
          }
          if (!is.null(tp$dropper)) {
            p$set_dropped(tp$dropper, time = i)
          }
          }
        }

        # Collect raw snapshots from all populations (as a list)
        locked_snapshot_list <- lapply(self$population, function(p) {
          subset(cbind(p$data, data.frame(
            #add measurement_time column
            enroll_time = rep(x=p$enrolled,times=dim(p$data)[1]),
            drop_time = rep(p$dropped,times=dim(p$data)[1])
          )), !is.na(p$enrolled))
        })


        combined <- do.call(rbind, locked_snapshot_list)


        # Add current time column for predicates like time >= 1
        combined$time <- i

        # Check all conditions on the combined snapshot
        results <- self$timer$check_conditions(
          locked_data  = combined,
          current_time = i
        )

        # Store only if there are results
        if (length(results) > 0) {
          self$locked_data[[paste0("time_", i)]] <- combined
          self$results[[paste0("time_", i)]] <- results
        }
      }
        }

  ) # end public
) # end class

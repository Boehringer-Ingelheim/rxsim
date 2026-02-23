#' Trial: Simulate a multi‑arm clinical trial
#'
#' @description
#' The `Trial` class coordinates one or more `Population` objects and a `Timer`
#' to simulate a clinical trial.
#'
#' At each unique time defined in the trial's `Timer`, the `Trial`:
#'
#' - applies enrollment and dropout updates to each `Population`
#' - builds a snapshot of all currently enrolled subjects
#' - evaluates all conditions in the `Timer`
#' - stores both the snapshot `locked_data` and the analysis outputs `results`
#'
#' Use `run()` to execute the simulation.
#'
#' @examples
#' # Create two populations
#' popA <- Population$new("A", data = vector_to_dataframe(rnorm(10)))
#' popB <- Population$new("B", data = vector_to_dataframe(rnorm(12)))
#'
#' # Create a timer and add timepoints
#' t <- Timer$new("Timer")
#' t$add_timepoint(time = 1, arm = "A", dropper = 0L, enroller = 4L)
#' t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 5L)
#' t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 2L)
#' t$add_timepoint(time = 2, arm = "B", dropper = 2L, enroller = 3L)
#'
#' # Create a trial
#' trial <- Trial$new(
#'   name = "ExampleTrial",
#'   seed = 123,
#'   timer = t,
#'   population = list(popA, popB)
#' )
#'
#' # Run the simulation
#' trial$run()
#'
#' prettify_results(trial$results)
#'
#' @export

Trial <- R6::R6Class(
  classname = "Trial",
  public = list(
    # --- fields ---


    #' @field name `character` Unique identifier for the trial
    name = NULL,


    #' @field seed `numeric` Optional seed for reproducibility. If provided,
    #' `set.seed()` is called during initialization.
    seed = NULL,


    #' @field timer `Timer` A `Timer` object describing timepoints and conditions.
    timer = NULL,        # a Timers object

    #' @field population `list` A list of `Population` objects, one per arm.
    population = NULL,   # list of Population objects

    #' @field locked_data `list` Snapshots of  subject‑level data at each timepoint.
    locked_data = NULL,  # list of snapshots per time

    #' @field results (`list`) Analysis outputs for each condition.
    results = NULL,      # list of results per time

    # --- constructorc ---
    #' @description
    #' Create a new `Trial` instance.
    #'
    #' @param name `character` Unique identifier for the trial.
    #' @param seed `numeric` or `NULL` Optional random seed. If not `NULL`,
    #' sets the RNG seed for reproducibility.
    #' @param timer `Timer` A `Timer` object defining timepoints and conditions.
    #' @param population `list` A list of `Population` objects, one for each arm.
    #' @param locked_data `list` Generated at each Trial$run() call.
    #' @param results (`list`) Generated at each Trial$run() call.
    #'
    #' @returns A new `Trial` instance.
    #'
    #' @examples
    #' t <- Timer$new(name="simple_timer")
    #' pop <- Population$new(
    #'   name = "simple_pop",
    #'   data = vector_to_dataframe(rnorm(5))
    #' )
    #' pop$set_enrolled(5, 1)
    #' Trial$new(name = "simple_trial", timer=t, population = list(pop))
    initialize = function(
    name,
    seed = NULL,
    timer = NULL,
    population = list(),   # default empty list
    locked_data = list(),
    results = list()
    ) {
      stopifnot(is.character(name))
      self$name <- name
      self$seed <- seed
      if (!is.null(seed)) set.seed(seed)


      stopifnot(is.list(population))   # enforce list of Population objects

      if (length(timer$timelist) == 0){
        if (all(sapply(population, function(x) all(is.na(x$enrolled))))) {
          stop("Neither Timer nor Population has enrollment data.")
        }
        else {
          timepoints <- data.frame(
            time = unlist(lapply(population, function(x) x$enrolled), recursive = FALSE),
            arm = rep(sapply(population, function(x) x$name), sapply(population, function(x) x$n)),
            enroller = 1L,
            dropper = 0L
          )
          add_timepoints(timer, timepoints)
          self$timer <- timer
        }
      } else self$timer <- timer

      self$population <- population
      self$locked_data <- locked_data
      self$results <- results

      # routine for empty time list

    },

    # --- methods ---

    #' @description
    #' Execute a trial simulation.
    #'
    #' At each unique time defined by the trial's `Timer`:
    #' - Apply enrollment and dropout actions to each `Population`
    #' - Build a combined snapshot of all currently enrolled subjects
    #' - Attach a `time` column to the snapshot
    #' - Evaluate all condition readers via `Timer$check_conditions()`
    #' - Store snapshots and condition outputs under time‑indexed list keys
    #'
    #' @returns Updates `locked_data` and `results` fields.
    #'
    #' @examples
    #' # Create two populations
    #' popA <- Population$new("A", data = vector_to_dataframe(rnorm(10)))
    #' popB <- Population$new("B", data = vector_to_dataframe(rnorm(12)))
    #'
    #' # Create a timer and add timepoints
    #' t <- Timer$new("Timer")
    #' t$add_timepoint(time = 1, arm = "A", dropper = 0L, enroller = 4L)
    #' t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 5L)
    #' t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 2L)
    #' t$add_timepoint(time = 2, arm = "B", dropper = 2L, enroller = 3L)
    #'
    #' # Create a trial
    #' trial <- Trial$new(
    #'   name = "ExampleTrial",
    #'   seed = 123,
    #'   timer = t,
    #'   population = list(popA, popB)
    #' )
    #'
    #' # Run the simulation
    #' trial$run()
    #'
    #' prettify_results(trial$results)
    run = function() {
      if (is.null(self$timer) || length(self$population) == 0) {
        stop("Timer and population list must be set before running run()")
      }

      # if( self$timer$get_n_arms() != length(self$population))
      # {
      #   stop("Need timers for the same amount of arms run()")
      #
      # }

        for (i in  sort(self$timer$get_unique_times())){

        # apply enrollment/dropout to each Population object in the list
        for (p in self$population) {
          #add an error statement if the time/population dne
          tp <- self$timer$get_timepoint(p$name, i)
          if (!is.null(tp)) {
            # check whether any one left to enroll?
            if (
              !is.null(tp$enroller) &
              any(is.na(p$enrolled))
              ) p$set_enrolled(tp$enroller, time = i)
            # check whether anyone left to drop
            if (
              !is.null(tp$dropper) &
              any(!is.na(p$enrolled)) &
              any(is.na(p$dropped))
              )p$set_dropped(tp$dropper,  time = i)
          }
        }

        # Collect raw snapshots from all populations (as a list)
          locked_snapshot_list <- lapply(self$population, function(p) {
            keep <- !is.na(p$enrolled)
            cbind(p$data[keep, , drop = FALSE],
            enroll_time = rep(x=p$enrolled[keep],times=p$m),
            drop_time   = rep(x=p$dropped[keep],times=p$m) )

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

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
#' popA <- Population$new("A", data = rnorm(10))
#' popB <- Population$new("B", data = rnorm(12))
#'
#' # Create a timer and add timepoints
#' t <- Timer$new("Timer")
#' t$add_timepoint(time = 1, arm = "A", dropper = 0, enroller = 4)
#' t$add_timepoint(time = 1, arm = "B", dropper = 0, enroller = 5)
#' t$add_timepoint(time = 2, arm = "A", dropper = 1, enroller = 2)
#' t$add_timepoint(time = 2, arm = "B", dropper = 2, enroller = 3)
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
    #' Trial$new(name = "Simple", population = list())
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

      self$timer <- timer
      stopifnot(is.list(population))   # enforce list of Population objects
      self$population <- population
      self$locked_data <- locked_data
      self$results <- results
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
    #' trial <- Trial$new(name = "Demo", timer = t, population = list(popA, popB))
    #' trial$run()
    #
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
          tp <- self$timer$get_timepoint(p$name, i)
          if (!is.null(tp)) {
            if (!is.null(tp$enroller)) p$set_enrolled(tp$enroller, time = i)
            if (!is.null(tp$dropper))  p$set_dropped(tp$dropper,  time = i)
          }
        }

        # Collect raw snapshots from all populations (as a list)
          locked_snapshot_list <- lapply(self$population, function(p) {
            keep <- !is.na(p$enrolled)
            cbind(p$data[keep, , drop = FALSE],
                  enroll_time = p$enrolled[keep],
                  drop_time   = p$dropped[keep])
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

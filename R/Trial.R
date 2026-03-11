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
#' Use `run()` to execute the simulation. Trigger conditions are best added with
#' helper functions [trigger_by_calendar()] or [trigger_by_fraction()].
#'
#' @seealso [Population], [Timer], [prettify_results()], [replicate_trial()], [clone_trial()].
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
#' # Add trigger for final analysis at time 2
#' trigger_by_calendar(2, t, analysis = function(df, current_time) {
#'   nrow(df)
#' })
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


    #' @field name `character` Unique trial identifier.
    name = NULL,

    #' @field seed `numeric` or `NULL` Random seed for reproducibility.
    seed = NULL,

    #' @field timer `Timer` object with timepoints and conditions.
    timer = NULL,

    #' @field population `list` of [Population] objects, one per arm.
    population = NULL,

    #' @field locked_data `list` Snapshots at each timepoint.
    locked_data = NULL,

    #' @field results `list` Analysis outputs per condition.
    results = NULL,

    # --- constructor ---
    #' @description
    #' Create a new `Trial` instance.
    #'
    #' @param name `character` Unique identifier for the trial.
    #' @param seed `numeric` or `NULL` Optional random seed for reproducibility.
    #' @param timer `Timer` object defining timepoints and conditions.
    #' @param population `list` of [Population] objects, one per arm.
    #' @param locked_data `list` Generated at each `$run()` call.
    #' @param results `list` Analysis outputs generated at each `$run()` call.
    #'
    #' @return A new `Trial` instance.
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
      population = list(), # default empty list
      locked_data = list(),
      results = list()
    ) {
      stopifnot(is.character(name))
      self$name <- name
      self$seed <- seed
      if (!is.null(seed)) set.seed(seed)
      if (!is.null(timer) && !inherits(timer, "Timer")) stop("`timer` must be a Timer instance.")
      stopifnot(is.list(population))

      if (is.null(timer) || length(timer$timelist) == 0) {
        # If timer has no timepoints, extract from population enrollment times
        if (all(sapply(population, function(x) all(is.na(x$enrolled))))) {
          stop("Neither Timer nor Population has enrollment data.")
        } else {
          if (is.null(timer)) timer <- Timer$new(name = paste0(name, "_timer"))
          timepoints <- data.frame(
            time = unlist(lapply(population, function(x) x$enrolled), recursive = FALSE),
            arm = rep(sapply(population, function(x) x$name), sapply(population, function(x) x$n)),
            enroller = 1L,
            dropper = 0L
          )
          add_timepoints(timer, timepoints)
          self$timer <- timer
        }
      } else {
        self$timer <- timer
      }

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
    #' @return Updates `locked_data` and `results` fields.
    #'
    #' @seealso [Timer], [prettify_results()].
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

      plan_df <- dplyr::bind_rows(self$timer$timelist)
      if (nrow(plan_df) == 0L) {
        return(invisible(self))
      }

      # if( self$timer$get_n_arms() != length(self$population))
      # {
      #   stop("Need timers for the same amount of arms run()")
      #
      # }

      for (i in sort(unique(plan_df$time))) {
        # Apply enrollment/dropout updates to each population at this timepoint
        for (p in self$population) {
          idx <- which(plan_df$arm == p$name & plan_df$time == i)
          if (length(idx) > 0L) {
            enroller_n <- as.integer(sum(plan_df$enroller[idx], na.rm = TRUE))
            dropper_n <- as.integer(sum(plan_df$dropper[idx], na.rm = TRUE))

            if (enroller_n > 0L && sum(is.na(p$enrolled)) > 0L) {
              p$set_enrolled(enroller_n, time = i)
            }

            if (dropper_n > 0L && sum(is.na(p$dropped) & !is.na(p$enrolled)) > 0L) {
              p$set_dropped(dropper_n, time = i)
            }
          }
        }

        # Create snapshots of enrolled subjects from all populations
        locked_snapshot_list <- lapply(self$population, function(p) {
          keep <- !is.na(p$enrolled)
          cbind(
            p$data[keep, , drop = FALSE],
            enroll_time = rep(x=p$enrolled[keep],times=p$n_readouts),
            drop_time   = rep(x=p$dropped[keep],times=p$n_readouts)
          )
        })

        combined <- do.call(rbind, locked_snapshot_list)
        combined$subject_id <- rep(
          x=seq_len(as.integer(dim(combined)[1] / p$n_readouts)),
          times = p$n_readouts
        )

        if (is.null(combined) || nrow(combined) == 0L) {
          next
        }

        # Add measurement and current time column
        combined$measurement_time <- combined$readout_time + combined$enroll_time
        combined$time <- rep(i, nrow(combined))

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
      invisible(self)
    }
  ) # end public
) # end class

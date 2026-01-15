Trials <- R6::R6Class(
  "Trials",

  public = list(
    # Fields
    name = NULL,
    description = NULL,
    silent = FALSE,
    seed = NULL,
    timer = NULL,        # a Timers object
    population = NULL,   # now a list of Population objects
    locked_data = NULL,
    results = NULL,

    # Constructor
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

    # Fit method: loop through timelist and apply reader conditions
    fit = function() {
      if (is.null(self$timer) || length(self$population) == 0) {
        stop("Timer and population list must be set before running fit()")
      }

      n_timepoints <- self$timer$get_n_timepoints()

      for (i in seq_len(n_timepoints)) {
        tp <- self$timer$get_timepoint(i)

        # Apply enrollment/dropout to each Population object in the list
        for (p in self$population) {
          if (!is.null(tp$enroller)) {
            p$set_enrolled(tp$enroller, time = i)
          }
          if (!is.null(tp$dropper)) {
            p$set_dropped(tp$dropper, time = i)
          }
        }

        # Collect snapshots from all populations
        locked_snapshot <- lapply(self$population, function(p) {
          cbind(p$data, data.frame(
            enroll_t = p$enrolled,
            drop_t = p$dropped
          ))
        })
        n_events <- sum(sapply(self$population, function(p) sum(!is.na(p$enrolled))))

        # Run all reader conditions
        reader_results <- self$timer$run_readers(
          locked_data = locked_snapshot,
          current_time = i,
          n_events = n_events
        )

        if (length(reader_results) > 0) {
          self$locked_data[[paste0("time_", i)]] <- locked_snapshot
          self$results[[paste0("time_", i)]] <- reader_results
        }
      }
    }
  )
)

# Two populations
pop1 <- Population$new("Arm A", data = rnorm(20, mean = 50))
pop2 <- Population$new("Arm B", data = rnorm(20, mean = 55))

# Timers with multiple reader conditions
t <- Timers$new(name = "TrialTimers")
t$add_timepoint(dropper = 2, enroller = 5)
t$add_timepoint(dropper = 1, enroller = 3)

# Add reader conditions
t$add_reader_condition(time = 1, n_events = 5, func = function(ld,time) lapply(ld, mean))
t$add_reader_condition(time = 2, threshold = 10)  # snapshot if >=10 events
# Add in condition flag of TTE

# Trial with list of populations
trial <- Trials$new(name = "Trial A", timer = t, population = list(pop1, pop2))

# Run
trial$fit()

# Inspect results
trial$results

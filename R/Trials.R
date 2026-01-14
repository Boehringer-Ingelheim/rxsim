
Trials <- R6::R6Class(
  "Trials",

  public = list(
    # Fields
    name = NULL,
    description = NULL,
    silent = FALSE,
    seed = NULL,
    timer = NULL,        # a Timers object
    population = NULL,   # list of Population objects
    locked_data = NULL,  # list of snapshots per time
    results = NULL,      # list of reader results per time

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

      for (i in seq_len(n_timepoints)[-1]) {
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


        # Collect raw snapshots from all populations (as a list)
        locked_snapshot_list <- lapply(self$population, function(p) p$get_data())

        # Try to use population names for an 'arm' column; fall back to list names or P1..Pn
        pop_names <- vapply(self$population, function(p) {
          nm <- tryCatch(p$name, error = function(e) NULL)
          if (is.null(nm) || isTRUE(nchar(nm) == 0)) NA_character_ else as.character(nm)
        }, character(1))
        if (all(is.na(pop_names))) {
          pop_names <- names(self$population)
          if (is.null(pop_names) || any(pop_names == "")) {
            pop_names <- paste0("P", seq_along(self$population))
          }
        }

        # Convert each snapshot to data.frame and tag with arm
        dfs <- Map(function(dat, arm_name) {
          df <- as.data.frame(dat, stringsAsFactors = FALSE)
          df$arm <- arm_name
          df
        }, locked_snapshot_list, pop_names)

        # Align columns across populations (union of cols) and row-bind
        all_cols <- Reduce(union, lapply(dfs, names))
        dfs_aligned <- lapply(dfs, function(df) {
          missing_cols <- setdiff(all_cols, names(df))
          if (length(missing_cols)) df[missing_cols] <- NA
          df[, all_cols, drop = FALSE]
        })
        combined <- do.call(rbind, dfs_aligned)
        rownames(combined) <- NULL

        # Add current time column for predicates like time >= 1
        combined$time <- i

        # Run all reader conditions on the combined snapshot
        reader_results <- self$timer$run_readers(
          locked_data  = combined,
          current_time = i
        )

        # Store only if there are results
        if (length(reader_results) > 0) {
          self$locked_data[[paste0("time_", i)]] <- combined
          self$results[[paste0("time_", i)]] <- reader_results
        }
      }
    }
  )
)



# --- Two populations with a common 'value' column ---
set.seed(123)
#long_format
pop1 <- Population$new("Arm A", data = data.frame(value = rnorm(20, mean = 50)))
pop2 <- Population$new("Arm B", data = data.frame(value = rnorm(20, mean = 55)))

# --- Timers with multiple timepoints ---
t <- Timers$new(name = "TrialTimers")  # Use your updated Timers from earlier
t$add_timepoint(dropper = 2, enroller = 3)
t$add_timepoint(dropper = 1, enroller = 3)
t$add_timepoint(dropper = 1, enroller = 3)
t$add_timepoint(dropper = 1, enroller = 3)
t$add_timepoint(dropper = 1, enroller = 3)
t$add_timepoint(dropper = 1, enroller = 3)

# --- Reader conditions (per-reader predicates like dplyr::filter) ---
# IMPORTANT: func must accept (data, current_time). Wrap base functions accordingly.

#event condition
t$add_reader_condition(
  length(value )>4,
  func = function(d, tt) mean(d$value),
  name = "overall_mean"
)

t$add_reader_condition(
  sum(value > 40)>10,
  func = function(d, tt) mean(d$value),
  name = "overall_mean"
)

t$add_reader_condition(
  time == 4,  # uses the 'time' column added by Trials$fit()
  func = function(d, tt) mean(d$value),
  name = "time_mean"
)
t$add_reader_condition(
  sum((value+ time) > 45)>10,
  func = function(d, tt) mean(d$value),
  name = "overall_mean_2"
)

# --- Trial with list of populations ---
trial <- Trials$new(name = "Trial A", timer = t, population = list(pop1, pop2))
self=trial
# --- Run ---
trial$fit()

Trial <- R6::R6Class(
  classname = "Trial",
  public = list(
    # fields
    name = NULL,
    description = NULL,
    silent = FALSE,
    seed = NULL,
    timer = NULL,        # a Timers object
    population = NULL,   # list of Population objects
    locked_data = NULL,  # list of snapshots per time
    results = NULL,      # list of results per time

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

      for (i in seq_len(n_timepoints)[-1]) {
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

        # Collect raw snapshots from all populations (as a list)
        locked_snapshot_list <- lapply(self$population, function(p) {
          subset(cbind(p$data, data.frame(
            enroll_time = p$enrolled,
            drop_time = p$dropped
          )), !is.na(p$enrolled))
        })

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

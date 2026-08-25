# Single source of truth for the columns Trial adds to a snapshot beyond the raw
# population data. NULL = assembled during snapshot construction (enroll_time /
# drop_time come from Population state via an alignment-sensitive cbind and are
# NOT recomputed here). function(df, ctx) = derived by .augment_snapshot from the
# assembled snapshot. List order matters: measurement_time reads enroll_time.
# Adding a derived column here propagates to get_col_names() and both run paths.
.ADDED_COLS <- list(
  enroll_time      = NULL,
  drop_time        = NULL,
  subject_id       = function(df, ctx) rep(seq_len(ctx$n_subj), each = ctx$nr),
  measurement_time = function(df, ctx) df$readout_time + df$enroll_time,
  time             = function(df, ctx) ctx$time
)

# Add the derived columns named in `cols` (default: all non-NULL entries) to `df`.
.augment_snapshot <- function(df, ctx,
                              cols = names(Filter(Negate(is.null), .ADDED_COLS))) {
  for (nm in cols) df[[nm]] <- .ADDED_COLS[[nm]](df, ctx)
  df
}

#' Trial: Simulate a multi‑arm clinical trial
#'
#' @description
#' The `Trial` class coordinates one or more `Population` objects, a `Timer`,
#' and a list of `Condition` objects to simulate a clinical trial.
#'
#' At each unique time defined in the trial's `Timer`, the `Trial`:
#'
#' - applies enrollment and dropout updates to each `Population`
#' - builds a snapshot of all currently enrolled subjects
#' - evaluates each [`Condition`] in `self$conditions` against the snapshot
#' - stores both the snapshot (`locked_data`) and the analysis outputs (`results`)
#'
#' Use `run()` to execute the simulation. Trigger conditions are built with
#' [`Condition`]`$new()` (or helpers [`condition_calendar_time()`] /
#' [`condition_enrollment_fraction()`]) and stored in `trial$conditions`.
#'
#' @seealso [Population], [Timer], [Condition], [collect_results()],
#'   [replicate_trial()], [clone_trial()].
#'
#' @examples
#' # Create two populations
#' popA <- Population$new("A", data = as_population_data(rnorm(10)))
#' popB <- Population$new("B", data = as_population_data(rnorm(12)))
#'
#' # Create a timer and add timepoints
#' t <- Timer$new("Timer")
#' t$add_schedule(data.frame(time = 1, arm = "A", drop = 0L, enroll = 4L))
#' t$add_schedule(data.frame(time = 1, arm = "B", drop = 0L, enroll = 5L))
#' t$add_schedule(data.frame(time = 2, arm = "A", drop = 1L, enroll = 2L))
#' t$add_schedule(data.frame(time = 2, arm = "B", drop = 2L, enroll = 3L))
#'
#' # Build a condition: fire at time >= 2 and count enrolled rows
#' cond <- Condition$new(
#'   where    = calendar_trigger(2),
#'   analysis = function(df, current_time) nrow(df),
#'   name     = "final"
#' )
#'
#' # Create a trial
#' trial <- Trial$new(
#'   name       = "ExampleTrial",
#'   seed       = 123,
#'   timer      = t,
#'   population = list(popA, popB),
#'   conditions = list(cond)
#' )
#'
#' # Run the simulation
#' trial$run()
#'
#' collect_results(trial)
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

    #' @field timer `Timer` object with timepoints.
    timer = NULL,

    #' @field population `list` of [Population] objects, one per arm.
    population = NULL,

    #' @field conditions `list` of [Condition] objects evaluated at each timepoint.
    conditions = NULL,

    #' @field locked_data `list` Snapshots at each timepoint.
    locked_data = NULL,

    #' @field results `list` Analysis outputs per condition.
    results = NULL,

    #' @field adaptive `logical` When `FALSE` (default), uses the fixed fast path:
    #'   enroll/drop times are precomputed deterministically before iteration,
    #'   and snapshots are cheap prefix slices. When `TRUE`, uses the adaptive
    #'   loop: enrollment and dropout are sampled incrementally at each
    #'   timepoint, supporting designs where the schedule may change mid-trial.
    adaptive = FALSE,

    # --- constructor ---
    #' @description
    #' Create a new `Trial` instance.
    #'
    #' @param name `character` Unique identifier for the trial.
    #' @param seed `numeric` or `NULL` Optional random seed for reproducibility.
    #' @param timer `Timer` object defining timepoints.
    #' @param population `list` of [Population] objects, one per arm.
    #' @param locked_data `list` Generated at each `$run()` call.
    #' @param results `list` Analysis outputs generated at each `$run()` call.
    #' @param conditions `list` of [Condition] objects to evaluate at each timepoint.
    #' @param adaptive `logical` When `FALSE` (default), uses the fixed fast path
    #'   (deterministic precompute, prefix snapshots). When `TRUE`, uses the
    #'   adaptive loop (incremental random sampling at each timepoint).
    #' @return A new `Trial` instance.
    #'
    #' @examples
    #' t <- Timer$new(name="simple_timer")
    #' pop <- Population$new(
    #'   name = "simple_pop",
    #'   data = as_population_data(rnorm(5))
    #' )
    #' pop$set_enrolled(5, 1)
    #' Trial$new(name = "simple_trial", timer=t, population = list(pop))
    initialize = function(
      name,
      seed = NULL,
      timer = NULL,
      population = list(), # default empty list
      locked_data = list(),
      conditions =list(),
      results = list(),
      adaptive = FALSE
    ) {
      stopifnot(is.character(name))
      self$name <- name
      self$seed <- seed
      if (!is.null(seed)) set.seed(seed)
      if (!is.null(timer) && !inherits(timer, "Timer")) stop("`timer` must be a Timer instance.")
      stopifnot(is.list(population))
      stopifnot(is.logical(adaptive), length(adaptive) == 1L)

      if (length(population) > 1) {
        readouts <- sapply(population, function(p) p$n_readouts)
        if (length(unique(readouts)) > 1) {
          stop(sprintf(
            "All populations must have the same n_readouts. Found: %s",
            paste(sprintf("%s=%d", sapply(population, function(p) p$name), readouts), collapse = ", ")
          ))
        }
      }

      if (is.null(timer) || nrow(timer$timelist) == 0L) {
        # If timer has no timepoints, extract from population enrollment times
        if (all(sapply(population, function(x) all(is.na(x$enrolled))))) {
          stop("Neither Timer nor Population has enrollment data.")
        } else {
          if (is.null(timer)) timer <- Timer$new(name = paste0(name, "_timer"))
          timepoints <- data.frame(
            time = unlist(lapply(population, function(x) x$enrolled), recursive = FALSE),
            arm = rep(sapply(population, function(x) x$name), sapply(population, function(x) x$n)),
            enroll = 1L,
            drop = 0L
          )
          timer$add_schedule(timepoints)
          self$timer <- timer
        }
      } else {
        self$timer <- timer
      }

      self$population <- population
      self$locked_data <- locked_data
      self$results <- results
      self$conditions <- conditions
      self$adaptive <- adaptive

    },

    # --- methods ---

    #' @description
    #' Execute a trial simulation.
    #'
    #' Dispatches to the fixed fast path (`adaptive = FALSE`, default) or the
    #' adaptive loop (`adaptive = TRUE`). Both paths update `locked_data` and
    #' `results` fields with the same structure.
    #'
    #' **Fixed path** (`adaptive = FALSE`): enroll/drop times are precomputed
    #' deterministically before iteration. Snapshots are cheap prefix slices of
    #' a single combined data frame. Conditions are evaluated with per-condition
    #' exhausted-skip.
    #'
    #' **Adaptive loop** (`adaptive = TRUE`): enrollment and dropout are sampled
    #' randomly at each timepoint, supporting designs where the schedule may
    #' change mid-trial based on interim results.
    #'
    #' @return Updates `locked_data` and `results` fields; returns `self`
    #'   invisibly.
    #'
    #' @seealso [Timer], [Condition], [collect_results()].
    #'
    #' @examples
    #' # Create two populations
    #' popA <- Population$new("A", data = as_population_data(rnorm(10)))
    #' popB <- Population$new("B", data = as_population_data(rnorm(12)))
    #'
    #' # Create a timer and add timepoints
    #' t <- Timer$new("Timer")
    #' t$add_schedule(data.frame(time = 1, arm = "A", drop = 0L, enroll = 4L))
    #' t$add_schedule(data.frame(time = 1, arm = "B", drop = 0L, enroll = 5L))
    #' t$add_schedule(data.frame(time = 2, arm = "A", drop = 1L, enroll = 2L))
    #' t$add_schedule(data.frame(time = 2, arm = "B", drop = 2L, enroll = 3L))
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
    #' collect_results(trial)
    run = function() {
      if (is.null(self$timer) || length(self$population) == 0) {
        stop("Timer and population list must be set before running run()")
      }
      if (self$adaptive) private$run_adaptive() else private$run_fixed()
      invisible(self)
    }
  ), # end public

  private = list(

    # Adaptive loop: incremental random sampling at each timepoint.
    # Used when adaptive = TRUE.
    run_adaptive = function() {
      plan_df <- self$timer$timelist
      if (nrow(plan_df) == 0L) return(invisible(NULL))

      for (i in sort(unique(plan_df$time))) {
        # Apply enrollment/dropout updates to each population at this timepoint
        for (p in self$population) {
          idx <- which(plan_df$arm == p$name & plan_df$time == i)
          if (length(idx) > 0L) {
            enroll_n <- as.integer(sum(plan_df$enroll[idx], na.rm = TRUE))
            drop_n <- as.integer(sum(plan_df$drop[idx], na.rm = TRUE))

            if (enroll_n > 0L && sum(is.na(p$enrolled)) > 0L) {
              p$set_enrolled(enroll_n, time = i)
            }

            if (drop_n > 0L && sum(is.na(p$dropped) & !is.na(p$enrolled)) > 0L) {
              p$set_dropped(drop_n, time = i)
            }
          }
        }

        # Create snapshots of enrolled subjects from all populations
        locked_snapshot_list <- lapply(self$population, function(p) {
          keep <- !is.na(p$enrolled)
          cbind(
            p$data[rep(keep, each = p$n_readouts), , drop = FALSE],
            enroll_time = rep(x=p$enrolled[keep],each=p$n_readouts),
            drop_time   = rep(x=p$dropped[keep],each=p$n_readouts)
          )
        })

        combined <- do.call(rbind, locked_snapshot_list)
        if (is.null(combined) || nrow(combined) == 0L) {
          next
        }

        nr     <- self$population[[1]]$n_readouts
        n_subj <- as.integer(nrow(combined) / nr)
        combined <- .augment_snapshot(combined, list(n_subj = n_subj, nr = nr, time = i))

        # Check all conditions on the combined snapshot
        results <- list()
        for (conds in self$conditions){
          results <- c(results, conds$check_conditions(
            locked_data  = combined,
            current_time = i
          ))
        }
        # Store only if there are results
        if (length(results) > 0) {
          self$locked_data[[paste0("time_", i)]] <- combined
          self$results[[paste0("time_", i)]] <- results
        }
      }
      invisible(NULL)
    },

    # Fixed fast path: deterministic precompute + prefix snapshots.
    # Used when adaptive = FALSE (default).
    run_fixed = function() {
      plan_df <- self$timer$timelist
      if (nrow(plan_df) == 0L) return(invisible(NULL))

      # Precompute enroll_time/drop_time for all populations deterministically.
      for (p in self$population) {
        private$precompute_population(p, plan_df)
      }

      # Build the full combined snapshot (all enrolled subjects).
      full_snap <- private$build_full_snapshot()
      if (is.null(full_snap) || nrow(full_snap) == 0L) return(invisible(NULL))

      unique_times <- sort(unique(plan_df$time))

      # Evaluate conditions at each timepoint via a prefix slice.
      # Per-condition exhausted-skip: skip once trigger_count >= max_triggers.
      for (i in unique_times) {
        # Prefix slice: enrolled subjects visible at time i
        prefix <- full_snap[full_snap$enroll_time <= i, , drop = FALSE]
        if (nrow(prefix) == 0L) next

        # Mask drop_time > i (right-censoring: future drops not yet observed)
        prefix$drop_time[!is.na(prefix$drop_time) & prefix$drop_time > i] <- NA_real_
        prefix <- .augment_snapshot(prefix, list(time = i), cols = "time")

        results <- list()
        for (cond in self$conditions) {
          # Skip exhausted conditions
          if (is.finite(cond$max_triggers) && cond$trigger_count >= cond$max_triggers) next
          results <- c(results, cond$check_conditions(
            locked_data  = prefix,
            current_time = i
          ))
        }

        if (length(results) > 0L) {
          self$locked_data[[paste0("time_", i)]] <- prefix
          self$results[[paste0("time_", i)]] <- results
        }

        # Break once all conditions are exhausted
        all_done <- all(vapply(self$conditions, function(cond)
          is.finite(cond$max_triggers) && cond$trigger_count >= cond$max_triggers,
          logical(1)
        ))
        if (all_done) break
      }
      invisible(NULL)
    },

    # Deterministically assign enroll_time and drop_time for one population.
    # Enroll: expand per-timepoint counts into a sorted per-subject vector
    #   (subject 1 gets earliest enroll_time, subject n gets latest).
    # Drop: assign to earliest-enrolled eligible subjects at each timepoint,
    #   so drop_time >= enroll_time is guaranteed.
    # NULL / zero drop column → all drop_time stay NA.
    precompute_population = function(p, plan_df) {
      arm_rows <- plan_df[plan_df$arm == p$name, , drop = FALSE]

      # --- enroll ---
      enroll_times <- rep(
        arm_rows$time,
        times = as.integer(arm_rows$enroll)
      )
      enroll_times <- sort(enroll_times)
      n_enroll <- min(length(enroll_times), p$n)

      p$enrolled <- rep(NA_real_, p$n)
      if (n_enroll > 0L) {
        p$enrolled[seq_len(n_enroll)] <- enroll_times[seq_len(n_enroll)]
      }

      # --- drop ---
      p$dropped <- rep(NA_real_, p$n)
      has_drops <- !is.null(arm_rows$drop) &&
        any(!is.na(arm_rows$drop) & arm_rows$drop > 0L)
      if (!has_drops) return(invisible(NULL))

      drop_rows <- arm_rows[!is.na(arm_rows$drop) & arm_rows$drop > 0L, , drop = FALSE]
      drop_rows <- drop_rows[order(drop_rows$time), , drop = FALSE]

      # Track which subjects are still drop-eligible at each timepoint:
      # eligible = enrolled (not NA) and not yet dropped (NA dropped).
      for (k in seq_len(nrow(drop_rows))) {
        t_drop <- drop_rows$time[k]
        n_drop  <- as.integer(drop_rows$drop[k])

        # eligible: enrolled at or before this time and not yet dropped
        eligible <- which(!is.na(p$enrolled) & p$enrolled <= t_drop & is.na(p$dropped))
        if (length(eligible) == 0L) next

        n_use <- min(n_drop, length(eligible))
        if (n_use < n_drop) {
          warning(sprintf(
            "Requested %d dropouts at time %s for arm '%s' but only %d eligible; dropping %d.",
            n_drop, t_drop, p$name, length(eligible), n_use
          ), call. = FALSE)
        }
        # Pick eligible subjects uniformly at random (matches set_dropped()).
        pick <- eligible[sample.int(length(eligible), n_use)]
        p$dropped[pick] <- t_drop
      }
      invisible(NULL)
    },

    # Build the full combined snapshot of all enrolled subjects across populations.
    # subject_id is globally unique across arms (same as adaptive path).
    build_full_snapshot = function() {
      nr <- self$population[[1]]$n_readouts
      parts <- lapply(self$population, function(p) {
        keep <- !is.na(p$enrolled)
        if (!any(keep)) return(NULL)
        cbind(
          p$data[rep(keep, each = nr), , drop = FALSE],
          enroll_time = rep(p$enrolled[keep], each = nr),
          drop_time   = rep(p$dropped[keep],  each = nr)
        )
      })
      combined <- do.call(rbind, Filter(Negate(is.null), parts))
      if (is.null(combined) || nrow(combined) == 0L) return(NULL)

      n_subj <- as.integer(nrow(combined) / nr)
      .augment_snapshot(combined, list(n_subj = n_subj, nr = nr),
                        cols = c("subject_id", "measurement_time"))
    }

  ) # end private
) # end class

#' Trial: Simulate a multi‑arm clinical trial
#'
#' @description
#' The `Trial` class coordinates one or more `Population` objects and a `Timer`
#' to simulate a clinical trial.
#'
#' ## Core simulation loop
#' At each timepoint, the `Trial`:
#' - applies enrollment and dropout updates to each `Population`
#' - builds a snapshot of all currently enrolled subjects ("locked" snapshot)
#' - evaluates all conditions in the `Timer`
#' - stores both the snapshot (`locked_data`) and analysis outputs (`results`)
#'
#' ## Adaptive scheduling
#' Unlike a static run-loop that freezes the plan at start, this implementation
#' re-reads the timer's timepoint plan (`timer$timelist`) at each iteration.
#' Therefore, if an analysis function modifies `timer$timelist`, those changes
#' can affect future timepoints during the same `run()` call.
#'
#' ## Audit trail
#' If `track_plan = TRUE`, the `Trial` stores:
#' - `plan_original`: the initial plan at the beginning of the run
#' - `conditions_original`: the condition list at the beginning of the run
#' - `plan_history`: entries recording when the plan changed, including signatures
#'   and (optionally) full before/after plan snapshots.
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
#' # Add a trigger for final analysis at time 2
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
#' @importFrom rlang .data
#' @importFrom dplyr bind_rows filter group_by summarise inner_join arrange mutate
#' @importFrom purrr pwalk map map_dfr
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

    #' @field locked_data `list` Snapshots at each timepoint (stored when conditions return output).
    locked_data = NULL,

    #' @field results `list` Analysis outputs per condition/timepoint.
    results = NULL,

    #' @field plan_original `data.frame` Plan snapshot captured at start of the first tracked run.
    plan_original = NULL,

    #' @field conditions_original `list` Conditions captured at start of the first tracked run.
    conditions_original = NULL,

    #' @field plan_history `list` Plan change log entries (signatures and optional snapshots).
    plan_history = NULL,

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
    population = list(),
    locked_data = list(),
    results = list()
    ) {

      stopifnot(is.character(name), length(name) == 1L)
      self$name <- name
      self$seed <- seed
      if (!is.null(seed)) set.seed(seed)

      if (!is.null(timer) && !inherits(timer, "Timer")) {
        stop("`timer` must be a Timer instance.")
      }
      stopifnot(is.list(population))

      # If timer missing or empty, try to infer timepoints from population enrollment times
      if (is.null(timer) || length(timer$timelist) == 0) {
        if (length(population) == 0) {
          stop("Timer is missing/empty and `population` is empty. Nothing to run.")
        }
        if (all(vapply(population, function(x) all(is.na(x$enrolled)), logical(1)))) {
          stop("Neither Timer nor Population has enrollment data.")
        }

        if (is.null(timer)) timer <- Timer$new(name = paste0(name, "_timer"))

        timepoints <- do.call(rbind, lapply(population, function(p) {
          data.frame(
            time = p$enrolled,
            arm = rep(p$name, length(p$enrolled)),
            enroller = 1L,
            dropper = 0L
          )
        }))
        # Filtering out NA enrollment times
        timepoints <- timepoints[!is.na(timepoints$time), , drop = FALSE]

        add_timepoints(timer, timepoints)
      }

      self$timer <- timer
      self$population <- population
      self$locked_data <- locked_data
      self$results <- results
      self$plan_history <- list()
    },

    # --- methods ---

    #' @description
    #' Execute a trial simulation.
    #' At each unique time defined by the trial's `Timer`:
    #' - Apply enrollment and dropout actions to each `Population`
    #' - Build a combined snapshot of all currently enrolled subjects
    #' - Attach a `time` column to the snapshot
    #' - Evaluate all condition readers via `Timer$check_conditions()`
    #' - Store snapshots and condition outputs under time‑indexed list keys
    #'
    #' ## Adaptive schedule following
    #' This method re-reads the plan (`timer$timelist`) at each iteration so that
    #' any plan updates made within analysis functions take effect immediately
    #' for subsequent timepoints.
    #'
    #' ## Plan auditing
    #' If `track_plan = TRUE`, the method:
    #' - stores `plan_original` and `conditions_original` on first tracked run
    #' - records plan changes into `plan_history`
    #'
    #' @param track_plan `logical` If `TRUE`, capture original plan and log plan changes.
    #' @param verbose `logical` Print progress messages.
    #' @param future_only `logical` If `TRUE`, newly added timepoints <= last processed time are ignored.
    #' @param keep_plan_snapshots `logical` If `TRUE`, store full before/after plan data.frames in `plan_history`.
    #'
    #' @return Invisibly returns `self`, with updated `locked_data` and `results`.
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
    run = function(
    track_plan = TRUE,
    verbose = FALSE,
    future_only = TRUE,
    keep_plan_snapshots = FALSE
    ) {

      if (is.null(self$timer) || length(self$population) == 0) {
        stop("Timer and population list must be set before running run()")
      }

      get_plan <- function() dplyr::bind_rows(self$timer$timelist)

      plan_signature <- function(df) {
        if (is.null(df) || nrow(df) == 0L) return("")
        df2 <- df[, c("time", "arm", "enroller", "dropper")]
        df2 <- df2[order(df2$arm, df2$time), , drop = FALSE]
        paste(apply(df2, 1, paste, collapse = "|"), collapse = ";;")
      }

      # Build a named lookup of populations by arm (and guard duplicates)
      arm_names <- vapply(self$population, function(x) x$name, character(1))
      if (anyDuplicated(arm_names)) {
        stop("Duplicate arm names in population list: ",
             paste(unique(arm_names[duplicated(arm_names)]), collapse = ", "))
      }
      pop_by_arm <- setNames(self$population, arm_names)

      if (track_plan && is.null(self$plan_original)) {
        self$plan_original <- get_plan()
        self$conditions_original <- self$timer$conditions
        self$plan_history <- list()
      }

      processed <- numeric(0)
      last_time <- -Inf

      last_plan <- get_plan()
      if (is.null(last_plan) || nrow(last_plan) == 0L) return(invisible(self))
      last_sig <- plan_signature(last_plan)

      repeat {

        plan_df <- get_plan()
        if (is.null(plan_df) || nrow(plan_df) == 0L) break

        if (track_plan) {
          sig <- plan_signature(plan_df)
          if (!identical(sig, last_sig)) {
            entry <- list(
              changed_after_time = if (is.finite(last_time)) last_time else NA,
              old_signature = last_sig,
              new_signature = sig
            )
            if (keep_plan_snapshots) {
              entry$old_plan <- last_plan
              entry$new_plan <- plan_df
            }
            self$plan_history[[length(self$plan_history) + 1]] <- entry
            if (verbose) message("Plan changed (logged).")
            last_sig <- sig
            last_plan <- plan_df
          }
        }

        remaining_times <- sort(setdiff(unique(plan_df$time), processed))
        if (future_only) remaining_times <- remaining_times[remaining_times > last_time]
        if (length(remaining_times) == 0) break

        i <- remaining_times[1]
        processed <- c(processed, i)
        last_time <- i

        if (verbose) cat("\n=== Trial loop time:", i, "===\n")

        actions_i <- plan_df |>
          dplyr::filter(.data$time == i) |>
          dplyr::group_by(.data$arm) |>
          dplyr::summarise(
            enroller_n = as.integer(sum(.data$enroller, na.rm = TRUE)),
            dropper_n  = as.integer(sum(.data$dropper,  na.rm = TRUE)),
            .groups = "drop"
          )

        # Apply actions (purrr)
        if (nrow(actions_i) > 0L) {
          purrr::pwalk(
            actions_i,
            function(arm, enroller_n, dropper_n) {
              p <- pop_by_arm[[arm]]
              if (is.null(p)) return(NULL)

              if (enroller_n > 0L && sum(is.na(p$enrolled)) > 0L) {
                p$set_enrolled(as.integer(enroller_n), time = i)
              }
              if (dropper_n > 0L && sum(is.na(p$dropped) & !is.na(p$enrolled)) > 0L) {
                p$set_dropped(as.integer(dropper_n), time = i)
              }
              NULL
            }
          )
        }

        # Build snapshot using explicit id_map
        locked_snapshot <- purrr::map_dfr(self$population, function(p) {

          ids <- p$id_map
          if (is.null(ids)) ids <- unique(p$data$id)

          status <- data.frame(
            id = ids,
            enroll_time = p$enrolled,
            drop_time   = p$dropped
          )

          status <- status[!is.na(status$enroll_time), , drop = FALSE]
          if (nrow(status) == 0L) return(NULL)

          dplyr::inner_join(p$data, status, by = "id") |>
            dplyr::mutate(subject_id = .data$id)
        })

        if (is.null(locked_snapshot) || nrow(locked_snapshot) == 0L) next

        locked_snapshot <- locked_snapshot |>
          dplyr::mutate(
            measurement_time = .data$readout_time + .data$enroll_time,
            time = i
          )

        results <- self$timer$check_conditions(
          locked_data  = locked_snapshot,
          current_time = i
        )

        if (length(results) > 0) {
          self$locked_data[[paste0("time_", i)]] <- locked_snapshot
          self$results[[paste0("time_", i)]] <- results
        }
      }

      invisible(self)
    }
  ) # end public
) # end class

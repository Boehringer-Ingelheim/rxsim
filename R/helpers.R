
#' Add Timepoints to a Timer
#'
#' Adds multiple enrollment and dropout events from a data frame.
#'
#' @param timer [`Timer`] instance.
#' @param df `data.frame` with columns: `time` (numeric), `arm` (character),
#'   `enroll` (integer), `drop` (integer).
#'
#' @seealso [Timer], [stochastic_schedule()], [deterministic_schedule()].
#'
#' @export
#'
#' @examples
#' t <- Timer$new(name = "Timer")
#'
#' timepoints <- data.frame(
#'   time = c(1, 2, 3.1, 4, 5, 6),
#'   arm = rep("Arm A", 6),
#'   drop = c(2L, rep(1L, 5)),
#'   enroll = rep(3L, 6)
#' )
#'
#' add_timepoints(t, timepoints)
add_timepoints <- function(timer, df) {
  if (!inherits(timer, "Timer")) stop("`timer` must be a Timer instance.")
  if (!is.data.frame(df)) stop("`df` must be a data.frame with columns: time, arm, enroll, drop")
  invisible(
    sapply(
      split(df, seq_len(nrow(df))),
      function(x) do.call(timer$add_timepoint, x)
    )
  )
  invisible(timer)
}


#' Collect Trial Results Across Replicates
#'
#' Gathers analysis outputs from one or more `Trial` objects into a single
#' tidy data frame. Works with any number of named analyses (e.g., both
#' an interim and a final) and any number of replicates.
#'
#' @param trials A `Trial` R6 object **or** a `list` of `Trial` objects
#'   (as returned by [replicate_trial()]).
#' @param analysis `character` or `NULL`. When supplied, only analyses whose
#'   name matches one of these values are included. Defaults to `NULL`
#'   (all analyses).
#'
#' @return A `data.frame` with columns:
#'   - `replicate` `integer` Index of the trial replicate (1-based).
#'   - `timepoint` `numeric` Calendar time at which the analysis fired.
#'   - `analysis`  `character` Name of the analysis (as given in
#'     `analysis_generators` or via a trigger helper).
#'   - Additional columns from the value returned by each analysis function.
#'
#' @details
#' Each analysis function may return either a `data.frame` (the standard
#' pattern) or a named `list`; both are coerced to a single-row data frame
#' and stacked. Analyses that return `NULL` or `NA` are silently skipped.
#'
#' When `trials` is a single `Trial` object (e.g., from a one-off
#' `Trial$new()` + `$run()` call), the `replicate` column is always `1`.
#'
#' @seealso [replicate_trial()], [run_trials()], [Trial].
#'
#' @export
#'
#' @examples
#' # --- replicated trial ---
#' pop_gens <- list(
#'   control   = function(n) as_population_data(rnorm(n)),
#'   treatment = function(n) as_population_data(rnorm(n, 0.5))
#' )
#' an_gens <- list(
#'   final = list(
#'     trigger  = count_trigger("enroll_time", ">=", 20L),
#'     analysis = function(df, current_time) {
#'       data.frame(mean_ctrl = mean(df$data[df$arm == "control"]))
#'     }
#'   )
#' )
#' trials <- replicate_trial(
#'   trial_name = "ex", sample_size = 20L,
#'   arms = c("control", "treatment"), allocation = c(1, 1),
#'   enrollment = function(n) rexp(n, 1), dropout = function(n) rexp(n, 0.01),
#'   analysis_generators = an_gens, population_generators = pop_gens, n = 3
#' )
#' run_trials(trials)
#' collect_results(trials)
#'
#' # --- filter to a specific analysis name ---
#' collect_results(trials, analysis = "final")
collect_results <- function(trials, analysis = NULL) {
  if (inherits(trials, "Trial")) trials <- list(trials)
  if (!is.list(trials)) stop("`trials` must be a Trial object or a list of Trial objects.")
  not_trials <- !vapply(trials, function(x) inherits(x, "Trial"), logical(1))
  if (any(not_trials)) {
    stop(sprintf(
      "`trials` must be a Trial object or a list of Trial objects. Element(s) %s are not Trial objects.",
      paste(which(not_trials), collapse = ", ")
    ))
  }

  rows <- lapply(seq_along(trials), function(i) {
    results <- trials[[i]]$results
    if (length(results) == 0L) return(NULL)

    tp_rows <- lapply(names(results), function(tp_name) {
      analyses <- results[[tp_name]]

      if (!is.null(analysis)) {
        analyses <- analyses[names(analyses) %in% analysis]
      }
      if (length(analyses) == 0L) return(NULL)

      an_rows <- lapply(names(analyses), function(an_name) {
        val <- analyses[[an_name]]
        if (is.null(val) || (length(val) == 1L && is.na(val))) return(NULL)

        df <- if (is.data.frame(val)) {
          val
        } else {
          as.data.frame(as.list(val), stringsAsFactors = FALSE)
        }

        cbind(
          data.frame(
            replicate = i,
            timepoint = as.numeric(sub("time_", "", tp_name)),
            analysis  = an_name,
            stringsAsFactors = FALSE,
            row.names = NULL
          ),
          df
        )
      })

      dplyr::bind_rows(an_rows)
    })

    dplyr::bind_rows(tp_rows)
  })

  result <- dplyr::bind_rows(rows)
  if (!is.null(result)) rownames(result) <- NULL
  result
}

#' Create a Population-Compatible Data Frame from a Vector
#'
#' Converts a numeric vector to a data frame with the columns required by
#' [`Population`]: `id`, `data`, and `readout_time`.
#'
#' @param data `numeric` vector of endpoint values (one per subject).
#'
#' @return `data.frame` with columns: `id` (integer), `data` (numeric),
#'   `readout_time` (0 for all subjects).
#'
#' @seealso [Population].
#'
#' @export
#'
#' @examples
#' as_population_data(rnorm(5))
as_population_data <- function(data) data.frame(
      id = seq_along(data),
      data = data,
      readout_time = 0
)

#' Extract Column Names from Populations
#'
#' Collects all data frame column names from one or more populations.
#'
#' @param populations [`Population`] object or `list` of [`Population`] objects.
#'
#' @return `character` vector of unique column names.
#'
#' @seealso [Population].
#'
#' @export
#'
#' @examples
#' pop1 <- Population$new(name = "P1", data = data.frame(
#' id = 1:10,
#' age = runif(10, 20, 60),
#' readout_time = 0
#' ))
#' pop2 <- Population$new(name = "P2", data = data.frame(
#' id = 1:10,
#' weight = runif(10, 150, 250),
#' readout_time = 0
#' ))
#' get_col_names(list(pop1, pop2))
get_col_names <- function(populations) {
  col_names <- NULL
  if (is.list(populations)) {
    for (p in populations) {
      col_names <- c(col_names, colnames(p$data))
    }
  } else {
    col_names <- c(col_names, colnames(populations$data))
  }

  col_names <- c(col_names, "time", "enroll_time", "drop_time", "measure_time")
  return(unique(col_names))
}

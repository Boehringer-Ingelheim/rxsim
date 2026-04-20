#' Trigger Analysis at a Calendar Time
#'
#' Builds a [`Condition`] that fires when the trial clock reaches a specified
#' calendar time. The returned `Condition` should be passed to
#' `Trial$new(conditions = list(...))`.
#'
#' @param cal_time `numeric` Calendar time(s) at which to trigger.
#' @param analysis `function` or `NULL` Optional analysis function called as
#'   `analysis(filtered_data, current_time)`. If `NULL`, the filtered snapshot
#'   is returned as-is with a warning.
#' @param name `character` or `NULL` Result key. Defaults to
#'   `"cal_time_<cal_time>"`.
#'
#' @return A [`Condition`] object.
#'
#' @seealso [Condition], [trigger_by_fraction()], [Trial].
#'
#' @export
#'
#' @examples
#' cond <- trigger_by_calendar(
#'   cal_time = 12,
#'   analysis = function(df, current_time) {
#'     data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
#'   }
#' )
#'
#' @importFrom rlang quos
#' @importFrom dplyr .data
trigger_by_calendar <- function(cal_time, analysis = NULL, name = NULL) {
  if (missing(cal_time)) stop("`cal_time` is required.")
  stopifnot(is.numeric(cal_time))
  if (is.null(name)) name <- paste0("cal_time_", paste(cal_time, collapse = "_"))

  Condition$new(
    where    = rlang::quos(.data$time %in% !!cal_time),
    analysis = analysis,
    name     = name
  )
}

#' Trigger Analysis at a Sample Fraction
#'
#' Builds a [`Condition`] that fires when a given fraction of the target sample
#' has been enrolled. The returned `Condition` should be passed to
#' `Trial$new(conditions = list(...))`.
#'
#' @param fraction `numeric` Sample fraction (0 < fraction <= 1).
#' @param sample_size `integer` Target sample size.
#' @param analysis `function` or `NULL` Optional analysis function called as
#'   `analysis(filtered_data, current_time)`. If `NULL`, the filtered snapshot
#'   is returned as-is with a warning.
#' @param name `character` or `NULL` Result key. Defaults to
#'   `"frac_<fraction>"`.
#'
#' @return A [`Condition`] object.
#'
#' @seealso [Condition], [trigger_by_calendar()], [Trial].
#'
#' @export
#'
#' @examples
#' cond <- trigger_by_fraction(
#'   fraction    = 0.5,
#'   sample_size = 100,
#'   analysis    = function(df, current_time) {
#'     data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
#'   }
#' )
#'
#' @importFrom rlang quos
#' @importFrom dplyr .data
trigger_by_fraction <- function(fraction, sample_size, analysis = NULL, name = NULL) {
  if (missing(fraction) || missing(sample_size)) stop("`fraction` and `sample_size` are required.")
  stopifnot(is.numeric(sample_size) && length(sample_size) == 1L)
  stopifnot(fraction > 0 && fraction <= 1)
  if (is.null(name)) name <- paste0("frac_", fraction)
  target_n <- fraction * sample_size

  Condition$new(
    where    = rlang::quos(sum(!is.na(.data$enroll_time)) >= !!target_n),
    analysis = analysis,
    name     = name
  )
}

#' Add Timepoints to a Timer
#'
#' Adds multiple enrollment and dropout events from a data frame.
#'
#' @param timer [`Timer`] instance.
#' @param df `data.frame` with columns: `time` (numeric), `arm` (character),
#'   `enroller` (integer), `dropper` (integer).
#'
#' @seealso [Timer], [gen_plan()], [gen_timepoints()].
#'
#' @export
#'
#' @examples
#' t <- Timer$new(name = "Timer")
#'
#' timepoints <- data.frame(
#'   time = c(1, 2, 3.1, 4, 5, 6),
#'   arm = rep("Arm A", 6),
#'   dropper = c(2L, rep(1L, 5)),
#'   enroller = rep(3L, 6)
#' )
#'
#' add_timepoints(t, timepoints)
add_timepoints <- function(timer, df) {
  if (!inherits(timer, "Timer")) stop("`timer` must be a Timer instance.")
  if (!is.data.frame(df)) stop("`df` must be a data.frame with columns: time, arm, enroller, dropper")
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
#'   control   = function(n) vector_to_dataframe(rnorm(n)),
#'   treatment = function(n) vector_to_dataframe(rnorm(n, 0.5))
#' )
#' an_gens <- list(
#'   final = list(
#'     trigger  = rlang::exprs(sum(!is.na(enroll_time)) >= 20L),
#'     analysis = function(df, timer) {
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

      do.call(rbind, an_rows)
    })

    do.call(rbind, tp_rows)
  })

  result <- do.call(rbind, rows)
  if (!is.null(result)) rownames(result) <- NULL
  result
}

#' Format Trial Results as a Data Frame
#'
#' Converts trial results to a single data frame with all measurements.
#'
#' @param results `list` Trial results (nested by time).
#'
#' @return `data.frame` with columns: `time` and measurement columns.
#'
#' @seealso [Trial] for generating results.
#'
#' @export
prettify_results <- function(results) {
  stopifnot(is.list(results))
  all_cols <- unique(unlist(lapply(results, names)))

  df <- do.call(rbind, lapply(names(results), function(nm) {
    row <- results[[nm]]
    row[setdiff(all_cols, names(row))] <- NA
    out <- as.data.frame(as.list(row), stringsAsFactors = FALSE)
    out$time <- as.numeric(sub("time_", "", nm))
    out
  }))

  df[c("time", setdiff(names(df), "time"))]
}

#' Convert Vector to Population Data Frame
#'
#' Formats a numeric vector as a population data frame.
#'
#' @param data `numeric` vector of population values.
#'
#' @return `data.frame` with columns: `id`, `data`, `readout_time`.
#'
#' @seealso [Population].
#'
#' @export
vector_to_dataframe <- function(data) data.frame(
      id = seq_along(data),
      data = data,
      readout_time = 0
)

#' Generate Piecewise-Linear Enrollment and Dropout Plan
#'
#' Creates a time-indexed plan with piecewise constant enrollment and dropout rates.
#'
#' @param sample_size `integer` Trial sample size.
#' @param arms `character` vector of arm identifiers.
#' @param allocation `numeric` vector of allocation ratios.
#' @param enrollment `list` with `end_time` (numeric endpoints) and
#'   `rate` (numeric subjects/unit time for each period).
#' @param dropout `list` with `end_time` and `rate` (same structure).
#'
#' @return `data.frame` with columns: `time`, `arm`, `enroller`, `dropper`.
#'
#' @seealso [gen_plan()] for random inter-event times, [add_timepoints()].
#'
#' @export
#'
#' @examples
#' gen_timepoints(
#'   sample_size = 100,
#'   arms = c("A", "B"),
#'   allocation = c(2, 1),
#'   enrollment = list(
#'     end_time = c(4, 8, 12),
#'     rate = c(6, 12, 18)
#'   ),
#'   dropout = list(
#'     end_time = c(5, 9, 13),
#'     rate = c(0, 3, 6)
#'   )
#' )
#'
#' @importFrom utils tail
#' @importFrom rlang :=
#' @importFrom dplyr .data
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr arrange
gen_timepoints <- function(sample_size, arms, allocation, enrollment, dropout) {
  # Input validation
  if (!is.numeric(sample_size) || length(sample_size) != 1L || sample_size <= 0) {
    stop("`sample_size` must be a single positive number.")
  }
  if (!is.character(arms) || length(arms) == 0L) {
    stop("`arms` must be a non-empty character vector.")
  }
  if (!is.numeric(allocation) || length(allocation) != length(arms)) {
    stop("`allocation` must be a numeric vector with same length as `arms`.")
  }
  if (!is.list(enrollment) || !all(c("end_time", "rate") %in% names(enrollment))) {
    stop("`enrollment` must be a list with 'end_time' and 'rate'.")
  }
  if (!is.list(dropout) || !all(c("end_time", "rate") %in% names(dropout))) {
    stop("`dropout` must be a list with 'end_time' and 'rate'.")
  }

  # Calculate arm allocation ratios
  n_arms <- length(arms)
  ratio <- allocation / sum(allocation)
  names(ratio) <- arms

  # Define timeline endpoint
  end <- max(utils::tail(enrollment$end_time, 1), utils::tail(dropout$end_time, 1))

  # Calculate target enrollment per arm
  target <- as.integer(round(ratio * sample_size))
  names(target) <- arms

  # Adjust for rounding: add missing subjects
  if (sample_size - sum(target) > 0) {
    addition <- table(sample(
      seq_len(n_arms),
      sample_size - sum(target),
      replace = TRUE,
      prob = ratio
    )) |> as.vector()
    target <- target + addition
  }

  # Adjust for rounding: remove excess subjects
  if (sample_size - sum(target) < 0) {
    remove <- table(sample(
      seq_len(n_arms),
      sum(target) - sample_size,
      replace = TRUE,
      prob = ratio
    )) |> as.vector()
    target <- target - remove
  }

  # Pad rate vectors to match timeline endpoint
  pad <- function(x, end) {
    if (tail(x$end_time, 1) != end) {
      x$end_time <- c(x$end_time, end)
      x$rate <- c(x$rate, 0)
      x
    } else {
      x
    }
  }

  enrollment <- pad(enrollment, end)
  dropout <- pad(dropout, end)

  # Calculate duration of each time period
  get_durations <- function(x) (c(0, x) - dplyr::lag(c(0, x)))[-1]

  # Create base schedule (may exceed target enrollment)
  df <- data.frame(
    time = rep(seq_len(end), n_arms),
    arm = rep(arms, each = end),
    enroller = rep(
      as.vector(outer(enrollment$rate, ratio)),
      rep(get_durations(enrollment$end_time), n_arms)
    ) |> as.integer(),
    dropper = rep(
      as.vector(outer(dropout$rate, ratio)),
      rep(get_durations(dropout$end_time), n_arms)
    ) |> as.integer()
  )

  # Identify undershooting periods (cumulative < target)
  checks <- df |>
    dplyr::group_by(.data$arm) |>
    dplyr::mutate(cum = cumsum(.data$enroller)) |>
    dplyr::ungroup() |>
    dplyr::mutate(under = .data$cum < target[.data$arm])

  # Find final undershooting timepoint per arm
  next_t <- checks |>
    dplyr::group_by(.data$arm) |>
    dplyr::filter(.data$under) |>
    dplyr::filter(dplyr::row_number() == dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::select(.data$time)
  next_t <- as.vector(next_t$time) + rep(1, n_arms)
  names(next_t) <- arms

  # Calculate enrollment gap per arm
  next_enroll <- checks |>
    dplyr::group_by(.data$arm) |>
    dplyr::filter(.data$under) |>
    dplyr::filter(dplyr::row_number() == dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::select(.data$cum)
  next_enroll <- target - as.vector(next_enroll$cum)
  names(next_enroll) <- arms

  # Create correction row(s) to reach target enrollment
  df_add <- data.frame(
    time = next_t,
    arm = arms,
    enroller = next_enroll,
    dropper = as.integer(round(dropout$rate[findInterval(next_t, dropout$end_time)] * ratio))
  )

  # Combine schedule and corrections, sort by arm and time
  checks |>
    dplyr::filter(.data$under) |>
    dplyr::bind_rows(df_add) |>
    dplyr::select(-c(.data$cum, .data$under)) |>
    dplyr::group_by(.data$arm) |>
    dplyr::arrange(.data$time, .by_group = TRUE) |>
    dplyr::ungroup()
}

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

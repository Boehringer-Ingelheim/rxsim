#' Trigger Analysis at a Calendar Time
#'
#' Adds an analysis trigger at a specified calendar time.
#'
#' @param cal_time `numeric` Calendar time(s) to trigger.
#' @param timer [`Timer`] instance to update.
#' @param analysis `function` or `NULL` Optional function to apply.
#'
#' @returns Invisible [`Timer`].
#'
#' @seealso [Timer], [trigger_by_fraction()].
#'
#' @export
#'
#' @examples
#' t <- Timer$new("timer")
#' trigger_by_calendar(2, t, analysis = function(df, current_time) nrow(df))
#' 
#' @importFrom rlang :=
#' @importFrom dplyr .data
trigger_by_calendar <- function(cal_time, timer, analysis = NULL) {
  timer$add_condition(
    .data$time %in% cal_time,
    analysis = analysis,
    name = paste0("cal_time_", do.call(paste, c(cal_time, sep = "_") |> as.list()))
  )
}

#' Trigger Analysis at a Sample Fraction
#'
#' Adds an analysis trigger when a fraction of the target sample enrolls.
#'
#' @param fraction `numeric` Sample fraction (0 < fraction <= 1).
#' @param timer [`Timer`] instance to update.
#' @param sample_size `integer` Target sample size.
#' @param analysis `function` or `NULL` Optional function to apply.
#'
#' @returns Invisible [`Timer`].
#'
#' @seealso [Timer], [trigger_by_calendar()].
#'
#' @export
#'
#' @examples
#' t <- Timer$new("timer")
#' trigger_by_fraction(0.5, t, sample_size = 100, analysis = function(df, current_time) nrow(df))
#' 
#' @importFrom rlang :=
#' @importFrom dplyr .data
trigger_by_fraction <- function(fraction, timer, sample_size, analysis = NULL) {
  stopifnot(fraction > 0 && fraction <= 1)

  timer$add_condition(
    sum(!is.na(.data$enroll_time)) >= fraction * sample_size,
    analysis = analysis,
    name = paste0("frac_", fraction)
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
  sapply(
    split(df, 1:nrow(df)),
    function(x) do.call(timer$add_timepoint, x)
  )
  invisible(timer)
}


#' Format Trial Results as a Data Frame
#'
#' Converts trial results to a single data frame with all measurements.
#'
#' @param results `list` Trial results (nested by time).
#'
#' @returns `data.frame` with columns: `time` and measurement columns.
#'
#' @seealso [Trial] for generating results.
#'
#' @export
prettify_results <- function(results) {
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
#' @returns `data.frame` with columns: `id`, `data`, `readout_time`.
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
#' @returns `data.frame` with columns: `time`, `arm`, `enroller`, `dropper`.
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
#' @returns `character` vector of unique column names.
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

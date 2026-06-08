#' Generate a Stochastic Enrollment and Dropout Schedule
#'
#' Creates a time-indexed schedule of enrollment and dropout events across
#' arms by sampling inter-event times from user-supplied distribution
#' functions. Each call produces a different realization, capturing
#' natural variability in study timelines.
#'
#' Use this when trial-duration variability is substantively important.
#' For a fixed, reproducible schedule see [deterministic_schedule()].
#'
#' @param sample_size `integer` Trial sample size.
#' @param arms `character` vector of arm identifiers.
#' @param allocation `numeric` vector of allocation ratios.
#' @param enrollment `function` that takes `n` and returns `n` inter-arrival
#'   times (e.g. `function(n) rexp(n, rate = 1)`).
#' @param dropout `function` that takes `n` and returns `n` inter-dropout
#'   times.
#'
#' @return `data.frame` with columns: `time`, `arm`, `enroll` (always 1),
#'   `drop` (always 0 or 1). One row per subject event, sorted by `time`.
#'
#' @seealso [deterministic_schedule()] for piecewise-constant rates,
#'   [add_timepoints()] to attach to a `Timer`.
#'
#' @export
#'
#' @examples
#' stochastic_schedule(
#'   sample_size = 100,
#'   arms = c("A", "B"),
#'   allocation = c(2, 1),
#'   enrollment = function(n) rexp(n, rate = 0.5),
#'   dropout    = function(n) rexp(n, rate = 0.1)
#' )
#'
#' @importFrom rlang :=
#' @importFrom dplyr .data
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr arrange
.validate_schedule_common_args <- function(sample_size, arms, allocation) {
  if (!is.numeric(sample_size) || length(sample_size) != 1L || sample_size <= 0) {
    stop("`sample_size` must be a single positive number.")
  }
  if (!is.character(arms) || length(arms) == 0L) {
    stop("`arms` must be a non-empty character vector.")
  }
  if (!is.numeric(allocation) || length(allocation) != length(arms)) {
    stop("`allocation` must be a numeric vector with same length as `arms`.")
  }
}

.allocate_targets <- function(sample_size, arms, ratio) {
  n_arms <- length(arms)
  target <- as.integer(round(ratio * sample_size))
  names(target) <- arms

  diff <- sample_size - sum(target)
  if (diff > 0L) {
    add_idx <- sample(seq_len(n_arms), diff, replace = TRUE, prob = ratio)
    target <- target + stats::setNames(tabulate(add_idx, nbins = n_arms), arms)
  } else if (diff < 0L) {
    remove_idx <- sample(seq_len(n_arms), -diff, replace = TRUE, prob = ratio)
    target <- target - stats::setNames(tabulate(remove_idx, nbins = n_arms), arms)
  }

  if (any(target < 0L) || sum(target) != sample_size) {
    stop("Enrollment target generation failed: arm targets do not sum to sample_size.")
  }
  target
}

stochastic_schedule <- function(sample_size, arms, allocation, enrollment, dropout) {

  # Input validation
  .validate_schedule_common_args(sample_size, arms, allocation)
  if (!is.function(enrollment) || !is.function(dropout)) {
    stop("`enrollment` and `dropout` must be functions.")
  }

  # Calculate arm allocation ratios
  n_arms <- length(arms)
  ratio <- allocation / sum(allocation)
  names(ratio) <- arms

  # Calculate target enrollment per arm
  target <- .allocate_targets(sample_size, arms, ratio)

  # Generate enrollment and dropout inter-event times
  enroll_events <- enrollment(sample_size)
  drop_events <- dropout(sample_size)

  # Shuffle arms to randomize allocation
  shuffled_arms <- sample(
    rep(arms, times = target),
    sample_size,
    replace = FALSE
  )

  # Create enrollment events (cumulative timing)
  df_enroll <- data.frame(
    time = cumsum(enroll_events),
    arm = shuffled_arms,
    enroll = 1L,
    drop = 0L
  )

  # Create dropout events (cumulative timing)
  df_drop <- data.frame(
    time = cumsum(drop_events),
    arm = sample(arms, sample_size, replace = TRUE, prob = ratio),
    enroll = 0L,
    drop = 1L
  )

  # Combine and sort by time
  rbind(df_enroll, df_drop) |> dplyr::arrange(.data$time)
}

#' Generate a Deterministic Enrollment and Dropout Schedule
#'
#' Creates a time-indexed schedule with piecewise-constant enrollment and
#' dropout rates. Every call with the same inputs returns the same schedule,
#' so all replicates follow an identical enrollment pattern.
#'
#' Use this when you have a well-characterized operational plan and want to
#' isolate endpoint and analysis variability from timing variability.
#' For a stochastic (random) schedule see [stochastic_schedule()].
#'
#' @param sample_size `integer` Trial sample size.
#' @param arms `character` vector of arm identifiers.
#' @param allocation `numeric` vector of allocation ratios.
#' @param enrollment `list` with `end_time` (numeric period endpoints) and
#'   `rate` (subjects/unit time for each period).
#' @param dropout `list` with `end_time` and `rate` (same structure).
#'
#' @return `data.frame` with columns: `time` (integer period), `arm`,
#'   `enroll` (subjects enrolled in that period), `drop` (subjects
#'   dropped). Aggregated counts  -  multiple subjects per row.
#'
#' @seealso [stochastic_schedule()] for random inter-event times,
#'   [add_timepoints()].
#'
#' @export
#'
#' @examples
#' deterministic_schedule(
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
deterministic_schedule <- function(sample_size, arms, allocation, enrollment, dropout) {
  # Input validation
  .validate_schedule_common_args(sample_size, arms, allocation)
  if (!is.list(enrollment) || !all(c("end_time", "rate") %in% names(enrollment))) {
    stop("`enrollment` must be a list with 'end_time' and 'rate'.")
  }
  if (!is.list(dropout) || !all(c("end_time", "rate") %in% names(dropout))) {
    stop("`dropout` must be a list with 'end_time' and 'rate'.")
  }

  # Calculate arm allocation ratios
  ratio <- allocation / sum(allocation)
  names(ratio) <- arms

  # Define timeline endpoint
  end <- max(utils::tail(enrollment$end_time, 1), utils::tail(dropout$end_time, 1))

  # Calculate target enrollment per arm
  target <- .allocate_targets(sample_size, arms, ratio)

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
  get_durations <- function(x) diff(c(0, x))
  n_arms <- length(arms)

  # Create base schedule (may exceed target enrollment)
  df <- data.frame(
    time = rep(seq_len(end), n_arms),
    arm = rep(arms, each = end),
    enroll = rep(
      as.vector(outer(enrollment$rate, ratio)),
      rep(get_durations(enrollment$end_time), n_arms)
    ) |> as.integer(),
    drop = rep(
      as.vector(outer(dropout$rate, ratio)),
      rep(get_durations(dropout$end_time), n_arms)
    ) |> as.integer()
  )

  # Identify undershooting periods (cumulative enrollment < target)
  checks <- df |>
    dplyr::group_by(.data$arm) |>
    dplyr::mutate(total_enrolled = cumsum(.data$enroll)) |>
    dplyr::ungroup() |>
    dplyr::mutate(below_target = .data$total_enrolled < target[.data$arm])

  last_under <- do.call(rbind, lapply(arms, function(a) {
    arm_rows <- checks[checks$arm == a & checks$below_target, , drop = FALSE]
    if (nrow(arm_rows) == 0L) {
      data.frame(arm = a, time = 0, total_enrolled = 0)
    } else {
      tail(arm_rows[, c("arm", "time", "total_enrolled"), drop = FALSE], 1)
    }
  }))
  rownames(last_under) <- NULL

  next_t <- as.integer(last_under$time + 1L)
  names(next_t) <- last_under$arm
  next_enroll <- target[last_under$arm] - as.integer(last_under$total_enrolled)
  names(next_enroll) <- last_under$arm

  # Create correction row(s) to reach target enrollment
  df_add <- data.frame(
    time = as.integer(next_t[arms]),
    arm = arms,
    enroll = as.integer(next_enroll[arms]),
    drop = as.integer(round(dropout$rate[findInterval(next_t, dropout$end_time)] * ratio))
  )

  # Combine schedule and corrections, sort by arm and time
  checks |>
    dplyr::filter(.data$below_target) |>
    dplyr::bind_rows(df_add) |>
    dplyr::filter(.data$enroll > 0L | .data$drop > 0L) |>
    dplyr::select(-dplyr::all_of(c("total_enrolled", "below_target"))) |>
    dplyr::group_by(.data$arm) |>
    dplyr::arrange(.data$time, .by_group = TRUE) |>
    dplyr::ungroup()
}
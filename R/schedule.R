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
stochastic_schedule <- function(sample_size, arms, allocation, enrollment, dropout) {

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
  if (!is.function(enrollment) || !is.function(dropout)) {
    stop("`enrollment` and `dropout` must be functions.")
  }

  # Calculate arm allocation ratios
  n_arms <- length(arms)
  ratio <- allocation / sum(allocation)
  names(ratio) <- arms

  # Calculate target enrollment per arm
  target <- as.integer(round(ratio * sample_size))
  names(target) <- arms

  # Adjust for rounding: remove excess subjects
  if (sample_size - sum(target) < 0) {
    remove_idx <- sample(
      seq_len(n_arms),
      sum(target) - sample_size,
      replace = TRUE,
      prob = ratio
    )
    remove <- tabulate(remove_idx, nbins = n_arms)
    names(remove) <- arms
    target <- target - remove
  }

  # Adjust for rounding: add missing subjects
  if (sample_size - sum(target) > 0) {
    add_idx <- sample(
      seq_len(n_arms),
      sample_size - sum(target),
      replace = TRUE,
      prob = ratio
    )
    addition <- tabulate(add_idx, nbins = n_arms)
    names(addition) <- arms
    target <- target + addition
  }

  if (any(target < 0L) || sum(target) != sample_size) {
    stop("Enrollment target generation failed: arm targets do not sum to sample_size.")
  }

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
    enroll = rep(
      as.vector(outer(enrollment$rate, ratio)),
      rep(get_durations(enrollment$end_time), n_arms)
    ) |> as.integer(),
    drop = rep(
      as.vector(outer(dropout$rate, ratio)),
      rep(get_durations(dropout$end_time), n_arms)
    ) |> as.integer()
  )

  # Identify undershooting periods (cumulative < target)
  checks <- df |>
    dplyr::group_by(.data$arm) |>
    dplyr::mutate(cum = cumsum(.data$enroll)) |>
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
    enroll = next_enroll,
    drop = as.integer(round(dropout$rate[findInterval(next_t, dropout$end_time)] * ratio))
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
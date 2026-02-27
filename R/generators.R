#' Generate Trial Enrollment and Dropout Plan
#'
#' Creates a time-indexed plan of enrollment and dropout events across arms.
#'
#' @param sample_size `integer` Trial sample size.
#' @param arms `character` vector of arm identifiers.
#' @param allocation `numeric` vector of allocation ratios.
#' @param enrollment `function` generating inter-enrollment times (PDF).
#' @param dropout `function` generating inter-dropout times (PDF).
#'
#' @returns `data.frame` with columns: `time`, `arm`, `enroller`, `dropper`.
#'
#' @seealso [gen_timepoints()] for piecewise-constant rates, [add_timepoints()]
#'   to attach generated plans to a `Timer`.
#'
#' @export
#'
#' @examples
#' gen_plan(
#'   sample_size = 100,
#'   arms = c("A", "B"),
#'   allocation = c(2,1),
#'   enrollment = function(n) rexp(n, rate = 0.5),
#'   dropout = function(n) rexp(n, rate = 0.1)
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
gen_plan <- function(sample_size, arms, allocation, enrollment, dropout) {

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
    enroller = 1L,
    dropper = 0L
  )

  # Create dropout events (cumulative timing)
  df_drop <- data.frame(
    time = cumsum(drop_events),
    arm = sample(arms, sample_size, replace = TRUE, prob = ratio),
    enroller = 0L,
    dropper = 1L
  )

  # Combine and sort by time
  rbind(df_enroll, df_drop) |> dplyr::arrange(.data$time)
}

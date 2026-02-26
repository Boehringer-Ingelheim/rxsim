#' Generate plan for trial with allocation, enrollment and dropout
#'
#' @param sample_size `integer` trial sample size
#' @param arms `character` vector of unique identifier of arms
#' @param allocation `numeric` vector of arm allocations
#' @param enrollment `function` PDF of inter-enrollment time
#' @param dropout `function` PDF of inter-dropout time
#'
#' @returns `data.frame` timepoints that may be passed to [add_timepoints()]
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
#' @importFrom utils tail
#' @importFrom rlang :=
#' @importFrom dplyr .data
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr arrange
gen_plan <- function(sample_size, arms, allocation, enrollment, dropout) {

  # determine number of arms
  n_arms <- length(arms)

  # get valid weights
  ratio <- allocation / sum(allocation)
  names(ratio) <- arms

  # get enrollment targets
  target <- as.integer(round(ratio * sample_size))
  names(target) <- arms

  # handle removal
  if (sample_size - sum(target) < 0) {
    remove <- table(sample(
      seq_len(n_arms),
      sum(target) - sample_size,
      replace = TRUE,
      prob = ratio
    )) |> as.vector()
    target <- target - remove
  }

    # handle additions
  if (sample_size - sum(target) > 0) {
    addition <- table(sample(
      seq_len(n_arms),
      sample_size - sum(target),
      replace = TRUE,
      prob = ratio
    )) |> as.vector()
  target <- target + addition
  }

    enroll_events <- enrollment(sample_size)
    drop_events <- dropout(sample_size)

  df_enroll <- data.frame(
    time = cumsum(enroll_events),
    arm = sample(arms, sample_size, replace = TRUE, prob = ratio),
    enroller = 1,
    dropper = 0
  )

    df_drop <- data.frame(
    time = cumsum(drop_events),
    arm = sample(arms, sample_size, replace = TRUE, prob = ratio),
    enroller = 0,
    dropper = 1
  )

rbind(df_enroll, df_drop) |> dplyr::arrange(.data$time)
}
#' Add multiple timepoints from a data frame
#'
#' @param timer an instance of [Timer]
#' @param df data frame with following columns:
#' - time (`numeric`) calendar time
#' - arm (`character`) unique arm key
#' - dropper (`integer`) # of subjects to drop
#' - enroller (`integer`) # of subjects to enroll
#'
#' @export
#'
#' @examples
#' t <- Timer$new(name = "Timer")
#'
#' timepoints <- data.frame(
#' time = c(1,2,3.1,4,5,6),
#' arm = rep("Arm A", 6),
#' dropper = c(2L, rep(1L, 5)),
#' enroller = rep(3L, 6)
#' )
#'
#' add_timepoints(t, timepoints)
add_timepoints <- function(timer,df){
  invisible(
    sapply(
      split(df, 1:nrow(df)),
      function(x) do.call(timer$add_timepoint, x)
    )
  )
}


#' Print trial results as a data.frame
#'
#' @param results pass results data field of your `Trial`
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

#' Convert vector to our dataframe format
#'
#' @param data `numeric` population data in vector format
#'
#' @returns `data.frame` population data that may be passed to [Population]

#' @export
vector_to_dataframe<-function(data)
  if (is.vector(data)) {
    data <- data.frame(
      subject_id = seq_along(data),
      data = data
    )
    return(data)
  }

#' Generate timepoints for trial with allocation, piece-wise linear enrollment and dropout
#'
#' @param n `integer` trial sample size
#' @param arms `character` vector of unique identifier of arms
#' @param allocation `numeric` vector of arm allocations
#' @param enrollment `list` named list with enrollment parameters
#' - `end_time` vector of end times of duration
#' - `rate` # of subjects to enroll per unit time within that duration
#' @param dropout `list` named list with dropout parameters
#' - `end_time` vector of end times of duration
#' - `rate` # of subjects to enroll per unit time within that duration
#'
#' @returns `data.frame` timepoints that may be passed to [add_timepoints()]
#' @export
#'
#' @examples
#' gen_timepoints(
#'   n = 100,
#'   arms = c("A", "B"),
#'   allocation = c(2,1),
#'   enrollment = list(
#'     end_time = c(4,8,12),
#'     rate = c(6,12,18)
#'   ),
#'   dropout = list(
#'     end_time = c(5,9,13),
#'     rate = c(0,3,6)
#'   )
#' )
gen_timepoints <- function(n, arms, allocation, enrollment, dropout) {

  # determine number of arms
  n_arms <- length(arms)
  # get valid weights
  ratio <- allocation / sum(allocation)
  # get the end time
  end <- max(tail(enrollment$end_time,1), tail(dropout$end_time, 1))

  # pad with 0 rates when necessary
  pad <- function(x, end) if (tail(x$end_time, 1) != end) {
    x$end_time <- c(x$end_time, end)
    x$rate <- c(x$rate, 0)
    x
  } else x

  enrollment <- pad(enrollment, end)
  dropout <- pad(dropout, end)

  # get duration from end times
  get_durations <- function(x) (c(0,x) - dplyr::lag(c(0,x)))[-1]

  df <- data.frame(
    time = rep(seq_len(end), n_arms),
    arm = rep(arms, each=end),
    enroller = rep(
      as.vector(outer(enrollment$rate, ratio)),
      rep(get_durations(enrollment$end_time), n_arms)
    ) |> as.integer(),
    dropper = rep(
      as.vector(outer(dropout$rate, ratio)),
      rep(get_durations(dropout$end_time), n_arms)
    ) |> as.integer()
  )
}

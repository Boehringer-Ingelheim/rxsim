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




#' Convert vector to our data standards
#'
#' @param results pass results data field of your `data`
#'
#' @export

vector_to_dataframe<-function(data)
  if (is.vector(data)) {
    data <- data.frame(
      subject_id = seq_along(data),
      data = data
    )
    return(data)
  }

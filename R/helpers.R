add_timepoints <- function(timer,df){
  invisible(
    sapply(
      split(df, 1:nrow(df)),
      function(x) do.call(timer$add_timepoint, x)
    )
  )
}

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

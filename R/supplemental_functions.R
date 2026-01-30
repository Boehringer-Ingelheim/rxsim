add_timepoints <- function(timer,df){
  sapply(
    split(df, 1:nrow(df)),
    function(x) do.call(timer$add_timepoint, x)
  )
}

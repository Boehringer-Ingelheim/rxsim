add_timepoints<-function(timer,df){
  sapply(
    split(time_matrix, 1:nrow(df)),
    function(x) do.call(timer$add_timepoint, x)
  )
}

add_timepoints(timer=t,df=time_matrix)

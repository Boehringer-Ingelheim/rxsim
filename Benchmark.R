
start_time <- Sys.time()

# --- Two populations with a common 'value' column ---
# long_format
pop1 <- Population$new("Arm A", data = data.frame(
  subject_id = 1:20,
  value = rnorm(20, mean = 50)
))
pop2 <- Population$new("Arm B", data = data.frame(
  subject_id = 21:40,
  value = rnorm(20, mean = 55)
))


timepoints <- data.frame(
  time = c(1,2,3.1,4,5,6),
  arm = rep("Arm A", 6),
  dropper = c(2, rep(1, 5)),
  enroller = rep(3, 6)
)

# --- Timers with multiple timepoints ---
t <- Timer$new(name = "TrialTimers") # Use your updated Timers from earlier
add_timepoints(t, timepoints)

timepoints$time <- 1:6
timepoints$arm <- rep("Arm B", 6)
add_timepoints(t, timepoints)

# --- Reader conditions (per-reader predicates like dplyr::filter) ---
# IMPORTANT: func must accept (data, current_time). Wrap base functions accordingly.



# event condition
t$add_condition(
  length(value) > 4,
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "overall_mean"
)

#


# --- Trial with list of populations ---
trials<-list()
for(i in 1:10000)
  (
    trials[[i]] <- Trial$new(
      name = paste("Trial",as.character(i)),
      timer = t$clone(),
      population = list(pop1$clone(), pop2$clone())
    )
  )

lapply(trials,function(x) x$run())
#lapply(trials,function(x) print(x$results))


end_time <- Sys.time()
end_time - start_time








start_time <- Sys.time()

# --- Two populations with a common 'value' column ---
# long_format
pop1 <- Population$new("Arm A", data = data.frame(
  subject_id = 1:20,
  value = rnorm(20, mean = 50)
))
pop2 <- Population$new("Arm B", data = data.frame(
  subject_id = 21:40,
  value = rnorm(20, mean = 55)
))


timepoints <- data.frame(
  time = c(1,2,3.1,4,5,6),
  arm = rep("Arm A", 6),
  dropper = c(2, rep(1, 5)),
  enroller = rep(3, 6)
)

# --- Timers with multiple timepoints ---
t <- Timer$new(name = "TrialTimers") # Use your updated Timers from earlier
add_timepoints(t, timepoints)

timepoints$time <- 1:6
timepoints$arm <- rep("Arm B", 6)
add_timepoints(t, timepoints)

# --- Reader conditions (per-reader predicates like dplyr::filter) ---
# IMPORTANT: func must accept (data, current_time). Wrap base functions accordingly.



# event condition
t$add_condition(
  length(value) > 4,
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "overall_mean"
)

#


# --- Trial with list of populations ---
trials<-list()
for(i in 1:10000)
  (
    trials[[i]] <- Trial$new(
      name = paste("Trial",as.character(i)),
      timer = t$clone(),
      population = list(pop1$clone(), pop2$clone())
    )
  )

lapply(trials,function(x) x$run())
#lapply(trials,function(x) print(x$results))


end_time <- Sys.time()
end_time - start_time


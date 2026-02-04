# Record the start time to measure total runtime
start_time <- Sys.time()

# ----------------------------------------------------
# Create two populations with a shared 'value' column
# ----------------------------------------------------

# Population for Arm A with 20 subjects
pop1 <- Population$new(
  "Arm A",
  data = data.frame(
    subject_id = 1:20,               # Unique subject IDs
    value = rnorm(20, mean = 50)     # Random values with mean 50
  )
)

# Population for Arm B with 20 subjects
pop2 <- Population$new(
  "Arm B",
  data = data.frame(
    subject_id = 21:40,              # Non-overlapping subject IDs
    value = rnorm(20, mean = 55)     # Random values with mean 55
  )
)

# ----------------------------------------------------
# Define timepoints for Arm A
# ----------------------------------------------------

timepoints <- data.frame(
  time = c(1, 2, 3.1, 4, 5, 6),       # Scheduled analysis times
  arm = rep("Arm A", 6),             # Associated treatment arm
  dropper = c(2, rep(1, 5)),          # Dropout pattern
  enroller = rep(3, 6)               # Enrollment pattern
)

# ----------------------------------------------------
# Create timer object and add timepoints
# ----------------------------------------------------

t <- Timer$new(name = "TrialTimers")

# Add Arm A timepoints to the timer
add_timepoints(t, timepoints)

# Modify the same timepoints for Arm B
timepoints$time <- 1:6
timepoints$arm <- rep("Arm B", 6)

# Add Arm B timepoints to the timer
add_timepoints(t, timepoints)

# ----------------------------------------------------
# Define an event-based analysis condition
# ----------------------------------------------------

# This condition triggers when the dataset has > 4 values
# The analysis computes the mean of 'value'
t$add_condition(
  length(value) > 4,
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "overall_mean"
)

# ----------------------------------------------------
# Create and run 10,000 independent trial objects
# ----------------------------------------------------

trials <- list()

for (i in 1:10000) {
  trials[[i]] <- Trial$new(
    name = paste("Trial", i),        # Unique trial name
    timer = t$clone(),               # Clone timer to avoid shared state
    population = list(
      pop1$clone(),                  # Clone populations for isolation
      pop2$clone()
    )
  )
}

# Execute all trials (results stored internally)
lapply(trials, function(x) x$run())

# ----------------------------------------------------
# Record end time and compute elapsed runtime
# ----------------------------------------------------

end_time <- Sys.time()
end_time - start_time

# Observed runtime: ~1.67 minutes



# Record the start time to measure total runtime
start_time <- Sys.time()

# ----------------------------------------------------
# Create two populations
# ----------------------------------------------------

pop1 <- Population$new(
  "Arm A",
  data = data.frame(
    subject_id = 1:20,
    value = rnorm(20, mean = 50)
  )
)

pop2 <- Population$new(
  "Arm B",
  data = data.frame(
    subject_id = 21:40,
    value = rnorm(20, mean = 55)
  )
)

# ----------------------------------------------------
# Define timepoints for Arm A
# ----------------------------------------------------

timepoints <- data.frame(
  time = c(1, 2, 3.1, 4, 5, 6),
  arm = rep("Arm A", 6),
  dropper = c(2, rep(1, 5)),
  enroller = rep(3, 6)
)

# ----------------------------------------------------
# Create timer and add timepoints for both arms
# ----------------------------------------------------

t <- Timer$new(name = "TrialTimers")
add_timepoints(t, timepoints)

timepoints$time <- 1:6
timepoints$arm <- rep("Arm B", 6)
add_timepoints(t, timepoints)

# ----------------------------------------------------
# Add multiple analysis conditions
# Each condition increases runtime
# ----------------------------------------------------

# Condition 1: minimum number of observations
t$add_condition(
  length(value) > 4,
  analysis = function(d, tt) mean(d$value),
  name = "overall_mean"
)

# Condition 2: threshold on values
t$add_condition(
  sum(value > 40) > 10,
  analysis = function(d, tt) mean(d$value),
  name = "overall_mean"
)

# Condition 3: time-specific condition
t$add_condition(
  time == 4,
  analysis = function(d, tt) mean(d$value),
  name = "time_mean"
)

# Condition 4: combined value and time condition
t$add_condition(
  sum((value + time) > 45) > 10,
  analysis = function(d, tt) mean(d$value),
  name = "overall_mean_2"
)

# ----------------------------------------------------
# Create and run 10,000 trials
# ----------------------------------------------------

trials <- list()

for (i in 1:10000) {
  trials[[i]] <- Trial$new(
    name = paste("Trial", i),
    timer = t$clone(),
    population = list(
      pop1$clone(),
      pop2$clone()
    )
  )
}

# Run all trials
lapply(trials, function(x) x$run())

# ----------------------------------------------------
# Measure elapsed runtime
# ----------------------------------------------------

end_time <- Sys.time()
end_time - start_time

# Observed runtime: ~3.91 minutes

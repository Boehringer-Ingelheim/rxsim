#problem 1 we want to have in our out out is good so say matt and saumil
grid=parm1 cross parm2
# Population for Arm A with 20 subjects
pop1 <- Population$new(
  "Arm A",
  data = function(n) makeData(n,parm1,parm2)     # Random values with mean 5  )
)

function
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

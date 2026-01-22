data= rbind(cbind(subject_id=1:20,measurement_time=1,value = rnorm(20, mean = 50)),cbind(subject_id=1:20,measurement_time=4,value = rnorm(20, mean = 50)))
pop1 <- Population$new("Arm A", data = data.frame(
data
))




# --- Timers with multiple timepoints ---
t <- Timer$new(name = "TrialTimers") # Use your updated Timers from earlier
t$add_timepoint(time=1.0,arm="Arm A",dropper = 2, enroller = 3)
t$add_timepoint(time=2.0,arm="Arm A",dropper = 1, enroller = 3)
t$add_timepoint(time=3.1,arm="Arm A",dropper = 1, enroller = 3)
t$add_timepoint(time=4.0,arm="Arm A",dropper = 1, enroller = 3)
t$add_timepoint(time=5.0,arm="Arm A",dropper = 1, enroller = 3)
t$add_timepoint(time=6.0,arm="Arm A",dropper = 1, enroller = 3)

t$add_condition(
  length(value) > 4,
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "overall_mean"
)


trial <- Trial$new(
  name = "Trial A",
  timer = t,
  population = list(pop1)
)
self <- trial

# --- Run ---
trial$run()

trial$results

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

# --- Timers with multiple timepoints ---
t <- Timer$new(name = "TrialTimers") # Use your updated Timers from earlier
t$add_timepoint(time=1.0,arm="Arm A",dropper = 2, enroller = 3)
t$add_timepoint(time=2.0,arm="Arm A",dropper = 1, enroller = 3)
t$add_timepoint(time=3.1,arm="Arm A",dropper = 1, enroller = 3)
t$add_timepoint(time=4.0,arm="Arm A",dropper = 1, enroller = 3)
t$add_timepoint(time=5.0,arm="Arm A",dropper = 1, enroller = 3)
t$add_timepoint(time=6.0,arm="Arm A",dropper = 1, enroller = 3)


t$add_timepoint(time=1.0,arm="Arm B",dropper = 2, enroller = 3)
t$add_timepoint(time=2.0,arm="Arm B",dropper = 1, enroller = 3)
t$add_timepoint(time=3.0,arm="Arm B",dropper = 1, enroller = 3)
t$add_timepoint(time=4.0,arm="Arm B",dropper = 1, enroller = 3)
t$add_timepoint(time=5.0,arm="Arm B",dropper = 1, enroller = 3)
t$add_timepoint(time=6.0,arm="Arm B",dropper = 1, enroller = 3)

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

t$add_condition(
  sum(value > 40) > 10,
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "overall_mean"
)

t$add_condition(
  time == 4.0,
  # uses the 'time' column added by Trials$fit()
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "time_mean"
)
t$add_condition(
  sum((value + time) > 45) > 10,
  analysis = function(d, tt) {
    mean(d$value)
  },
  name = "overall_mean_2"
)

# --- Trial with list of populations ---
trial <- Trial$new(
  name = "Trial A",
  timer = t,
  population = list(pop1, pop2)
)
 self <- trial

# --- Run ---
trial$run()

trial$results


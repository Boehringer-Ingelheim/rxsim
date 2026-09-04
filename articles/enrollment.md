# Enrollment and Dropout

``` r
library(rxsim)
```

## Introduction: why enrollment timing matters

Enrollment pace determines when each subject enters the trial and
therefore when the trial clock can reach the timepoints that trigger
analyses. In most real programs, the rate at which sites activate and
patients are screened is uncertain - enrollment is itself a random
process. When that randomness is ignored (for example, by assuming all
subjects enroll on a fixed schedule), the simulated distribution of
study duration is artificially narrow, and operating characteristics
such as power or expected decision error rates can be materially
overestimated. Capturing stochastic enrollment propagates uncertainty in
study duration through to every downstream quantity the simulation is
designed to evaluate.

rxsim offers two complementary strategies, each suited to a different
modeling philosophy:

- **[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)** -
  random. You supply a function for inter-enrollment times. Every call
  draws a new realization, so each replicate gets its own unique
  enrollment timeline. Use this when trial-duration variability is
  substantively important and you want operating characteristics to
  reflect it.

- **[`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)** -
  piecewise-constant. You specify how many subjects enroll per unit time
  in each time period. The schedule is fixed. Every replicate uses
  exactly the same enrollment pattern. Use this when you have a
  well-characterized enrollment plan and want to isolate variability in
  endpoints and analyses from variability in timing.

Both functions return a data frame that is passed directly to
`Timer$add_schedule()` to populate a `Timer`.

## Stochastic enrollment with `stochastic_schedule()`

### Concept

In the stochastic approach, enrollment is modeled as a sequence of
waiting times between successive subject arrivals. You provide a
function that, given `n`, returns a vector of `n` independent
inter-arrival durations. rxsim then takes the cumulative sum of those
durations to produce the calendar time at which each subject enrolls.
This is the standard structure of a Poisson process - a simple and
widely used model for patient arrival - although any non-negative
distribution can be substituted.

For example, `function(n) rexp(n, rate = 1)` draws waiting times from an
exponential distribution with rate 1, which means the average gap
between consecutive enrollments is 1 time unit. Over the course of the
trial these gaps accumulate to produce enrollment calendar times that,
in expectation, grow linearly but vary around that line across
replicates.

### Example

Below we generate a plan for 12 subjects across two arms in balanced
allocation. The dropout function uses a much lower rate (0.05),
reflecting a trial where most subjects are expected to complete.

``` r
set.seed(101)

plan_a <- stochastic_schedule(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = function(n) rexp(n, rate = 1),
  dropout     = function(n) rexp(n, rate = 0.05)
)

plan_a
#>          time arm enroll drop
#> 1    1.181941 pbo      1    0
#> 2    4.776292 pbo      1    0
#> 3    6.685949 pbo      1    0
#> 4    7.175917 pbo      1    0
#> 5    9.545312 pbo      1    0
#> 6   10.094619 pbo      1    0
#> 7   28.718580 pbo      0    1
#> 8   39.586323 pbo      0    1
#> 9   82.273269 pbo      0    1
#> 10 122.125954 pbo      0    1
#> 11 122.874100 pbo      0    1
#> 12 134.037356 pbo      0    1
#> 13 146.620943 pbo      0    1
#> 14 159.022740 pbo      0    1
#> 15 163.420177 pbo      0    1
#> 16   4.356924 trt      1    0
#> 17   5.091672 trt      1    0
#> 18   7.491383 trt      1    0
#> 19   8.132255 trt      1    0
#> 20   9.946735 trt      1    0
#> 21  10.143980 trt      1    0
#> 22  89.464239 trt      0    1
#> 23  98.838488 trt      0    1
#> 24 130.506330 trt      0    1
```

### Output columns

Each row in the returned data frame represents a single event - either
an enrollment or a dropout:

| Column   | Meaning                                               |
|----------|-------------------------------------------------------|
| `time`   | Calendar time of the event                            |
| `arm`    | Arm label (`"pbo"` or `"trt"`)                        |
| `enroll` | `1` if this row is an enrollment event, `0` otherwise |
| `drop`   | `1` if this row is a dropout event, `0` otherwise     |

Because the enrollment times are cumulative sums of exponential random
variables, they are continuous-valued and always positive.

### Every call produces a different schedule

Because the inter-arrival times are drawn freshly each time, repeating
the call yields a different timeline. This is the core mechanism by
which
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
produces independent replicates - it calls
[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)
once per replicate.

``` r
set.seed(202)
plan_b <- stochastic_schedule(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = function(n) rexp(n, rate = 1),
  dropout     = function(n) rexp(n, rate = 0.05)
)

# Enrollment times differ across the two calls
head(plan_a, 4)
#>       time arm enroll drop
#> 1 1.181941 pbo      1    0
#> 2 4.776292 pbo      1    0
#> 3 6.685949 pbo      1    0
#> 4 7.175917 pbo      1    0
head(plan_b, 4)
#>        time arm enroll drop
#> 1  5.952201 pbo      1    0
#> 2  6.516808 pbo      1    0
#> 3  8.490974 pbo      1    0
#> 4 12.483522 pbo      1    0
```

### Passing the plan to a Timer

Once you have the data frame, pass it to `Timer$add_schedule()` to
register all events on a `Timer` object.

``` r
set.seed(404)
plan <- stochastic_schedule(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = function(n) rexp(n, rate = 1),
  dropout     = function(n) rexp(n, rate = 0.05)
)

tmr <- Timer$new(name = "stochastic_timer")
tmr$add_schedule(plan)

tmr$get_end_timepoint()   # last event time
#> [1] 197.1576
tmr$get_n_arms()          # number of unique arms
#> [1] 2
```

## Deterministic enrollment with `deterministic_schedule()`

### Concept

In many programs, enrollment is planned in phases: sites start slowly,
ramp up as more centers open, and may plateau near the end. The
piecewise-constant approach lets you encode this structure directly. You
supply a list with two elements:

- `end_time`: a numeric vector of period endpoints (time boundaries).
- `rate`: a numeric vector of enrollment rates (subjects per unit time)
  for each period.

rxsim interprets these as a step function. In the period from time 0 up
to `end_time[1]`, `rate[1]` subjects are enrolled per unit time; from
`end_time[1]` to `end_time[2]` the rate switches to `rate[2]`; and so
on. The schedule is deterministic - every call with the same inputs
returns the same data frame.

### Example

Here enrollment ramps up across three four-time-unit periods: 3, 6, and
9 subjects per unit time. Dropout rate is low and constant.

``` r
enrollment_schedule <- list(
  end_time = c(4, 8, 12),
  rate     = c(3, 6, 9)
)

dropout_schedule <- list(
  end_time = c(4, 8, 12),
  rate     = c(0, 1, 1)
)

tp <- deterministic_schedule(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = enrollment_schedule,
  dropout     = dropout_schedule
)

tp
#>    time arm enroll drop
#> 1     1 pbo      1    0
#> 2     2 pbo      1    0
#> 3     3 pbo      1    0
#> 4     4 pbo      1    0
#> 5     5 pbo      2    0
#> 6     1 trt      1    0
#> 7     2 trt      1    0
#> 8     3 trt      1    0
#> 9     4 trt      1    0
#> 10    5 trt      2    0
```

### Output columns

The output format uses the same column names as
[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md),
but the semantics differ:

| Column   | Meaning                                                |
|----------|--------------------------------------------------------|
| `time`   | Integer calendar time of the event                     |
| `arm`    | Arm label                                              |
| `enroll` | Count of subjects enrolling at this time (can be \> 1) |
| `drop`   | Count of subjects dropping out at this time            |

All times are positive integers.

### Every call returns the same schedule

Calling
[`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
twice with the same arguments produces identical output - confirming the
deterministic nature of this approach.

``` r
tp2 <- deterministic_schedule(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = enrollment_schedule,
  dropout     = dropout_schedule
)

identical(tp, tp2)
#> [1] TRUE
```

### Passing the schedule to a Timer

As with
[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md),
pass the result directly to `Timer$add_schedule()`.

``` r
tmr2 <- Timer$new(name = "deterministic_timer")
tmr2$add_schedule(tp)

tmr2$get_end_timepoint()
#> [1] 5
tmr2$get_unique_times()
#> [1] 1 2 3 4 5
```

## When to use which

Running both functions with the same target sample size and the same
two-arm design makes the structural differences clear.

``` r
set.seed(555)

# Stochastic
plan_stoch <- stochastic_schedule(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = function(n) rexp(n, rate = 1),
  dropout     = function(n) rexp(n, rate = 0.05)
)

# Deterministic
plan_det <- deterministic_schedule(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = list(end_time = c(4, 8, 12), rate = c(3, 6, 9)),
  dropout     = list(end_time = c(4, 8, 12), rate = c(0, 1, 1))
)

# stochastic_schedule: continuous (fractional) enrollment times
head(plan_stoch[plan_stoch$enroll == 1, ], 5)
#>       time arm enroll drop
#> 1 1.840486 pbo      1    0
#> 2 4.792944 pbo      1    0
#> 3 4.818969 pbo      1    0
#> 4 9.206899 pbo      1    0
#> 5 9.709935 pbo      1    0

# deterministic_schedule: integer times, counts >= 1
head(plan_det[plan_det$enroll > 0, ], 5)
#>   time arm enroll drop
#> 1    1 pbo      1    0
#> 2    2 pbo      1    0
#> 3    3 pbo      1    0
#> 4    4 pbo      1    0
#> 5    5 pbo      2    0
```

The key trade-offs:

|                             | [`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md) | [`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md) |
|-----------------------------|----------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------|
| Enrollment times            | Continuous (real-valued)                                                                                 | Integer                                                                                                        |
| Replicate-to-replicate      | Different each call                                                                                      | Identical every call                                                                                           |
| Study-duration distribution | Fully propagated                                                                                         | Fixed by design                                                                                                |
| Best for                    | Uncertainty in timing                                                                                    | Sensitivity to endpoints only                                                                                  |

For large-scale simulations where study duration is the main
uncertainty,
[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)
is the natural choice. When you need a reproducible baseline enrollment
profile - for example, to match a pre-specified operational plan or to
run sensitivity analyses against a fixed schedule -
[`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
is preferable.

## Dropout modeling

Dropout uses exactly the same two-approach structure as enrollment,
controlled by the `dropout` argument.

### Stochastic dropout with `stochastic_schedule()`

Pass a function to `dropout`. A low rate means most subjects complete
the study:

``` r
# rate = 0.05 => average time-to-dropout is 20 time units
dropout_fn <- function(n) rexp(n, rate = 0.05)
```

Each subject is assigned an independent dropout time drawn from this
distribution. Dropout events appear in the returned data frame as rows
where `drop == 1` and `enroll == 0`. In the locked data snapshot
produced during a trial run, subjects who have not yet dropped out have
`drop_time = NA`.

### Piecewise-constant dropout with `deterministic_schedule()`

Pass a list to `dropout` using the same `end_time` / `rate` structure:

``` r
# No dropouts in the first period; 1 per unit time thereafter
dropout_stepped <- list(
  end_time = c(4, 8, 12),
  rate     = c(0, 1,  1)
)
```

The dropout rate applies independently of enrollment: it is not
conditioned on how many subjects are currently enrolled. Both processes
are generated separately and then combined into a single schedule.

## The Timer class: `add_schedule`, `timelist`

A `Timer` is the trial clock. It stores a data frame of timepoint
records, where each record says when something happens in one arm and
how many subjects are affected. Each record has four fields: `time`,
`arm`, `enroll`, and `drop`.

`Trial$run()` does not step through subjects one by one. Instead, it
collects all distinct event times from the timer, sorts them, and
updates the trial at those times. That makes `Timer` the object that
defines the calendar of the simulation.

### Building a Timer manually with `add_schedule()`

`add_schedule()` takes a data frame with columns `time`, `arm`,
`enroll`, and `drop` - one row per event. That is the same shape
[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)
and
[`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
return, so hand-written and generated schedules load through the same
call.

Manual construction is useful when you want complete control over the
event schedule:

``` r
manual_timer <- Timer$new(name = "manual_timer")

manual_timer$add_schedule(data.frame(
  time   = c(1, 2, 4, 1, 3, 4),
  arm    = c("pbo", "pbo", "pbo", "trt", "trt", "trt"),
  drop   = c(0L, 1L, 1L, 0L, 0L, 1L),
  enroll = c(3L, 2L, 0L, 2L, 3L, 0L)
))

# Inspect the timelist structure
manual_timer$get_end_timepoint()
#> [1] 4
manual_timer$get_n_arms()
#> [1] 2
sort(manual_timer$get_unique_times())
#> [1] 1 2 3 4
subset(manual_timer$timelist, arm == "trt" & time == 3)
#>   time arm drop enroll
#> 5    3 trt    0      3
```

`enroll` and `drop` are subject counts, so they must be integers (`3L`,
not `3`); fractional counts are rejected. Calls accumulate, so you can
add events in several batches.

### Adding a generated schedule

The output of
[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)
or
[`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
goes straight in, with no reshaping:

``` r
set.seed(2024)
plan_stoch2 <- stochastic_schedule(
  sample_size = 8,
  arms = c("pbo", "trt"),
  allocation = c(1, 1),
  enrollment = function(n) rexp(n, rate = 1),
  dropout = function(n) rexp(n, rate = 0.15)
)

stochastic_timer <- Timer$new(name = "stochastic_timer2")
stochastic_timer$add_schedule(plan_stoch2)

head(stochastic_timer$get_unique_times())
#> [1] 2.405811 2.407575 3.388871 4.377373 4.734888 5.293186
stochastic_timer$get_end_timepoint()
#> [1] 87.41329
```

### The `timelist` structure

`timer$timelist` is a `data.frame` with columns `time`, `arm`, `enroll`,
and `drop`. Each row is one event registered through `add_schedule()`.

``` r
trial_clock <- Timer$new(name = "trial_clock")
trial_clock$add_schedule(data.frame(
  time   = c(1, 1, 4),
  arm    = c("pbo", "trt", "pbo"),
  drop   = c(0L, 0L, 1L),
  enroll = c(3L, 3L, 0L)
))

trial_clock$timelist
#>   time arm drop enroll
#> 1    1 pbo    0      3
#> 2    1 trt    0      3
#> 3    4 pbo    1      0
sort(unique(trial_clock$timelist$time))
#> [1] 1 4
```

## Putting it together: full example

This example builds a complete trial manually using
[`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md),
a `Timer`, a `Condition`, and two `Population` objects.

1.  Call
    [`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
    to get the schedule data frame.
2.  Create a `Timer$new()` and register timepoints with
    `add_schedule()`.
3.  Build `Condition` objects and pass them to
    `Trial$new(conditions = ...)`.
4.  Build populations and construct a `Trial$new()`.

``` r
# Fixed enrollment schedule
enroll_sched  <- list(end_time = c(4, 8, 12), rate = c(3, 6, 9))
dropout_sched <- list(end_time = c(4, 8, 12), rate = c(0, 1, 1))

fixed_plan <- deterministic_schedule(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = enroll_sched,
  dropout     = dropout_sched
)

tmr3 <- Timer$new(name = "fixed_timer")
tmr3$add_schedule(fixed_plan)

# Condition: fire at the final calendar time
final_t <- tmr3$get_end_timepoint()
cond_final <- condition_calendar_time(
  cal_time = final_t,
  analysis = function(df, current_time) {
    enrolled <- subset(df, !is.na(enroll_time))
    data.frame(n_enrolled = nrow(enrolled), time = current_time)
  },
  name = "final"
)

# Build populations (sized to match the fixed plan)
n_pbo <- sum(fixed_plan$enroll[fixed_plan$arm == "pbo"])
n_trt <- sum(fixed_plan$enroll[fixed_plan$arm == "trt"])

pop_pbo <- Population$new(
  name = "pbo",
  data = data.frame(
    id = seq_len(n_pbo),
    response = rnorm(n_pbo, 0),
    readout_time = 0
  )
)
pop_trt <- Population$new(
  name = "trt",
  data = data.frame(
    id = seq_len(n_trt),
    response = rnorm(n_trt, 0.5),
    readout_time = 0
  )
)

trial_fixed <- Trial$new(
  name       = "fixed_schedule_trial",
  timer      = tmr3,
  seed       = 123,
  population = list(pop_pbo, pop_trt),
  conditions = list(cond_final)
)

trial_fixed$run()
collect_results(trial_fixed)
#>   replicate timepoint analysis n_enrolled time
#> 1         1         5    final         12    5
```

Because every replicate cloned from `trial_fixed` uses the same
`fixed_plan`, enrollment timing is constant - only subject-level
endpoint variability contributes to simulation noise.

This same timer-centric pattern scales to seamless and platform-style
designs. You can stage arm activation by time, define multiple
`Condition` objects for interim and final looks, and carry decision
logic across looks inside the trial state. See [Example
5](https://boehringer-ingelheim.github.io/rxsim-gallery/examples/example-5.html)
for a complete seamless Ph2a/2b implementation.

## Next steps

- [Core
  Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.md) -
  detailed explanations of `Population`, `Condition`, `Trial`, and how
  they interact during a simulation run.
- [Conditions and
  Triggers](https://boehringer-ingelheim.github.io/rxsim/articles/conditions.md) -
  trigger API, interim analyses, `cooldown`, and `max_triggers`.
- [Example
  2](https://boehringer-ingelheim.github.io/rxsim-gallery/examples/example-2.html) -
  a worked example using
  [`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
  with two correlated continuous endpoints.

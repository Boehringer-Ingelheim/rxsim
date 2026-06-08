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
[`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md)
to populate a `Timer`.

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
#> 1    1.181941 trt      1    0
#> 2    4.356924 pbo      1    0
#> 3    4.776292 pbo      1    0
#> 4    5.091672 trt      1    0
#> 5    6.685949 trt      1    0
#> 6    7.175917 trt      1    0
#> 7    7.491383 pbo      1    0
#> 8    8.132255 pbo      1    0
#> 9    9.545312 trt      1    0
#> 10   9.946735 pbo      1    0
#> 11  10.094619 trt      1    0
#> 12  10.143980 pbo      1    0
#> 13  25.015565 pbo      0    1
#> 14  51.394606 pbo      0    1
#> 15  57.768809 pbo      0    1
#> 16  74.372883 pbo      0    1
#> 17 108.205985 pbo      0    1
#> 18 109.138428 trt      0    1
#> 19 112.004907 trt      0    1
#> 20 187.133638 pbo      0    1
#> 21 193.163489 pbo      0    1
#> 22 203.591067 trt      0    1
#> 23 224.090914 pbo      0    1
#> 24 249.502430 pbo      0    1
```

### Output columns

Each row in the returned data frame represents a single event - either
an enrollment or a dropout:

| Column        | Meaning                                                            |
|---------------|--------------------------------------------------------------------|
| `arm`         | Arm label (`"pbo"` or `"trt"`)                                     |
| `enroll_time` | Calendar time of enrollment (cumulative sum of inter-arrival gaps) |
| `drop_time`   | Calendar time of dropout; `NA` if the subject completes            |
| `enroll`      | `1` if this row is an enrollment event, `0` otherwise              |
| `drop`        | `1` if this row is a dropout event, `0` otherwise                  |

Because the times are cumulative sums of exponential random variables,
they are continuous-valued and always positive.

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
#> 1 1.181941 trt      1    0
#> 2 4.356924 pbo      1    0
#> 3 4.776292 pbo      1    0
#> 4 5.091672 trt      1    0
head(plan_b, 4)
#>       time arm enroll drop
#> 1 1.417195 trt      1    0
#> 2 1.757497 pbo      1    0
#> 3 2.769804 trt      1    0
#> 4 4.785916 pbo      1    0
```

### Passing the plan to a Timer

Once you have the data frame, pass it to
[`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md)
to register all events on a `Timer` object.

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
add_timepoints(tmr, plan)

tmr$get_end_timepoint()   # last event time
#> [1] 206.6357
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
#> # A tibble: 10 × 4
#>     time arm   enroll  drop
#>    <int> <chr>  <int> <int>
#>  1     1 pbo        1     0
#>  2     2 pbo        1     0
#>  3     3 pbo        1     0
#>  4     4 pbo        1     0
#>  5     5 pbo        2     0
#>  6     1 trt        1     0
#>  7     2 trt        1     0
#>  8     3 trt        1     0
#>  9     4 trt        1     0
#> 10     5 trt        2     0
```

### Output columns

The output format uses the same column names as
[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md),
but the semantics differ:

| Column        | Meaning                                                    |
|---------------|------------------------------------------------------------|
| `arm`         | Arm label                                                  |
| `enroll_time` | Integer calendar time at which one or more subjects enroll |
| `drop_time`   | Integer calendar time of dropout event                     |
| `enroll`      | Count of subjects enrolling at this time (can be \> 1)     |
| `drop`        | Count of subjects dropping out at this time                |

Because multiple subjects may enroll in the same time unit, `enroll` can
be greater than 1. All times are positive integers.

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
pass the result directly to
[`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md).

``` r
tmr2 <- Timer$new(name = "deterministic_timer")
add_timepoints(tmr2, tp)

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
#> 1 1.176158 trt      1    0
#> 2 1.548616 pbo      1    0
#> 3 1.840486 trt      1    0
#> 4 4.792944 pbo      1    0
#> 5 4.818969 pbo      1    0

# deterministic_schedule: integer times, counts >= 1
head(plan_det[plan_det$enroll > 0, ], 5)
#> # A tibble: 5 × 4
#>    time arm   enroll  drop
#>   <int> <chr>  <int> <int>
#> 1     1 pbo        1     0
#> 2     2 pbo        1     0
#> 3     3 pbo        1     0
#> 4     4 pbo        1     0
#> 5     5 pbo        2     0
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

## The Timer class: `add_timepoint`, `add_timepoints`, `timelist`

A `Timer` is the trial clock. It stores a list of timepoint records,
where each record says when something happens in one arm and how many
subjects are affected. Each record has four fields: `time`, `arm`,
`enroll`, and `drop`.

`Trial$run()` does not step through subjects one by one. Instead, it
collects all distinct event times from the timer, sorts them, and
updates the trial at those times. That makes `Timer` the object that
defines the calendar of the simulation.

### Building a Timer manually with `add_timepoint()`

Manual construction is useful when you want complete control over the
event schedule. You create a `Timer`, add one event at a time with
`add_timepoint()`, and inspect it with the query methods.

``` r
manual_timer <- Timer$new(name = "manual_timer")

manual_timer$add_timepoint(time = 1, arm = "pbo", drop = 0L, enroll = 3L)
manual_timer$add_timepoint(time = 2, arm = "pbo", drop = 1L, enroll = 2L)
manual_timer$add_timepoint(time = 4, arm = "pbo", drop = 1L, enroll = 0L)

manual_timer$add_timepoint(time = 1, arm = "trt", drop = 0L, enroll = 2L)
manual_timer$add_timepoint(time = 3, arm = "trt", drop = 0L, enroll = 3L)
manual_timer$add_timepoint(time = 4, arm = "trt", drop = 1L, enroll = 0L)

# Inspect the timelist structure
manual_timer$get_end_timepoint()
#> [1] 4
manual_timer$get_n_arms()
#> [1] 2
sort(manual_timer$get_unique_times())
#> [1] 1 2 3 4
manual_timer$get_timepoint("trt", 3)
#> $time
#> [1] 3
#> 
#> $arm
#> [1] "trt"
#> 
#> $drop
#> [1] 0
#> 
#> $enroll
#> [1] 3
```

### Adding a whole schedule with `add_timepoints()`

`add_timepoints(timer, plan_df)` is the bulk-loading function. It reads
every row of a schedule data frame produced by
[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)
or
[`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
and registers each row as a timepoint.

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
add_timepoints(stochastic_timer, plan_stoch2)

head(stochastic_timer$get_unique_times())
#> [1] 0.6738851 1.6505024 2.0112290 2.4075753 3.3888714 4.2948551
stochastic_timer$get_end_timepoint()
#> [1] 60.10328
```

### The `timelist` structure

`timer$timelist` is a named list of individual timepoint records. Each
element corresponds to one call to `add_timepoint()` and contains the
four fields `time`, `arm`, `enroll`, and `drop`.

``` r
trial_clock <- Timer$new(name = "trial_clock")
trial_clock$add_timepoint(time = 1, arm = "pbo", drop = 0L, enroll = 3L)
trial_clock$add_timepoint(time = 1, arm = "trt", drop = 0L, enroll = 3L)
trial_clock$add_timepoint(time = 4, arm = "pbo", drop = 1L, enroll = 0L)

trial_clock$timelist
#> [[1]]
#> [[1]]$time
#> [1] 1
#> 
#> [[1]]$arm
#> [1] "pbo"
#> 
#> [[1]]$drop
#> [1] 0
#> 
#> [[1]]$enroll
#> [1] 3
#> 
#> 
#> [[2]]
#> [[2]]$time
#> [1] 1
#> 
#> [[2]]$arm
#> [1] "trt"
#> 
#> [[2]]$drop
#> [1] 0
#> 
#> [[2]]$enroll
#> [1] 3
#> 
#> 
#> [[3]]
#> [[3]]$time
#> [1] 4
#> 
#> [[3]]$arm
#> [1] "pbo"
#> 
#> [[3]]$drop
#> [1] 1
#> 
#> [[3]]$enroll
#> [1] 0
sort(unique(vapply(trial_clock$timelist, function(x) x$time, numeric(1))))
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
    [`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md).
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
add_timepoints(tmr3, fixed_plan)

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
5](https://boehringer-ingelheim.github.io/rxsim/articles/example-5.md)
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
  2](https://boehringer-ingelheim.github.io/rxsim/articles/example-2.md) -
  a worked example using
  [`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
  with two correlated continuous endpoints.

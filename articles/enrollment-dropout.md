# Enrollment & Dropout Modeling

``` r
library(rxsim)
```

## Why enrollment modeling matters

Enrollment pace determines when each subject enters the trial and
therefore when the trial clock can reach the timepoints that trigger
analyses. In most real programs, the rate at which sites activate and
patients are screened is uncertain — enrollment is itself a random
process. When that randomness is ignored (for example, by assuming all
subjects enroll on a fixed schedule), the simulated distribution of
study duration is artificially narrow, and operating characteristics
such as power or expected decision error rates can be materially
overestimated. Capturing stochastic enrollment propagates uncertainty in
study duration through to every downstream quantity the simulation is
designed to evaluate.

## The two approaches

rxsim offers two complementary strategies, each suited to a different
modeling philosophy:

- **[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)**
  — stochastic. You supply a random function for inter-enrollment times.
  Every call draws a new realization, so each replicate gets its own
  unique enrollment timeline. Use this when trial-duration variability
  is substantively important and you want operating characteristics to
  reflect it.

- **[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)**
  — piecewise-constant. You specify how many subjects enroll per unit
  time in each time period. The schedule is fixed and deterministic.
  Every replicate uses exactly the same enrollment pattern. Use this
  when you have a well-characterized enrollment plan and want to isolate
  variability in endpoints and analyses from variability in timing.

Both functions return a data frame that is passed directly to
[`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md)
to populate a `Timer`.

## Stochastic enrollment with gen_plan()

### Concept

In the stochastic approach, enrollment is modeled as a sequence of
waiting times between successive subject arrivals. You provide a
function that, given `n`, returns a vector of `n` independent
inter-arrival durations. rxsim then takes the cumulative sum of those
durations to produce the calendar time at which each subject enrolls.
This is the standard structure of a Poisson process — a simple and
widely used model for patient arrival — although any non-negative
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

plan_a <- gen_plan(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = function(n) rexp(n, rate = 1),
  dropout     = function(n) rexp(n, rate = 0.05)
)

plan_a
#>          time arm enroller dropper
#> 1    1.181941 trt        1       0
#> 2    4.356924 pbo        1       0
#> 3    4.776292 pbo        1       0
#> 4    5.091672 trt        1       0
#> 5    6.685949 trt        1       0
#> 6    7.175917 trt        1       0
#> 7    7.491383 pbo        1       0
#> 8    8.132255 pbo        1       0
#> 9    9.545312 trt        1       0
#> 10   9.946735 pbo        1       0
#> 11  10.094619 trt        1       0
#> 12  10.143980 pbo        1       0
#> 13  25.015565 pbo        0       1
#> 14  51.394606 pbo        0       1
#> 15  57.768809 pbo        0       1
#> 16  74.372883 pbo        0       1
#> 17 108.205985 pbo        0       1
#> 18 109.138428 trt        0       1
#> 19 112.004907 trt        0       1
#> 20 187.133638 pbo        0       1
#> 21 193.163489 pbo        0       1
#> 22 203.591067 trt        0       1
#> 23 224.090914 pbo        0       1
#> 24 249.502430 pbo        0       1
```

### Output columns

Each row in the returned data frame represents a single event — either
an enrollment or a dropout:

| Column        | Meaning                                                            |
|---------------|--------------------------------------------------------------------|
| `arm`         | Arm label (`"pbo"` or `"trt"`)                                     |
| `enroll_time` | Calendar time of enrollment (cumulative sum of inter-arrival gaps) |
| `drop_time`   | Calendar time of dropout; `NA` if the subject completes            |
| `enroller`    | `1` if this row is an enrollment event, `0` otherwise              |
| `dropper`     | `1` if this row is a dropout event, `0` otherwise                  |

Because the times are cumulative sums of exponential random variables,
they are continuous-valued and always positive.

### Every call produces a different schedule

Because the inter-arrival times are drawn freshly each time, repeating
the call yields a different timeline. This is the core mechanism by
which
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
produces independent replicates — it calls
[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
once per replicate.

``` r
set.seed(202)
plan_b <- gen_plan(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = function(n) rexp(n, rate = 1),
  dropout     = function(n) rexp(n, rate = 0.05)
)

# Enrollment times differ across the two calls
head(plan_a, 4)
#>       time arm enroller dropper
#> 1 1.181941 trt        1       0
#> 2 4.356924 pbo        1       0
#> 3 4.776292 pbo        1       0
#> 4 5.091672 trt        1       0
head(plan_b, 4)
#>       time arm enroller dropper
#> 1 1.417195 trt        1       0
#> 2 1.757497 pbo        1       0
#> 3 2.769804 trt        1       0
#> 4 4.785916 pbo        1       0
```

### Passing the plan to a Timer

Once you have the data frame, pass it to
[`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md)
to register all events on a `Timer` object. The timer can then be
attached to a `Trial`.

``` r
set.seed(404)
plan <- gen_plan(
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

## Piecewise-constant enrollment with gen_timepoints()

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
on. The schedule is deterministic — every call with the same inputs
returns the same data frame. This makes it straightforward to verify the
exact enrollment profile being fed into the simulation.

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

tp <- gen_timepoints(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = enrollment_schedule,
  dropout     = dropout_schedule
)

tp
#> # A tibble: 10 × 4
#>     time arm   enroller dropper
#>    <dbl> <chr>    <int>   <int>
#>  1     1 pbo          1       0
#>  2     2 pbo          1       0
#>  3     3 pbo          1       0
#>  4     4 pbo          1       0
#>  5     5 pbo          2       0
#>  6     1 trt          1       0
#>  7     2 trt          1       0
#>  8     3 trt          1       0
#>  9     4 trt          1       0
#> 10     5 trt          2       0
```

### Output columns

The output format uses the same column names as
[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md),
but the semantics differ:

| Column        | Meaning                                                    |
|---------------|------------------------------------------------------------|
| `arm`         | Arm label                                                  |
| `enroll_time` | Integer calendar time at which one or more subjects enroll |
| `drop_time`   | Integer calendar time of dropout event                     |
| `enroller`    | Count of subjects enrolling at this time (can be \> 1)     |
| `dropper`     | Count of subjects dropping out at this time                |

Because multiple subjects may enroll in the same time unit, `enroller`
can be greater than 1. All times are positive integers.

### Every call returns the same schedule

Calling
[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
twice with the same arguments produces identical output — confirming the
deterministic nature of this approach.

``` r
tp2 <- gen_timepoints(
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
[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md),
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

## Side-by-side comparison

Running both functions with the same target sample size and the same
two-arm design makes the structural differences clear.

``` r
set.seed(555)

# Stochastic
plan_stoch <- gen_plan(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = function(n) rexp(n, rate = 1),
  dropout     = function(n) rexp(n, rate = 0.05)
)

# Deterministic
plan_det <- gen_timepoints(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = list(end_time = c(4, 8, 12), rate = c(3, 6, 9)),
  dropout     = list(end_time = c(4, 8, 12), rate = c(0, 1, 1))
)

# gen_plan: continuous (fractional) enrollment times
head(plan_stoch[plan_stoch$enroller == 1, ], 5)
#>       time arm enroller dropper
#> 1 1.176158 trt        1       0
#> 2 1.548616 pbo        1       0
#> 3 1.840486 trt        1       0
#> 4 4.792944 pbo        1       0
#> 5 4.818969 pbo        1       0

# gen_timepoints: integer times, counts >= 1
head(plan_det[plan_det$enroller > 0, ], 5)
#> # A tibble: 5 × 4
#>    time arm   enroller dropper
#>   <dbl> <chr>    <int>   <int>
#> 1     1 pbo          1       0
#> 2     2 pbo          1       0
#> 3     3 pbo          1       0
#> 4     4 pbo          1       0
#> 5     5 pbo          2       0
```

The key trade-offs:

|                             | [`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md) | [`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md) |
|-----------------------------|------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------|
| Enrollment times            | Continuous (real-valued)                                                           | Integer                                                                                        |
| Replicate-to-replicate      | Different each call                                                                | Identical every call                                                                           |
| Study-duration distribution | Fully propagated                                                                   | Fixed by design                                                                                |
| Best for                    | Uncertainty in timing                                                              | Sensitivity to endpoints only                                                                  |

For large-scale simulations where study duration is the main
uncertainty,
[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
is the natural choice. When you need a reproducible baseline enrollment
profile — for example, to match a pre-specified operational plan or to
run sensitivity analyses against a fixed schedule —
[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
is preferable.

## Dropout modeling

Dropout uses exactly the same two-approach structure as enrollment,
controlled by the `dropout` argument.

### Stochastic dropout with gen_plan()

Pass a function to `dropout`. A low rate means most subjects complete
the study:

``` r
# rate = 0.05 => average time-to-dropout is 20 time units
dropout_fn <- function(n) rexp(n, rate = 0.05)
```

Each subject is assigned an independent dropout time drawn from this
distribution. Dropout events appear in the returned data frame as rows
where `dropper == 1` and `enroller == 0`. In the locked data snapshot
produced during a trial run, subjects who have not yet dropped out have
`drop_time = NA`.

### Piecewise-constant dropout with gen_timepoints()

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

## Allocation ratios

The `allocation` argument is a numeric vector of relative weights — one
entry per arm. rxsim converts these to proportions internally:

``` r
# Three-arm trial: 2:1:1 allocation
plan_3arm <- gen_plan(
  sample_size = 12,
  arms        = c("pbo", "d1", "d2"),
  allocation  = c(2, 1, 1),
  enrollment  = function(n) rexp(n, rate = 1),
  dropout     = function(n) rexp(n, rate = 0.05)
)

# Per-arm enrollment counts
table(plan_3arm$arm[plan_3arm$enroller == 1])
#> 
#>  d1  d2 pbo 
#>   3   3   6
```

For three arms with `allocation = c(2, 1, 1)`, this yields proportions
0.50, 0.25, 0.25, giving a 2:1:1 split. The target per-arm sample sizes
are computed as `round(ratio * sample_size)`, and any rounding
discrepancy is corrected automatically by adding or removing subjects
from arms sampled proportionally to their weights — so the total always
equals `sample_size` exactly.

The same `allocation` argument works identically for
[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md).

## Using with replicate_trial()

When you call
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
with functions for `enrollment` and `dropout`, it calls
[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
internally once per replicate and constructs an independent `Timer` for
each trial. You do not need to call
[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
directly in this workflow.

``` r
set.seed(42)

population_generators <- list(
  pbo = function(n) data.frame(
    id = seq_len(n), response = rnorm(n, 0.0), readout_time = 0
  ),
  trt = function(n) data.frame(
    id = seq_len(n), response = rnorm(n, 0.5), readout_time = 0
  )
)

analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(sum(!is.na(enroll_time)) >= 12L),
    analysis = function(df, current_time) {
      enrolled <- subset(df, !is.na(enroll_time))
      data.frame(
        n         = nrow(enrolled),
        mean_pbo  = mean(enrolled$response[enrolled$arm == "pbo"]),
        mean_trt  = mean(enrolled$response[enrolled$arm == "trt"])
      )
    }
  )
)

trials <- replicate_trial(
  trial_name            = "demo",
  sample_size           = 12,
  arms                  = c("pbo", "trt"),
  allocation            = c(1, 1),
  enrollment            = function(n) rexp(n, rate = 1),
  dropout               = function(n) rexp(n, rate = 0.05),
  analysis_generators   = analysis_generators,
  population_generators = population_generators,
  n                     = 3
)

results <- run_trials(trials)

# Collect results across replicates
collect_results(trials)
#>    replicate  timepoint analysis  n   mean_pbo   mean_trt
#> 1          1   9.501693    final 12  0.1510848  0.4637581
#> 2          1  28.144979    final 12  0.1510848  0.4637581
#> 3          1  34.341294    final 12  0.1510848  0.4637581
#> 4          1  43.845119    final 12  0.1510848  0.4637581
#> 5          1  56.287325    final 12  0.1510848  0.4637581
#> 6          1  81.198931    final 12  0.1510848  0.4637581
#> 7          1  88.605720    final 12  0.1510848  0.4637581
#> 8          1 185.861820    final 12  0.1510848  0.4637581
#> 9          1 199.178463    final 12  0.1510848  0.4637581
#> 10         1 299.097833    final 12  0.1510848  0.4637581
#> 11         1 303.568979    final 12  0.1510848  0.4637581
#> 12         2  19.250321    final 12 -0.2867542 -0.1167607
#> 13         2  22.585208    final 12 -0.2867542 -0.1167607
#> 14         2  28.410483    final 12 -0.2867542 -0.1167607
#> 15         2  39.443418    final 12 -0.2867542 -0.1167607
#> 16         2  41.989291    final 12 -0.2867542 -0.1167607
#> 17         2  70.902014    final 12 -0.2867542 -0.1167607
#> 18         2  89.182427    final 12 -0.2867542 -0.1167607
#> 19         2  95.879488    final 12 -0.2867542 -0.1167607
#> 20         2 265.126204    final 12 -0.2867542 -0.1167607
#> 21         2 274.489977    final 12 -0.2867542 -0.1167607
#> 22         2 294.998679    final 12 -0.2867542 -0.1167607
#> 23         2 295.601212    final 12 -0.2867542 -0.1167607
#> 24         2 305.360198    final 12 -0.2867542 -0.1167607
#> 25         3  12.210078    final 12 -0.5907032  0.6425829
#> 26         3  28.564963    final 12 -0.5907032  0.6425829
#> 27         3  39.717343    final 12 -0.5907032  0.6425829
#> 28         3  65.135567    final 12 -0.5907032  0.6425829
#> 29         3  72.279282    final 12 -0.5907032  0.6425829
#> 30         3  83.312284    final 12 -0.5907032  0.6425829
#> 31         3 121.097218    final 12 -0.5907032  0.6425829
#> 32         3 192.293329    final 12 -0.5907032  0.6425829
#> 33         3 199.924797    final 12 -0.5907032  0.6425829
#> 34         3 239.169802    final 12 -0.5907032  0.6425829
#> 35         3 240.579267    final 12 -0.5907032  0.6425829
#> 36         3 277.034072    final 12 -0.5907032  0.6425829
```

Each element of `trials` holds an independent `Timer` built from a
freshly drawn enrollment plan, so study duration differs across
replicates.

## Building a Timer manually with gen_timepoints()

When you need a fixed enrollment schedule and want more control than
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
exposes, you can assemble the `Timer` directly. The pattern is:

1.  Call
    [`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
    to get the schedule data frame.
2.  Create a `Timer$new()`.
3.  Register timepoints with
    [`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md).
4.  Add analysis triggers (e.g., with
    [`trigger_by_calendar()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_calendar.md)).
5.  Build populations and construct a `Trial$new()`.
6.  Use
    [`clone_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/clone_trial.md)
    if you want multiple independent replicates sharing the same fixed
    schedule.

``` r
# Fixed enrollment schedule
enroll_sched  <- list(end_time = c(4, 8, 12), rate = c(3, 6, 9))
dropout_sched <- list(end_time = c(4, 8, 12), rate = c(0, 1, 1))

fixed_plan <- gen_timepoints(
  sample_size = 12,
  arms        = c("pbo", "trt"),
  allocation  = c(1, 1),
  enrollment  = enroll_sched,
  dropout     = dropout_sched
)

tmr3 <- Timer$new(name = "fixed_timer")
add_timepoints(tmr3, fixed_plan)

# Add an analysis at the final calendar time
final_t <- tmr3$get_end_timepoint()
trigger_by_calendar(
  cal_time = final_t,
  timer    = tmr3,
  analysis = function(df, current_time) {
    enrolled <- subset(df, !is.na(enroll_time))
    data.frame(n_enrolled = nrow(enrolled), time = current_time)
  }
)

# Build populations (sized to match the fixed plan)
pop_pbo <- gen_population(
  name        = "pbo",
  generator   = function(n) data.frame(
    id = seq_len(n), response = rnorm(n, 0), readout_time = 0
  ),
  sample_size = 6
)
pop_trt <- gen_population(
  name        = "trt",
  generator   = function(n) data.frame(
    id = seq_len(n), response = rnorm(n, 0.5), readout_time = 0
  ),
  sample_size = 6
)

trial_fixed <- Trial$new(
  name       = "fixed_schedule_trial",
  timer      = tmr3,
  seed       = 123,
  population = list(pop_pbo, pop_trt)
)

trial_fixed$run()
prettify_results(trial_fixed$results)
#>   time cal_time_5.n_enrolled cal_time_5.time
#> 1    5                    12               5
```

Because every replicate cloned from `trial_fixed` uses the same
`fixed_plan`, enrollment timing is constant — only subject-level
endpoint variability contributes to simulation noise.

## Next steps

- **[Getting
  Started](https://boehringer-ingelheim.github.io/rxsim/articles/getting-started.md)**
  — end-to-end walkthrough of a complete simulation using stochastic
  enrollment via
  [`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md).
- **[Core
  Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.md)**
  — detailed explanations of `Timer`, `Population`, `Trial`, and how
  they interact during a simulation run.
- **[Example
  3](https://boehringer-ingelheim.github.io/rxsim/articles/example-3.md)**
  — a worked example using
  [`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
  with two correlated continuous endpoints and a piecewise-linear
  enrollment ramp.

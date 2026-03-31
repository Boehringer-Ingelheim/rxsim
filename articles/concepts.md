# Core Concepts

## Overview

rxsim organises a clinical trial simulation around three collaborating
objects. A `Population` owns the subject-level data and tracks each
subject’s enrollment and dropout times. A `Timer` drives the trial
clock. It stores discrete timepoints per arm and evaluates
condition-triggered analyses when the data meet a criterion. A `Trial`
orchestrates the simulation by iterating over timepoints, updating
populations, snapshotting the enrolled cohort, and collecting results.

``` mermaid
graph LR
  P1(Control) --> TR(Trial)
  P2(Treatment) --> TR
  TI(Timer) --> TR
  TR --> LD(Locked Data)
  TR --> RS(Results)
#> <div class="mermaid">
#> graph LR
#>   P1(Control) --> TR(Trial)
#>   P2(Treatment) --> TR
#>   TI(Timer) --> TR
#>   TR --> LD(Locked Data)
#>   TR --> RS(Results)
#> </div>
```

The `run()` loop that `Trial` executes at each timepoint follows this
sequence:

``` mermaid
graph TD
  A([Next timepoint]) --> B[Enroll and drop subjects]
  B --> C[Snapshot enrolled subjects]
  C --> D{Condition triggered?}
  D -- Yes --> E[Run analysis]
  D -- No --> F[Advance clock]
  E --> F
  F --> G{More timepoints?}
  G -- Yes --> A
  G -- No --> H([Done])
#> <div class="mermaid">
#> graph TD
#>   A([Next timepoint]) --> B[Enroll and drop subjects]
#>   B --> C[Snapshot enrolled subjects]
#>   C --> D{Condition triggered?}
#>   D -- Yes --> E[Run analysis]
#>   D -- No --> F[Advance clock]
#>   E --> F
#>   F --> G{More timepoints?}
#>   G -- Yes --> A
#>   G -- No --> H([Done])
#> </div>
```

In most workflows you will never construct these objects by hand.
Instead you use the high-level entry point
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md) +
[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md),
which build and execute `n` independent `Trial` objects from your
generator functions. Understanding the three classes directly is useful
when you want to:

- inspect the locked snapshot mid-simulation for debugging
- write custom multi-timepoint designs that
  [`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
  cannot express
- use `Trial$new()` directly for a one-off single-run simulation (as in
  [Example
  8](https://boehringer-ingelheim.github.io/rxsim/articles/example-8.md))

The rest of this vignette unpacks each building block in turn.

## Population

A `Population` object represents a single arm’s worth of subjects. It
stores the endpoint data alongside two tracking vectors (`enrolled`,
`dropped`) that record the calendar time each subject was enrolled or
dropped — or `NA` if neither event has occurred yet.

### Data structure

The `data` field must be a data frame with at least four columns:

| Column         | Type      | Meaning                                                 |
|----------------|-----------|---------------------------------------------------------|
| `id`           | integer   | Unique subject identifier within the arm                |
| `arm`          | character | Arm label (auto-filled from `name` if absent)           |
| `readout_time` | numeric   | Lag from enrollment to endpoint observation             |
| *(endpoint)*   | numeric   | At least one endpoint column (e.g., `data`, `y`, `tte`) |

The `readout_time` column is not a calendar time — it is the lag between
a subject’s enrollment date and the date their endpoint is observed. Use
`0` for immediately available data (e.g., a baseline measure) and a
positive value for delayed readouts (e.g., `readout_time = 12` means the
endpoint is read 12 weeks after enrollment).

### vector_to_dataframe()

For the common case of a single continuous endpoint,
[`vector_to_dataframe()`](https://boehringer-ingelheim.github.io/rxsim/reference/vector_to_dataframe.md)
wraps a numeric vector into the required format with columns `id`,
`data`, and `readout_time = 0`:

``` r
ep <- rnorm(6, mean = 1, sd = 0.5)
df <- vector_to_dataframe(ep)
df
#>   id       data readout_time
#> 1  1  0.2999782            0
#> 2  2  1.1276585            0
#> 3  3 -0.2186318            0
#> 4  4  0.9972144            0
#> 5  5  1.3107764            0
#> 6  6  1.5742058            0
```

You can also build the data frame yourself — useful when you have
covariates, time-varying endpoints, or different readout lags per
subject.

### Repeated measurements and n_readouts

When a subject contributes multiple measurements (a longitudinal or
pharmacokinetic design), each measurement appears as its own row with a
distinct `readout_time`. `Population` detects this automatically and
stores the number of rows per subject in `n_readouts`:

``` r
# Two timepoints per subject: baseline (0) and week 12 (12)
long_df <- data.frame(
  id           = rep(1:4, each = 2),
  readout_time = rep(c(0, 12), times = 4),
  response     = rnorm(8)
)

pop_long <- Population$new(name = "treatment", data = long_df)
pop_long$n          # 4 unique subjects
#> [1] 4
pop_long$n_readouts # 2 rows per subject
#> [1] 2
```

### Enrollment and dropout vectors

`enrolled` and `dropped` are numeric vectors of length `n`, initialised
to `NA` for every subject. `NA` means “not yet acted on”: the subject
exists in the pool but has not been enrolled or dropped. The `Trial`
assigns calendar times to these vectors as the simulation clock
advances.

`set_enrolled(n, time)` marks `n` randomly chosen currently unenrolled
subjects as enrolled at `time`. `set_dropped(n, time)` similarly marks
`n` randomly chosen enrolled-and-not-yet-dropped subjects as dropped at
`time`.

``` r
set.seed(1)
pop <- Population$new(name = "control", data = vector_to_dataframe(rnorm(8)))
pop$enrolled  # all NA to start
#> [1] NA NA NA NA NA NA NA NA

pop$set_enrolled(5, time = 2)
pop$enrolled  # 5 subjects enrolled at t=2, 3 still NA
#> [1]  2  2  2 NA NA  2 NA  2

pop$set_dropped(2, time = 7)
pop$dropped   # 2 of the enrolled subjects dropped at t=7
#> [1]  7 NA NA NA NA NA NA  7
```

The `Trial`’s `run()` method calls these setters automatically based on
the `Timer`’s schedule — you rarely need to call them directly.

## Timer

A `Timer` drives the trial clock. It holds two data structures: a
`timelist` that specifies how many subjects to enroll or drop in each
arm at each time unit, and a `conditions` list that defines
analysis-triggering criteria. At each unique time in the `timelist`,
`Trial$run()` asks the `Timer` to evaluate its conditions against the
current snapshot.

### Timepoints

Each entry in the `timelist` has four fields:

| Field      | Type      | Meaning                                   |
|------------|-----------|-------------------------------------------|
| `time`     | numeric   | Calendar time of the event                |
| `arm`      | character | Which arm the event applies to            |
| `enroller` | integer   | Number of subjects to enroll at this time |
| `dropper`  | integer   | Number of subjects to drop at this time   |

Timepoints are per-arm because the two arms in a trial may enroll
subjects on different schedules, and enrollment events in one arm do not
affect the other.

``` r
t <- Timer$new(name = "my_timer")

t$add_timepoint(time = 1,  arm = "control",   enroller = 5L, dropper = 0L)
t$add_timepoint(time = 1,  arm = "treatment", enroller = 5L, dropper = 0L)
t$add_timepoint(time = 2,  arm = "control",   enroller = 3L, dropper = 1L)
t$add_timepoint(time = 2,  arm = "treatment", enroller = 4L, dropper = 0L)
t$add_timepoint(time = 10, arm = "control",   enroller = 0L, dropper = 0L)
t$add_timepoint(time = 10, arm = "treatment", enroller = 0L, dropper = 0L)

t$get_unique_times()   # c(1, 2, 10)
#> [1]  1  2 10
t$get_end_timepoint()  # 10
#> [1] 10
t$get_n_arms()         # 2
#> [1] 2
```

Rather than adding timepoints one by one,
[`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md)
accepts a data frame with the four columns above — exactly what
[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
and
[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
return. See [Enrollment & Dropout
Modeling](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment-dropout.md)
for details.

### Conditions

A condition pairs a filter expression with an optional analysis
function. When the `Trial` calls `check_conditions()` at each timepoint,
each condition’s expression is evaluated against the current snapshot
using
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).
If the filtered result is non-empty (the condition is satisfied), the
associated analysis function is called and its return value is stored.

Conditions are added with `add_condition()`. The `...` arguments are
[`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html)-style
boolean expressions, evaluated lazily against the snapshot columns:

``` r
# A toy snapshot data frame
snapshot <- data.frame(
  id          = 1:8,
  arm         = rep(c("control", "treatment"), 4),
  enroll_time = c(1, 1, 2, 2, NA, NA, NA, NA),
  data        = rnorm(8)
)

t2 <- Timer$new(name = "demo")

# Fire when at least 4 subjects are enrolled
t2$add_condition(
  sum(!is.na(enroll_time)) >= 4,
  analysis = function(df, current_time) {
    data.frame(n_enrolled = sum(!is.na(df$enroll_time)),
               fired_at   = current_time)
  },
  name = "interim"
)

res <- t2$check_conditions(locked_data = snapshot, current_time = 5)
res[["interim"]]
#>   n_enrolled fired_at
#> 1          4        5
```

If no analysis function is provided, `check_conditions()` returns the
filtered subset instead of calling a function — convenient for
inspection or debugging.

## Trial

A `Trial` wires one or more `Population`s and a `Timer` together and
runs the simulation.

### Constructor

``` r
Trial$new(
  name       = "my_trial",
  seed       = 42,           # optional; set for reproducibility
  timer      = my_timer,
  population = list(pop_a, pop_b)
)
```

The `population` argument is a list of `Population` objects, one per
arm. The arm labels in each `Population`’s `name` field must match the
`arm` identifiers in the `Timer`’s timelist.

### What run() does

Calling `trial$run()` executes the following loop:

1.  Iterate over each unique time `t` in the `Timer`’s timelist, in
    ascending order.
2.  For each arm, apply `set_enrolled()` and `set_dropped()` according
    to that arm’s timepoint entry.
3.  Build a combined snapshot by row-binding all arms’ enrolled
    subjects. Four columns are appended by `Trial`:

| Column             | Meaning                                                           |
|--------------------|-------------------------------------------------------------------|
| `enroll_time`      | Calendar time the subject was enrolled (`NA` if not yet enrolled) |
| `drop_time`        | Calendar time the subject dropped out (`NA` if still active)      |
| `measurement_time` | `enroll_time + readout_time`                                      |
| `time`             | The current clock time `t`                                        |

4.  Call `timer$check_conditions(snapshot, t)` — evaluate all conditions
    and collect analysis results.
5.  If any condition fired, store the snapshot in
    `locked_data[["time_t"]]` and the analysis outputs in
    `results[["time_t"]]`.

### locked_data and results

`locked_data` is a named list of snapshots, one per unique timepoint at
which at least one analysis fired. Each snapshot is the full
subject-level data frame at that moment — useful for auditing which
subjects were enrolled, computing post-hoc statistics, or debugging an
analysis function.

`results` has the same time-indexed structure, but each element is
itself a named list — one entry per condition that fired at that
timepoint. The value is whatever your analysis function returned.

### A minimal end-to-end example

``` r
set.seed(7)

# Two small populations
pop_a <- Population$new("A", data = vector_to_dataframe(rnorm(10)))
pop_b <- Population$new("B", data = vector_to_dataframe(rnorm(10, mean = 0.5)))

# Timer: enroll in two waves, final readout at time 15
tm <- Timer$new("tm")
tm$add_timepoint(time = 1,  arm = "A", enroller = 5L, dropper = 0L)
tm$add_timepoint(time = 1,  arm = "B", enroller = 5L, dropper = 0L)
tm$add_timepoint(time = 3,  arm = "A", enroller = 5L, dropper = 0L)
tm$add_timepoint(time = 3,  arm = "B", enroller = 5L, dropper = 0L)
tm$add_timepoint(time = 15, arm = "A", enroller = 0L, dropper = 0L)
tm$add_timepoint(time = 15, arm = "B", enroller = 0L, dropper = 0L)

# Trigger: fire at calendar time 15
trigger_by_calendar(15, tm, analysis = function(df, current_time) {
  enrolled <- subset(df, !is.na(enroll_time))
  data.frame(
    n       = nrow(enrolled),
    mean_A  = mean(enrolled$data[enrolled$arm == "A"]),
    mean_B  = mean(enrolled$data[enrolled$arm == "B"])
  )
})

trial <- Trial$new(name = "demo", seed = 7, timer = tm,
                   population = list(pop_a, pop_b))
trial$run()

# Inspect results
prettify_results(trial$results)
#>   time cal_time_15.n cal_time_15.mean_A cal_time_15.mean_B
#> 1   15            20          0.1039757           1.282517
```

``` r
# How many subjects were enrolled at each timepoint?
sapply(trial$locked_data, function(snap) sum(!is.na(snap$enroll_time)))
#> time_15 
#>      20
```

## Analysis triggers in depth

Triggers are the mechanism by which rxsim knows when to analyse the data
and what to compute. Understanding how they are stored and evaluated is
key to writing correct and flexible simulation code.

### Why rlang::exprs()?

A trigger condition needs to be a *stored* expression, not an
immediately evaluated one. When you write
`sum(!is.na(enroll_time)) >= 20`, R would evaluate it at the point of
definition — before any data exists.
[`rlang::exprs()`](https://rlang.r-lib.org/reference/defusing-advanced.html)
wraps the expression in a quoted form so it is only evaluated later,
inside `check_conditions()`, against the actual snapshot:

``` r
# This is a stored expression, not an evaluated boolean:
rlang::exprs(sum(!is.na(enroll_time)) >= 20)
```

### The !! (bang-bang) operator

When you want to inject a value from the current R environment into a
stored expression, use `!!`. Without it, the variable name is treated as
a column name in the snapshot data frame — almost certainly not what you
want:

``` r
target_n <- 40

# CORRECT: !! injects the value 40 at definition time
trigger <- rlang::exprs(sum(!is.na(enroll_time)) >= !!target_n)

# WRONG: "target_n" would be looked for as a column name in the snapshot
trigger_bad <- rlang::exprs(sum(!is.na(enroll_time)) >= target_n)
```

This distinction matters when you loop over scenarios with different
sample sizes: each iteration should bake in its own `target_n` value via
`!!`.

### Columns available in a trigger expression

The trigger expression is evaluated against the snapshot data frame,
which contains all columns from the `Population`’s `data` plus the four
columns appended by `Trial`:

- `enroll_time` — calendar enrollment time (`NA` if not enrolled)
- `drop_time` — calendar dropout time (`NA` if not dropped)
- `measurement_time` — `enroll_time + readout_time`
- `time` — current clock time
- `arm` — arm label

Any user-defined endpoint column (e.g., `data`, `response`, `tte`) is
also available.

### trigger_by_calendar() and trigger_by_fraction()

These two helpers wrap the most common trigger patterns:

`trigger_by_calendar(cal_time, timer, analysis)` fires at a specific
calendar time — useful for a pre-planned final analysis or a scheduled
interim review:

``` r
trigger_by_calendar(24, timer = tm, analysis = function(df, current_time) {
  data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
})
```

`trigger_by_fraction(fraction, timer, sample_size, analysis)` fires once
a given fraction of the total target sample has enrolled — natural for
information-fraction-based interim rules:

``` r
# Interim at 50% enrollment
trigger_by_fraction(0.5, timer = tm, sample_size = 100,
                    analysis = function(df, current_time) {
                      enrolled <- subset(df, !is.na(enroll_time))
                      data.frame(n = nrow(enrolled), time = current_time)
                    })
```

### Custom conditions

You can condition on any expression involving the snapshot columns. For
time-to-event endpoints this commonly means waiting for a target number
of events rather than a target enrollment count:

``` r
# Fire when 30 events have been observed (event = 1, censored = 0)
timer$add_condition(
  sum(event == 1 & !is.na(enroll_time)) >= !!n_events,
  analysis = my_tte_analysis,
  name = "event_driven_interim"
)
```

Multiple conditions can be combined with commas (they are ANDed
together, just as in
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)):

``` r
# Fire only for the treatment arm when 15 subjects are enrolled there
timer$add_condition(
  arm == "treatment",
  sum(!is.na(enroll_time[arm == "treatment"])) >= 15,
  analysis = my_analysis,
  name = "trt_interim"
)
```

### The analysis function signature

Every analysis function receives two arguments:

- `df`: the **full snapshot** — all arms, all currently enrolled
  subjects. The condition filter determined *whether* to fire; the
  analysis function decides *what to compute* from the full data.
- `current_time`: the numeric clock time at which the condition fired.

``` r
my_analysis <- function(df, current_time) {
  enrolled <- subset(df, !is.na(enroll_time))
  data.frame(
    fired_at  = current_time,
    n_ctrl    = sum(enrolled$arm == "control"),
    n_trt     = sum(enrolled$arm == "treatment"),
    mean_diff = mean(enrolled$data[enrolled$arm == "treatment"]) -
                mean(enrolled$data[enrolled$arm == "control"])
  )
}
```

### Return values

The analysis function’s return value is stored verbatim in `results`.
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
handles three cases:

- **`data.frame`** (recommended) — the standard pattern; one row per
  trigger event is the convention.
- **named `list`** — coerced to a single-row data frame.
- **`NULL`** — silently skipped by
  [`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md).

Returning `NULL` is useful for side-effects-only analyses (e.g., writing
a log entry) or for marking that a trigger condition was not yet
satisfied.

## Putting it together: replicate_trial() and run_trials()

Building a single `Trial` by hand is useful for exploration, but the
purpose of rxsim is operating-characteristic simulation across many
replicates.
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md) +
[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md)
handle this ergonomically.

### How replicate_trial() works

For each of the `n` replicates,
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md):

1.  Calls
    [`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
    with your `enrollment` and `dropout` functions to generate a fresh,
    stochastic enrollment/dropout schedule.
2.  Builds a new `Timer`, loads the schedule via
    [`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md),
    and registers all analysis triggers.
3.  Calls each population generator function to draw fresh subject-level
    endpoint data.
4.  Constructs a `Trial$new()` wiring the new `Timer` and `Population`
    objects together.

Each replicate therefore has independent endpoint data *and* independent
enrollment timing — both sources of variability that operating
characteristic simulations are designed to characterise.

### A complete six-step example

``` r
set.seed(42)

# Step 1 — design parameters
sample_size <- 30
arms        <- c("control", "treatment")
allocation  <- c(1, 1)
true_delta  <- 0.5

scenario <- data.frame(sample_size = sample_size, true_delta = true_delta)

# Step 2 — population generators (called fresh per replicate)
population_generators <- list(
  control   = function(n) vector_to_dataframe(rnorm(n)),
  treatment = function(n) vector_to_dataframe(rnorm(n, mean = true_delta))
)

# Step 3 — enrollment and dropout (inter-arrival functions)
enrollment <- function(n) rexp(n, rate = 2)
dropout    <- function(n) rexp(n, rate = 0.02)

# Step 4 — analysis trigger: fire when full enrollment is reached
analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, current_time) {
      enrolled <- subset(df, !is.na(enroll_time))
      data.frame(
        scenario,
        fired_at  = current_time,
        n_ctrl    = sum(enrolled$arm == "control"),
        n_trt     = sum(enrolled$arm == "treatment"),
        mean_ctrl = mean(enrolled$data[enrolled$arm == "control"]),
        mean_trt  = mean(enrolled$data[enrolled$arm == "treatment"])
      )
    }
  )
)

# Step 5 — create and run replicates
trials <- replicate_trial(
  trial_name            = "concepts_demo",
  sample_size           = sample_size,
  arms                  = arms,
  allocation            = allocation,
  enrollment            = enrollment,
  dropout               = dropout,
  analysis_generators   = analysis_generators,
  population_generators = population_generators,
  n                     = 10
)

run_trials(trials)
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_1
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[2]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_2
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[3]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_3
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[4]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_4
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[5]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_5
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[6]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_6
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[7]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_7
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[8]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_8
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[9]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_9
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[10]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: concepts_demo_10
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6

# Step 6 — collect and inspect
results <- collect_results(trials)
results
#>     replicate  timepoint analysis sample_size true_delta   fired_at n_ctrl
#> 1           1   14.85864    final          30        0.5   14.85864     15
#> 2           1   28.46928    final          30        0.5   28.46928     15
#> 3           1  179.39710    final          30        0.5  179.39710     15
#> 4           1  204.27664    final          30        0.5  204.27664     15
#> 5           1  222.00432    final          30        0.5  222.00432     15
#> 6           1  309.82477    final          30        0.5  309.82477     15
#> 7           1  346.69972    final          30        0.5  346.69972     15
#> 8           1  348.14102    final          30        0.5  348.14102     15
#> 9           1  365.70174    final          30        0.5  365.70174     15
#> 10          1  385.33474    final          30        0.5  385.33474     15
#> 11          1  480.43903    final          30        0.5  480.43903     15
#> 12          1  519.35371    final          30        0.5  519.35371     15
#> 13          1  552.16956    final          30        0.5  552.16956     15
#> 14          1  571.49004    final          30        0.5  571.49004     15
#> 15          1  642.29464    final          30        0.5  642.29464     15
#> 16          1  718.45108    final          30        0.5  718.45108     15
#> 17          1  740.38666    final          30        0.5  740.38666     15
#> 18          1  948.79314    final          30        0.5  948.79314     15
#> 19          1 1008.54849    final          30        0.5 1008.54849     15
#> 20          1 1009.98927    final          30        0.5 1009.98927     15
#> 21          1 1352.31787    final          30        0.5 1352.31787     15
#> 22          1 1360.47827    final          30        0.5 1360.47827     15
#> 23          1 1442.95507    final          30        0.5 1442.95507     15
#> 24          1 1499.41809    final          30        0.5 1499.41809     15
#> 25          1 1513.98128    final          30        0.5 1513.98128     15
#> 26          1 1541.56362    final          30        0.5 1541.56362     15
#> 27          1 1547.92830    final          30        0.5 1547.92830     15
#> 28          1 1620.21011    final          30        0.5 1620.21011     15
#> 29          1 1665.91114    final          30        0.5 1665.91114     15
#> 30          1 1682.65379    final          30        0.5 1682.65379     15
#> 31          1 2105.77059    final          30        0.5 2105.77059     15
#> 32          2   13.65704    final          30        0.5   13.65704     15
#> 33          2   40.12907    final          30        0.5   40.12907     15
#> 34          2   48.87861    final          30        0.5   48.87861     15
#> 35          2   53.08020    final          30        0.5   53.08020     15
#> 36          2   68.06779    final          30        0.5   68.06779     15
#> 37          2  120.00897    final          30        0.5  120.00897     15
#> 38          2  274.27965    final          30        0.5  274.27965     15
#> 39          2  329.60709    final          30        0.5  329.60709     15
#> 40          2  334.34969    final          30        0.5  334.34969     15
#> 41          2  340.29801    final          30        0.5  340.29801     15
#> 42          2  344.60111    final          30        0.5  344.60111     15
#> 43          2  355.47256    final          30        0.5  355.47256     15
#> 44          2  401.54080    final          30        0.5  401.54080     15
#> 45          2  495.03695    final          30        0.5  495.03695     15
#> 46          2  661.02123    final          30        0.5  661.02123     15
#> 47          2  735.84724    final          30        0.5  735.84724     15
#> 48          2  915.28868    final          30        0.5  915.28868     15
#> 49          2  918.28933    final          30        0.5  918.28933     15
#> 50          2  987.22238    final          30        0.5  987.22238     15
#> 51          2 1095.30391    final          30        0.5 1095.30391     15
#> 52          2 1101.44770    final          30        0.5 1101.44770     15
#> 53          2 1211.99730    final          30        0.5 1211.99730     15
#> 54          2 1315.86879    final          30        0.5 1315.86879     15
#> 55          2 1320.83083    final          30        0.5 1320.83083     15
#> 56          2 1360.67220    final          30        0.5 1360.67220     15
#> 57          2 1401.86020    final          30        0.5 1401.86020     15
#> 58          2 1404.93461    final          30        0.5 1404.93461     15
#> 59          2 1596.86213    final          30        0.5 1596.86213     15
#> 60          2 1626.73816    final          30        0.5 1626.73816     15
#> 61          2 1749.41822    final          30        0.5 1749.41822     15
#> 62          2 1761.31338    final          30        0.5 1761.31338     15
#> 63          3   17.64924    final          30        0.5   17.64924     15
#> 64          3   41.82855    final          30        0.5   41.82855     15
#> 65          3  174.06645    final          30        0.5  174.06645     15
#> 66          3  177.04987    final          30        0.5  177.04987     15
#> 67          3  195.23773    final          30        0.5  195.23773     15
#> 68          3  281.56056    final          30        0.5  281.56056     15
#> 69          3  300.66545    final          30        0.5  300.66545     15
#> 70          3  415.11680    final          30        0.5  415.11680     15
#> 71          3  611.04363    final          30        0.5  611.04363     15
#> 72          3  612.36759    final          30        0.5  612.36759     15
#> 73          3  625.44020    final          30        0.5  625.44020     15
#> 74          3  693.85189    final          30        0.5  693.85189     15
#> 75          3  697.59443    final          30        0.5  697.59443     15
#> 76          3  735.24512    final          30        0.5  735.24512     15
#> 77          3  856.63921    final          30        0.5  856.63921     15
#> 78          3  918.47801    final          30        0.5  918.47801     15
#> 79          3 1027.50759    final          30        0.5 1027.50759     15
#> 80          3 1107.95781    final          30        0.5 1107.95781     15
#> 81          3 1144.45971    final          30        0.5 1144.45971     15
#> 82          3 1172.24833    final          30        0.5 1172.24833     15
#> 83          3 1208.98681    final          30        0.5 1208.98681     15
#> 84          3 1271.76920    final          30        0.5 1271.76920     15
#> 85          3 1276.06048    final          30        0.5 1276.06048     15
#> 86          3 1387.92622    final          30        0.5 1387.92622     15
#> 87          3 1434.36382    final          30        0.5 1434.36382     15
#> 88          3 1435.77955    final          30        0.5 1435.77955     15
#> 89          3 1437.27185    final          30        0.5 1437.27185     15
#> 90          3 1441.82329    final          30        0.5 1441.82329     15
#> 91          3 1479.38788    final          30        0.5 1479.38788     15
#> 92          3 1655.30887    final          30        0.5 1655.30887     15
#> 93          4   11.80808    final          30        0.5   11.80808     15
#> 94          4   74.32182    final          30        0.5   74.32182     15
#> 95          4  181.51216    final          30        0.5  181.51216     15
#> 96          4  221.58161    final          30        0.5  221.58161     15
#> 97          4  248.36705    final          30        0.5  248.36705     15
#> 98          4  349.95121    final          30        0.5  349.95121     15
#> 99          4  515.68996    final          30        0.5  515.68996     15
#> 100         4  666.42489    final          30        0.5  666.42489     15
#> 101         4  707.99707    final          30        0.5  707.99707     15
#> 102         4  762.41640    final          30        0.5  762.41640     15
#> 103         4  786.16173    final          30        0.5  786.16173     15
#> 104         4  821.15181    final          30        0.5  821.15181     15
#> 105         4  822.89551    final          30        0.5  822.89551     15
#> 106         4  848.83998    final          30        0.5  848.83998     15
#> 107         4  862.44843    final          30        0.5  862.44843     15
#> 108         4  963.33942    final          30        0.5  963.33942     15
#> 109         4  963.35346    final          30        0.5  963.35346     15
#> 110         4 1141.97410    final          30        0.5 1141.97410     15
#> 111         4 1335.03435    final          30        0.5 1335.03435     15
#> 112         4 1426.71867    final          30        0.5 1426.71867     15
#> 113         4 1529.34757    final          30        0.5 1529.34757     15
#> 114         4 1598.19506    final          30        0.5 1598.19506     15
#> 115         4 1643.82361    final          30        0.5 1643.82361     15
#> 116         4 1677.19237    final          30        0.5 1677.19237     15
#> 117         4 1739.77701    final          30        0.5 1739.77701     15
#> 118         4 1776.23710    final          30        0.5 1776.23710     15
#> 119         4 1820.69432    final          30        0.5 1820.69432     15
#> 120         4 1848.33883    final          30        0.5 1848.33883     15
#> 121         4 1893.76674    final          30        0.5 1893.76674     15
#> 122         4 1895.32186    final          30        0.5 1895.32186     15
#> 123         4 1935.66056    final          30        0.5 1935.66056     15
#> 124         5   12.91792    final          30        0.5   12.91792     15
#> 125         5   59.27306    final          30        0.5   59.27306     15
#> 126         5  106.74839    final          30        0.5  106.74839     15
#> 127         5  138.26688    final          30        0.5  138.26688     15
#> 128         5  196.58786    final          30        0.5  196.58786     15
#> 129         5  246.46181    final          30        0.5  246.46181     15
#> 130         5  311.90055    final          30        0.5  311.90055     15
#> 131         5  343.59568    final          30        0.5  343.59568     15
#> 132         5  463.86672    final          30        0.5  463.86672     15
#> 133         5  510.27762    final          30        0.5  510.27762     15
#> 134         5  523.92621    final          30        0.5  523.92621     15
#> 135         5  543.81839    final          30        0.5  543.81839     15
#> 136         5  547.34530    final          30        0.5  547.34530     15
#> 137         5  561.56543    final          30        0.5  561.56543     15
#> 138         5  730.24406    final          30        0.5  730.24406     15
#> 139         5  758.07099    final          30        0.5  758.07099     15
#> 140         5  766.55931    final          30        0.5  766.55931     15
#> 141         5  890.64011    final          30        0.5  890.64011     15
#> 142         5  953.51917    final          30        0.5  953.51917     15
#> 143         5 1045.85846    final          30        0.5 1045.85846     15
#> 144         5 1118.62309    final          30        0.5 1118.62309     15
#> 145         5 1154.77303    final          30        0.5 1154.77303     15
#> 146         5 1160.09774    final          30        0.5 1160.09774     15
#> 147         5 1238.64017    final          30        0.5 1238.64017     15
#> 148         5 1271.37483    final          30        0.5 1271.37483     15
#> 149         5 1291.99589    final          30        0.5 1291.99589     15
#> 150         5 1367.62195    final          30        0.5 1367.62195     15
#> 151         5 1583.05490    final          30        0.5 1583.05490     15
#> 152         5 1601.70445    final          30        0.5 1601.70445     15
#> 153         5 1619.82627    final          30        0.5 1619.82627     15
#> 154         5 1633.98107    final          30        0.5 1633.98107     15
#> 155         6   15.97475    final          30        0.5   15.97475     15
#> 156         6   30.58640    final          30        0.5   30.58640     15
#> 157         6   71.09775    final          30        0.5   71.09775     15
#> 158         6  118.61027    final          30        0.5  118.61027     15
#> 159         6  179.73142    final          30        0.5  179.73142     15
#> 160         6  182.16761    final          30        0.5  182.16761     15
#> 161         6  202.91827    final          30        0.5  202.91827     15
#> 162         6  222.21796    final          30        0.5  222.21796     15
#> 163         6  353.10786    final          30        0.5  353.10786     15
#> 164         6  500.78247    final          30        0.5  500.78247     15
#> 165         6  565.62785    final          30        0.5  565.62785     15
#> 166         6  797.41025    final          30        0.5  797.41025     15
#> 167         6  804.19935    final          30        0.5  804.19935     15
#> 168         6  812.67301    final          30        0.5  812.67301     15
#> 169         6  816.84071    final          30        0.5  816.84071     15
#> 170         6  932.16346    final          30        0.5  932.16346     15
#> 171         6  958.58778    final          30        0.5  958.58778     15
#> 172         6  984.39527    final          30        0.5  984.39527     15
#> 173         6  989.02383    final          30        0.5  989.02383     15
#> 174         6  993.37137    final          30        0.5  993.37137     15
#> 175         6  998.20217    final          30        0.5  998.20217     15
#> 176         6 1080.36898    final          30        0.5 1080.36898     15
#> 177         6 1131.86547    final          30        0.5 1131.86547     15
#> 178         6 1158.67270    final          30        0.5 1158.67270     15
#> 179         6 1286.41357    final          30        0.5 1286.41357     15
#> 180         6 1461.67934    final          30        0.5 1461.67934     15
#> 181         6 1513.03861    final          30        0.5 1513.03861     15
#> 182         6 1599.55679    final          30        0.5 1599.55679     15
#> 183         6 1711.45037    final          30        0.5 1711.45037     15
#> 184         6 1765.13555    final          30        0.5 1765.13555     15
#> 185         6 1790.76805    final          30        0.5 1790.76805     15
#> 186         7   14.76681    final          30        0.5   14.76681     15
#> 187         7   17.22987    final          30        0.5   17.22987     15
#> 188         7  191.74989    final          30        0.5  191.74989     15
#> 189         7  193.53281    final          30        0.5  193.53281     15
#> 190         7  250.97947    final          30        0.5  250.97947     15
#> 191         7  270.52755    final          30        0.5  270.52755     15
#> 192         7  278.99352    final          30        0.5  278.99352     15
#> 193         7  291.00901    final          30        0.5  291.00901     15
#> 194         7  297.45871    final          30        0.5  297.45871     15
#> 195         7  452.33066    final          30        0.5  452.33066     15
#> 196         7  495.31727    final          30        0.5  495.31727     15
#> 197         7  513.13144    final          30        0.5  513.13144     15
#> 198         7  514.71117    final          30        0.5  514.71117     15
#> 199         7  648.96664    final          30        0.5  648.96664     15
#> 200         7  715.55492    final          30        0.5  715.55492     15
#> 201         7  728.16967    final          30        0.5  728.16967     15
#> 202         7  930.86318    final          30        0.5  930.86318     15
#> 203         7  949.35259    final          30        0.5  949.35259     15
#> 204         7 1110.39403    final          30        0.5 1110.39403     15
#> 205         7 1111.53232    final          30        0.5 1111.53232     15
#> 206         7 1114.14650    final          30        0.5 1114.14650     15
#> 207         7 1152.02824    final          30        0.5 1152.02824     15
#> 208         7 1152.52184    final          30        0.5 1152.52184     15
#> 209         7 1155.77137    final          30        0.5 1155.77137     15
#> 210         7 1279.23964    final          30        0.5 1279.23964     15
#> 211         7 1281.29241    final          30        0.5 1281.29241     15
#> 212         7 1298.56498    final          30        0.5 1298.56498     15
#> 213         7 1311.29725    final          30        0.5 1311.29725     15
#> 214         7 1423.55366    final          30        0.5 1423.55366     15
#> 215         7 1480.44513    final          30        0.5 1480.44513     15
#> 216         8   13.66482    final          30        0.5   13.66482     15
#> 217         8   13.79376    final          30        0.5   13.79376     15
#> 218         8   51.92968    final          30        0.5   51.92968     15
#> 219         8  117.08581    final          30        0.5  117.08581     15
#> 220         8  144.34627    final          30        0.5  144.34627     15
#> 221         8  153.54673    final          30        0.5  153.54673     15
#> 222         8  168.91284    final          30        0.5  168.91284     15
#> 223         8  198.76434    final          30        0.5  198.76434     15
#> 224         8  292.53964    final          30        0.5  292.53964     15
#> 225         8  318.01040    final          30        0.5  318.01040     15
#> 226         8  341.60244    final          30        0.5  341.60244     15
#> 227         8  366.78101    final          30        0.5  366.78101     15
#> 228         8  417.96127    final          30        0.5  417.96127     15
#> 229         8  558.03142    final          30        0.5  558.03142     15
#> 230         8  572.38772    final          30        0.5  572.38772     15
#> 231         8  616.96598    final          30        0.5  616.96598     15
#> 232         8  680.30956    final          30        0.5  680.30956     15
#> 233         8  703.84945    final          30        0.5  703.84945     15
#> 234         8  725.56311    final          30        0.5  725.56311     15
#> 235         8  809.04486    final          30        0.5  809.04486     15
#> 236         8  829.65483    final          30        0.5  829.65483     15
#> 237         8  884.58893    final          30        0.5  884.58893     15
#> 238         8  920.04397    final          30        0.5  920.04397     15
#> 239         8  941.45513    final          30        0.5  941.45513     15
#> 240         8 1026.13573    final          30        0.5 1026.13573     15
#> 241         8 1043.48614    final          30        0.5 1043.48614     15
#> 242         8 1135.01653    final          30        0.5 1135.01653     15
#> 243         8 1158.19105    final          30        0.5 1158.19105     15
#> 244         8 1208.16100    final          30        0.5 1208.16100     15
#> 245         8 1263.30328    final          30        0.5 1263.30328     15
#> 246         8 1292.20260    final          30        0.5 1292.20260     15
#> 247         9   18.81521    final          30        0.5   18.81521     15
#> 248         9   45.25489    final          30        0.5   45.25489     15
#> 249         9  149.20340    final          30        0.5  149.20340     15
#> 250         9  151.30597    final          30        0.5  151.30597     15
#> 251         9  372.61226    final          30        0.5  372.61226     15
#> 252         9  713.06154    final          30        0.5  713.06154     15
#> 253         9  719.69989    final          30        0.5  719.69989     15
#> 254         9  756.09751    final          30        0.5  756.09751     15
#> 255         9  771.37198    final          30        0.5  771.37198     15
#> 256         9  820.63185    final          30        0.5  820.63185     15
#> 257         9  903.32705    final          30        0.5  903.32705     15
#> 258         9  985.54524    final          30        0.5  985.54524     15
#> 259         9 1067.51736    final          30        0.5 1067.51736     15
#> 260         9 1132.92803    final          30        0.5 1132.92803     15
#> 261         9 1170.10010    final          30        0.5 1170.10010     15
#> 262         9 1224.67732    final          30        0.5 1224.67732     15
#> 263         9 1333.22805    final          30        0.5 1333.22805     15
#> 264         9 1369.75168    final          30        0.5 1369.75168     15
#> 265         9 1432.93113    final          30        0.5 1432.93113     15
#> 266         9 1467.08644    final          30        0.5 1467.08644     15
#> 267         9 1531.31958    final          30        0.5 1531.31958     15
#> 268         9 1592.86738    final          30        0.5 1592.86738     15
#> 269         9 1637.76986    final          30        0.5 1637.76986     15
#> 270         9 1729.62034    final          30        0.5 1729.62034     15
#> 271         9 1756.87582    final          30        0.5 1756.87582     15
#> 272         9 1767.72488    final          30        0.5 1767.72488     15
#> 273         9 1779.95328    final          30        0.5 1779.95328     15
#> 274         9 1797.51059    final          30        0.5 1797.51059     15
#> 275         9 1835.94733    final          30        0.5 1835.94733     15
#> 276         9 1840.15397    final          30        0.5 1840.15397     15
#> 277         9 1853.03260    final          30        0.5 1853.03260     15
#> 278        10   11.84402    final          30        0.5   11.84402     15
#> 279        10   65.93802    final          30        0.5   65.93802     15
#> 280        10  114.79415    final          30        0.5  114.79415     15
#> 281        10  118.56286    final          30        0.5  118.56286     15
#> 282        10  249.50260    final          30        0.5  249.50260     15
#> 283        10  315.79657    final          30        0.5  315.79657     15
#> 284        10  356.10059    final          30        0.5  356.10059     15
#> 285        10  360.26315    final          30        0.5  360.26315     15
#> 286        10  364.76216    final          30        0.5  364.76216     15
#> 287        10  382.11500    final          30        0.5  382.11500     15
#> 288        10  441.63897    final          30        0.5  441.63897     15
#> 289        10  519.80005    final          30        0.5  519.80005     15
#> 290        10  522.19961    final          30        0.5  522.19961     15
#> 291        10  528.61700    final          30        0.5  528.61700     15
#> 292        10  606.97257    final          30        0.5  606.97257     15
#> 293        10  639.38183    final          30        0.5  639.38183     15
#> 294        10  689.94437    final          30        0.5  689.94437     15
#> 295        10  731.17913    final          30        0.5  731.17913     15
#> 296        10  766.21511    final          30        0.5  766.21511     15
#> 297        10  856.71461    final          30        0.5  856.71461     15
#> 298        10  979.28740    final          30        0.5  979.28740     15
#> 299        10 1169.52332    final          30        0.5 1169.52332     15
#> 300        10 1249.13102    final          30        0.5 1249.13102     15
#> 301        10 1254.76477    final          30        0.5 1254.76477     15
#> 302        10 1280.21756    final          30        0.5 1280.21756     15
#> 303        10 1296.35175    final          30        0.5 1296.35175     15
#> 304        10 1322.93170    final          30        0.5 1322.93170     15
#> 305        10 1399.38323    final          30        0.5 1399.38323     15
#> 306        10 1407.98763    final          30        0.5 1407.98763     15
#> 307        10 1428.71451    final          30        0.5 1428.71451     15
#> 308        10 1497.39038    final          30        0.5 1497.39038     15
#>     n_trt  mean_ctrl    mean_trt
#> 1      15 -0.2426125  0.73924372
#> 2      15 -0.2426125  0.73924372
#> 3      15 -0.2426125  0.73924372
#> 4      15 -0.2426125  0.73924372
#> 5      15 -0.2426125  0.73924372
#> 6      15 -0.2426125  0.73924372
#> 7      15 -0.2426125  0.73924372
#> 8      15 -0.2426125  0.73924372
#> 9      15 -0.2426125  0.73924372
#> 10     15 -0.2426125  0.73924372
#> 11     15 -0.2426125  0.73924372
#> 12     15 -0.2426125  0.73924372
#> 13     15 -0.2426125  0.73924372
#> 14     15 -0.2426125  0.73924372
#> 15     15 -0.2426125  0.73924372
#> 16     15 -0.2426125  0.73924372
#> 17     15 -0.2426125  0.73924372
#> 18     15 -0.2426125  0.73924372
#> 19     15 -0.2426125  0.73924372
#> 20     15 -0.2426125  0.73924372
#> 21     15 -0.2426125  0.73924372
#> 22     15 -0.2426125  0.73924372
#> 23     15 -0.2426125  0.73924372
#> 24     15 -0.2426125  0.73924372
#> 25     15 -0.2426125  0.73924372
#> 26     15 -0.2426125  0.73924372
#> 27     15 -0.2426125  0.73924372
#> 28     15 -0.2426125  0.73924372
#> 29     15 -0.2426125  0.73924372
#> 30     15 -0.2426125  0.73924372
#> 31     15 -0.2426125  0.73924372
#> 32     15  0.1855949  0.77785489
#> 33     15  0.1855949  0.77785489
#> 34     15  0.1855949  0.77785489
#> 35     15  0.1855949  0.77785489
#> 36     15  0.1855949  0.77785489
#> 37     15  0.1855949  0.77785489
#> 38     15  0.1855949  0.77785489
#> 39     15  0.1855949  0.77785489
#> 40     15  0.1855949  0.77785489
#> 41     15  0.1855949  0.77785489
#> 42     15  0.1855949  0.77785489
#> 43     15  0.1855949  0.77785489
#> 44     15  0.1855949  0.77785489
#> 45     15  0.1855949  0.77785489
#> 46     15  0.1855949  0.77785489
#> 47     15  0.1855949  0.77785489
#> 48     15  0.1855949  0.77785489
#> 49     15  0.1855949  0.77785489
#> 50     15  0.1855949  0.77785489
#> 51     15  0.1855949  0.77785489
#> 52     15  0.1855949  0.77785489
#> 53     15  0.1855949  0.77785489
#> 54     15  0.1855949  0.77785489
#> 55     15  0.1855949  0.77785489
#> 56     15  0.1855949  0.77785489
#> 57     15  0.1855949  0.77785489
#> 58     15  0.1855949  0.77785489
#> 59     15  0.1855949  0.77785489
#> 60     15  0.1855949  0.77785489
#> 61     15  0.1855949  0.77785489
#> 62     15  0.1855949  0.77785489
#> 63     15  0.2729714  0.52646154
#> 64     15  0.2729714  0.52646154
#> 65     15  0.2729714  0.52646154
#> 66     15  0.2729714  0.52646154
#> 67     15  0.2729714  0.52646154
#> 68     15  0.2729714  0.52646154
#> 69     15  0.2729714  0.52646154
#> 70     15  0.2729714  0.52646154
#> 71     15  0.2729714  0.52646154
#> 72     15  0.2729714  0.52646154
#> 73     15  0.2729714  0.52646154
#> 74     15  0.2729714  0.52646154
#> 75     15  0.2729714  0.52646154
#> 76     15  0.2729714  0.52646154
#> 77     15  0.2729714  0.52646154
#> 78     15  0.2729714  0.52646154
#> 79     15  0.2729714  0.52646154
#> 80     15  0.2729714  0.52646154
#> 81     15  0.2729714  0.52646154
#> 82     15  0.2729714  0.52646154
#> 83     15  0.2729714  0.52646154
#> 84     15  0.2729714  0.52646154
#> 85     15  0.2729714  0.52646154
#> 86     15  0.2729714  0.52646154
#> 87     15  0.2729714  0.52646154
#> 88     15  0.2729714  0.52646154
#> 89     15  0.2729714  0.52646154
#> 90     15  0.2729714  0.52646154
#> 91     15  0.2729714  0.52646154
#> 92     15  0.2729714  0.52646154
#> 93     15 -0.2590078  0.50676828
#> 94     15 -0.2590078  0.50676828
#> 95     15 -0.2590078  0.50676828
#> 96     15 -0.2590078  0.50676828
#> 97     15 -0.2590078  0.50676828
#> 98     15 -0.2590078  0.50676828
#> 99     15 -0.2590078  0.50676828
#> 100    15 -0.2590078  0.50676828
#> 101    15 -0.2590078  0.50676828
#> 102    15 -0.2590078  0.50676828
#> 103    15 -0.2590078  0.50676828
#> 104    15 -0.2590078  0.50676828
#> 105    15 -0.2590078  0.50676828
#> 106    15 -0.2590078  0.50676828
#> 107    15 -0.2590078  0.50676828
#> 108    15 -0.2590078  0.50676828
#> 109    15 -0.2590078  0.50676828
#> 110    15 -0.2590078  0.50676828
#> 111    15 -0.2590078  0.50676828
#> 112    15 -0.2590078  0.50676828
#> 113    15 -0.2590078  0.50676828
#> 114    15 -0.2590078  0.50676828
#> 115    15 -0.2590078  0.50676828
#> 116    15 -0.2590078  0.50676828
#> 117    15 -0.2590078  0.50676828
#> 118    15 -0.2590078  0.50676828
#> 119    15 -0.2590078  0.50676828
#> 120    15 -0.2590078  0.50676828
#> 121    15 -0.2590078  0.50676828
#> 122    15 -0.2590078  0.50676828
#> 123    15 -0.2590078  0.50676828
#> 124    15  0.1494857  0.90866496
#> 125    15  0.1494857  0.90866496
#> 126    15  0.1494857  0.90866496
#> 127    15  0.1494857  0.90866496
#> 128    15  0.1494857  0.90866496
#> 129    15  0.1494857  0.90866496
#> 130    15  0.1494857  0.90866496
#> 131    15  0.1494857  0.90866496
#> 132    15  0.1494857  0.90866496
#> 133    15  0.1494857  0.90866496
#> 134    15  0.1494857  0.90866496
#> 135    15  0.1494857  0.90866496
#> 136    15  0.1494857  0.90866496
#> 137    15  0.1494857  0.90866496
#> 138    15  0.1494857  0.90866496
#> 139    15  0.1494857  0.90866496
#> 140    15  0.1494857  0.90866496
#> 141    15  0.1494857  0.90866496
#> 142    15  0.1494857  0.90866496
#> 143    15  0.1494857  0.90866496
#> 144    15  0.1494857  0.90866496
#> 145    15  0.1494857  0.90866496
#> 146    15  0.1494857  0.90866496
#> 147    15  0.1494857  0.90866496
#> 148    15  0.1494857  0.90866496
#> 149    15  0.1494857  0.90866496
#> 150    15  0.1494857  0.90866496
#> 151    15  0.1494857  0.90866496
#> 152    15  0.1494857  0.90866496
#> 153    15  0.1494857  0.90866496
#> 154    15  0.1494857  0.90866496
#> 155    15  0.4779952 -0.03207007
#> 156    15  0.4779952 -0.03207007
#> 157    15  0.4779952 -0.03207007
#> 158    15  0.4779952 -0.03207007
#> 159    15  0.4779952 -0.03207007
#> 160    15  0.4779952 -0.03207007
#> 161    15  0.4779952 -0.03207007
#> 162    15  0.4779952 -0.03207007
#> 163    15  0.4779952 -0.03207007
#> 164    15  0.4779952 -0.03207007
#> 165    15  0.4779952 -0.03207007
#> 166    15  0.4779952 -0.03207007
#> 167    15  0.4779952 -0.03207007
#> 168    15  0.4779952 -0.03207007
#> 169    15  0.4779952 -0.03207007
#> 170    15  0.4779952 -0.03207007
#> 171    15  0.4779952 -0.03207007
#> 172    15  0.4779952 -0.03207007
#> 173    15  0.4779952 -0.03207007
#> 174    15  0.4779952 -0.03207007
#> 175    15  0.4779952 -0.03207007
#> 176    15  0.4779952 -0.03207007
#> 177    15  0.4779952 -0.03207007
#> 178    15  0.4779952 -0.03207007
#> 179    15  0.4779952 -0.03207007
#> 180    15  0.4779952 -0.03207007
#> 181    15  0.4779952 -0.03207007
#> 182    15  0.4779952 -0.03207007
#> 183    15  0.4779952 -0.03207007
#> 184    15  0.4779952 -0.03207007
#> 185    15  0.4779952 -0.03207007
#> 186    15  0.1935956  0.05437164
#> 187    15  0.1935956  0.05437164
#> 188    15  0.1935956  0.05437164
#> 189    15  0.1935956  0.05437164
#> 190    15  0.1935956  0.05437164
#> 191    15  0.1935956  0.05437164
#> 192    15  0.1935956  0.05437164
#> 193    15  0.1935956  0.05437164
#> 194    15  0.1935956  0.05437164
#> 195    15  0.1935956  0.05437164
#> 196    15  0.1935956  0.05437164
#> 197    15  0.1935956  0.05437164
#> 198    15  0.1935956  0.05437164
#> 199    15  0.1935956  0.05437164
#> 200    15  0.1935956  0.05437164
#> 201    15  0.1935956  0.05437164
#> 202    15  0.1935956  0.05437164
#> 203    15  0.1935956  0.05437164
#> 204    15  0.1935956  0.05437164
#> 205    15  0.1935956  0.05437164
#> 206    15  0.1935956  0.05437164
#> 207    15  0.1935956  0.05437164
#> 208    15  0.1935956  0.05437164
#> 209    15  0.1935956  0.05437164
#> 210    15  0.1935956  0.05437164
#> 211    15  0.1935956  0.05437164
#> 212    15  0.1935956  0.05437164
#> 213    15  0.1935956  0.05437164
#> 214    15  0.1935956  0.05437164
#> 215    15  0.1935956  0.05437164
#> 216    15  0.2760705  0.16748723
#> 217    15  0.2760705  0.16748723
#> 218    15  0.2760705  0.16748723
#> 219    15  0.2760705  0.16748723
#> 220    15  0.2760705  0.16748723
#> 221    15  0.2760705  0.16748723
#> 222    15  0.2760705  0.16748723
#> 223    15  0.2760705  0.16748723
#> 224    15  0.2760705  0.16748723
#> 225    15  0.2760705  0.16748723
#> 226    15  0.2760705  0.16748723
#> 227    15  0.2760705  0.16748723
#> 228    15  0.2760705  0.16748723
#> 229    15  0.2760705  0.16748723
#> 230    15  0.2760705  0.16748723
#> 231    15  0.2760705  0.16748723
#> 232    15  0.2760705  0.16748723
#> 233    15  0.2760705  0.16748723
#> 234    15  0.2760705  0.16748723
#> 235    15  0.2760705  0.16748723
#> 236    15  0.2760705  0.16748723
#> 237    15  0.2760705  0.16748723
#> 238    15  0.2760705  0.16748723
#> 239    15  0.2760705  0.16748723
#> 240    15  0.2760705  0.16748723
#> 241    15  0.2760705  0.16748723
#> 242    15  0.2760705  0.16748723
#> 243    15  0.2760705  0.16748723
#> 244    15  0.2760705  0.16748723
#> 245    15  0.2760705  0.16748723
#> 246    15  0.2760705  0.16748723
#> 247    15 -0.5170736  0.59219044
#> 248    15 -0.5170736  0.59219044
#> 249    15 -0.5170736  0.59219044
#> 250    15 -0.5170736  0.59219044
#> 251    15 -0.5170736  0.59219044
#> 252    15 -0.5170736  0.59219044
#> 253    15 -0.5170736  0.59219044
#> 254    15 -0.5170736  0.59219044
#> 255    15 -0.5170736  0.59219044
#> 256    15 -0.5170736  0.59219044
#> 257    15 -0.5170736  0.59219044
#> 258    15 -0.5170736  0.59219044
#> 259    15 -0.5170736  0.59219044
#> 260    15 -0.5170736  0.59219044
#> 261    15 -0.5170736  0.59219044
#> 262    15 -0.5170736  0.59219044
#> 263    15 -0.5170736  0.59219044
#> 264    15 -0.5170736  0.59219044
#> 265    15 -0.5170736  0.59219044
#> 266    15 -0.5170736  0.59219044
#> 267    15 -0.5170736  0.59219044
#> 268    15 -0.5170736  0.59219044
#> 269    15 -0.5170736  0.59219044
#> 270    15 -0.5170736  0.59219044
#> 271    15 -0.5170736  0.59219044
#> 272    15 -0.5170736  0.59219044
#> 273    15 -0.5170736  0.59219044
#> 274    15 -0.5170736  0.59219044
#> 275    15 -0.5170736  0.59219044
#> 276    15 -0.5170736  0.59219044
#> 277    15 -0.5170736  0.59219044
#> 278    15  0.4187062  0.42788315
#> 279    15  0.4187062  0.42788315
#> 280    15  0.4187062  0.42788315
#> 281    15  0.4187062  0.42788315
#> 282    15  0.4187062  0.42788315
#> 283    15  0.4187062  0.42788315
#> 284    15  0.4187062  0.42788315
#> 285    15  0.4187062  0.42788315
#> 286    15  0.4187062  0.42788315
#> 287    15  0.4187062  0.42788315
#> 288    15  0.4187062  0.42788315
#> 289    15  0.4187062  0.42788315
#> 290    15  0.4187062  0.42788315
#> 291    15  0.4187062  0.42788315
#> 292    15  0.4187062  0.42788315
#> 293    15  0.4187062  0.42788315
#> 294    15  0.4187062  0.42788315
#> 295    15  0.4187062  0.42788315
#> 296    15  0.4187062  0.42788315
#> 297    15  0.4187062  0.42788315
#> 298    15  0.4187062  0.42788315
#> 299    15  0.4187062  0.42788315
#> 300    15  0.4187062  0.42788315
#> 301    15  0.4187062  0.42788315
#> 302    15  0.4187062  0.42788315
#> 303    15  0.4187062  0.42788315
#> 304    15  0.4187062  0.42788315
#> 305    15  0.4187062  0.42788315
#> 306    15  0.4187062  0.42788315
#> 307    15  0.4187062  0.42788315
#> 308    15  0.4187062  0.42788315
```

### collect_results() and the analysis= filter

When a trial has multiple named analyses (e.g., an interim and a final),
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
returns rows for all of them. Use the `analysis` argument to restrict to
a specific name:

``` r
# Only the final analysis rows
collect_results(trials, analysis = "final")
```

Each row in the output corresponds to one replicate firing one named
analysis. The `replicate`, `timepoint`, and `analysis` columns identify
the provenance of every result row, making it straightforward to compare
interim and final results within the same replicate or aggregate across
replicates.

## Next steps

- [Getting
  Started](https://boehringer-ingelheim.github.io/rxsim/articles/getting-started.md)
  — end-to-end walkthrough of the six-step pattern shown above, with
  commentary on each design choice.
- [Enrollment & Dropout
  Modeling](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment-dropout.md)
  — choose between stochastic (`gen_plan`) and piecewise-constant
  (`gen_timepoints`) schedules.
- [Example
  1](https://boehringer-ingelheim.github.io/rxsim/articles/example-1.md)
  through [Example
  8](https://boehringer-ingelheim.github.io/rxsim/articles/example-8.md)
  — progressively complex designs: correlated endpoints, time-to-event,
  multi-arm dose-finding, subgroup analyses, and Bayesian Go/No-Go
  rules.

# Getting Started

This vignette walks through a complete rxsim simulation from start to
finish: defining a trial scenario, building arm populations, registering
an analysis trigger, running replicates, and collecting results. By the
end you will have seen the full rxsim workflow in one place and will be
ready to adapt it to your own design. For deeper explanations of each
building block, see the [Core
Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.md)
vignette.

## Load the package

``` r
library(rxsim)
```

## Step 1 — Define the scenario

The `scenario` tibble captures the design parameters you want to track
in your results — sample size, allocation ratios, or any other factors
you might vary across simulation runs. Wrapping them in
[`tidyr::expand_grid()`](https://tidyr.tidyverse.org/reference/expand_grid.html)
produces a one-row data frame that gets embedded in every analysis
result. When you later stack results across many scenarios or
sensitivity analyses, this metadata keeps each row traceable back to its
design.

``` r
sample_size <- 30
arms        <- c("control", "treatment")
allocation  <- c(1, 1)
true_delta  <- 0.5

scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation  = list(allocation),
  true_delta  = true_delta
)
```

## Step 2 — Define arm populations

Each trial arm is represented by a **generator function** — a plain R
function that takes `n` (the number of subjects planned for that arm)
and returns a data frame of subject-level endpoint data. rxsim calls
these functions once per replicate to draw a fresh population, so every
replicate gets independent data.

The data frame must contain at least three columns:

- `id` — a unique integer identifier per subject
- `readout_time` — the time (in trial-time units) after enrollment when
  the endpoint is observed (use `0` for baseline or
  immediately-available data)
- at least one endpoint column

[`vector_to_dataframe()`](https://boehringer-ingelheim.github.io/rxsim/reference/vector_to_dataframe.md)
is a convenience helper that wraps a numeric vector into this standard
format with `id`, `data`, and `readout_time = 0`.

Here we simulate a simple continuous endpoint: the control arm is
standard-normal, the treatment arm has a mean shift of `true_delta`.

``` r
population_generators <- list(
  control   = function(n) vector_to_dataframe(rnorm(n)),
  treatment = function(n) vector_to_dataframe(rnorm(n, mean = true_delta))
)
```

## Step 3 — Define enrollment and dropout

Enrollment and dropout are modelled as random processes. Each function
takes `n` and returns a vector of **inter-event times** — the waiting
times between successive enrollments (or dropouts). These times are
drawn independently for every replicate, giving each trial its own
enrollment trajectory.

`rexp(n, rate = 1)` generates exponentially-distributed inter-arrival
times with a mean gap of 1 time unit between enrolments — a common
approximation for Poisson arrivals. A lower dropout rate (`0.05`)
reflects a trial where most subjects complete the study.

``` r
enrollment <- function(n) rexp(n, rate = 1.0)
dropout    <- function(n) rexp(n, rate = 0.05)
```

## Step 4 — Define analysis triggers

Analysis triggers are the heart of rxsim. Each trigger pairs a
**condition** (when to fire) with an **analysis function** (what to
compute when it fires).

The condition is written as a `dplyr`-style boolean expression inside
[`rlang::exprs()`](https://rlang.r-lib.org/reference/defusing-advanced.html).
It is evaluated against a snapshot of all currently enrolled subjects at
each timepoint. Here the condition fires once full enrollment is
reached, i.e., `sample_size` subjects have accumulated a non-`NA`
`enroll_time`.

The `!!` operator (pronounced “bang-bang”) **injects the current value**
of `sample_size` into the expression at definition time, rather than
looking it up at evaluation time. This is necessary because the
expression is stored and evaluated later inside the simulation loop.

When the condition is met, rxsim calls the analysis function with two
arguments:

- `df`: a data frame snapshot of all enrolled subjects at the triggering
  timepoint, with columns from the population data plus `enroll_time`,
  `drop_time`, and `arm`
- `time`: the current trial clock time at which the trigger fired

The function should return a data frame (one row per trigger event is
the conventional pattern).

``` r
analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer) {
      enrolled <- subset(df, !is.na(enroll_time))
      data.frame(
        scenario,
        n_enrolled  = nrow(enrolled),
        mean_ctrl   = mean(enrolled$data[enrolled$arm == "control"]),
        mean_trt    = mean(enrolled$data[enrolled$arm == "treatment"]),
        stringsAsFactors = FALSE
      )
    }
  )
)
```

## Step 5 — Create replicates and run

[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
generates `n` fully independent `Trial` objects. For each replicate it:

1.  Samples a fresh enrollment/dropout timeline from your functions
2.  Calls each population generator to draw new subject-level data
3.  Registers your analysis triggers on the trial’s internal timer

[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md)
then executes every replicate’s simulation loop in sequence.

``` r
set.seed(42)

trials <- replicate_trial(
  trial_name            = "getting_started",
  sample_size           = sample_size,
  arms                  = arms,
  allocation            = allocation,
  enrollment            = enrollment,
  dropout               = dropout,
  analysis_generators   = analysis_generators,
  population_generators = population_generators,
  n                     = 5
)

run_trials(trials)
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: getting_started_1
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
#>     name: getting_started_2
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
#>     name: getting_started_3
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
#>     name: getting_started_4
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
#>     name: getting_started_5
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
```

## Step 6 — Collect and interpret results

### Analysis results

After running, each `Trial` object exposes a `results` list. It is
indexed first by the **timepoint** at which an analysis fired (e.g.,
`"time_30.4"`), then by the **analysis name** you gave it (here
`"final"`). The value is whatever your analysis function returned.

The
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
helper gathers every analysis from every timepoint across all replicates
into one tidy data frame:

``` r
replicate_results <- collect_results(trials)
replicate_results
#>   replicate timepoint analysis sample_size allocation true_delta n_enrolled
#> 1         1  29.71727    final          30       1, 1        0.5         30
#> 2         2  27.31408    final          30       1, 1        0.5         30
#> 3         3  35.29848    final          30       1, 1        0.5         30
#> 4         4  23.61616    final          30       1, 1        0.5         30
#> 5         5  25.83584    final          30       1, 1        0.5         30
#>    mean_ctrl   mean_trt
#> 1 -0.1657292 0.61748974
#> 2 -0.0978443 0.64303730
#> 3 -0.2573593 0.07050019
#> 4 -0.4199802 0.25120329
#> 5 -0.4170091 0.59879107
```

Each row is one replicate. The `mean_ctrl` and `mean_trt` columns show
the arm-level sample means at the moment the final analysis fired. The
`timepoint` and `analysis` columns identify when and which analysis
produced each row — essential when a trial has both interim and final
analyses. Variation across replicates reflects both the stochastic
endpoint data and the randomness in enrollment timing; exactly what
operating characteristic simulations are designed to characterise.

### Locked data snapshots

Every time an analysis fires, rxsim also saves a **locked data
snapshot**, the full subject-level data frame at that timepoint. You can
inspect it directly to debug your analysis function, audit which
subjects were enrolled, or compute additional statistics post-hoc.

``` r
# One snapshot per timepoint that fired in replicate 1
names(trials[[1]]$locked_data)
#> [1] "time_29.7172711929598"

# First six rows of the snapshot
head(trials[[1]]$locked_data[[1]])
#>   id       data readout_time     arm enroll_time drop_time subject_id
#> 1  1 -2.6882473            0 control   0.8592321        NA          1
#> 2  2  0.8666502            0 control  11.6939492        NA          2
#> 3  3  0.1687397            0 control  26.6106187        NA          3
#> 4  4 -1.0908241            0 control   1.6540916        NA          4
#> 5  5 -0.3803481            0 control   9.5016932        NA          5
#> 6  6 -0.9480918            0 control  13.9319792        NA          6
#>   measurement_time     time
#> 1        0.8592321 29.71727
#> 2       11.6939492 29.71727
#> 3       26.6106187 29.71727
#> 4        1.6540916 29.71727
#> 5        9.5016932 29.71727
#> 6       13.9319792 29.71727
```

The locked data contains the population columns (`id`, `data`, `arm`,
`readout_time`) plus three columns added by rxsim:

| Column        | Meaning                                                           |
|---------------|-------------------------------------------------------------------|
| `enroll_time` | Calendar time the subject was enrolled (`NA` if not yet enrolled) |
| `drop_time`   | Calendar time the subject dropped out (`NA` if still active)      |
| `time`        | The trial clock time at which this snapshot was taken             |

## Next steps

Now that you have seen the full workflow, here is where to go next:

- **[Core
  Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.md)**
  — understand how `Population`, `Timer`, and `Trial` compose, and how
  to write more advanced trigger expressions
- **[Enrollment & Dropout
  Modeling](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment-dropout.md)**
  — choose between stochastic (`gen_plan`) and piecewise-constant
  (`gen_timepoints`) schedules
- **[Example
  1](https://boehringer-ingelheim.github.io/rxsim/articles/example-1.md)**
  through **[Example
  7](https://boehringer-ingelheim.github.io/rxsim/articles/example-7.md)**
  — progressively complex designs: correlated endpoints, time-to-event,
  multi-arm dose-finding, subgroup analyses, and Bayesian Go/No-Go rules

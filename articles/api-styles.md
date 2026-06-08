# Two API Styles: Direct and Generator

``` r
library(rxsim)
set.seed(42)
```

## Introduction

rxsim supports two equivalent ways to build the same trial.

- **Style A: Direct instantiation** builds `Population`, `Timer`,
  `Condition`, and `Trial` objects yourself.
- **Style B: Generator API** gives
  [`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
  the population and analysis generators, then lets rxsim build each
  `Trial` for you.

The two styles are equivalent when the design itself is the same: same
arms, endpoint distributions, enrollment, dropout, and trigger logic.
The difference is where you want control.

Prefer the direct style when you want to inspect or edit the `Timer`,
work with fixed schedules, or debug one trial step by step. Prefer the
generator style when you want many stochastic replicates, quick scenario
grids, or minimal boilerplate.

The sections below use the same two-arm design in both styles:

- placebo (`pbo`) versus treatment (`trt`)
- `sample_size = 40` with `allocation = c(1, 1)`
- stochastic enrollment from `rexp(n, rate = 1)`
- stochastic dropout from `rexp(n, rate = 0.05)`
- continuous endpoint with mean `0` for placebo and `0.5` for treatment
- one final analysis at full enrollment
- `100` replicates for the operating characteristics

``` r
sample_size <- 40L
arms <- c("pbo", "trt")
allocation <- c(1, 1)
delta <- 0.5
n_reps <- 100L

enrollment <- function(n) rexp(n, rate = 1)
dropout <- function(n) rexp(n, rate = 0.05)

make_outcome_data <- function(n, mean_shift) {
  data.frame(
    id = seq_len(n),
    outcome = rnorm(n, mean = mean_shift, sd = 1),
    readout_time = 0,
    stringsAsFactors = FALSE
  )
}

final_analysis <- function(df, current_time) {
  enrolled <- subset(df, !is.na(enroll_time))
  fit <- stats::t.test(outcome ~ arm, data = enrolled)

  data.frame(
    n = nrow(enrolled),
    mean_pbo = mean(enrolled$outcome[enrolled$arm == "pbo"]),
    mean_trt = mean(enrolled$outcome[enrolled$arm == "trt"]),
    p_value = unname(fit$p.value),
    stringsAsFactors = FALSE
  )
}

make_final_condition <- function() {
  Condition$new(
    where = enroll_trigger(1.0, sample_size),
    analysis = final_analysis,
    name = "final"
  )
}
```

## Style A: Direct instantiation

Use the direct style when you want to see every moving part. You create
the plan, build the timer, size each arm, define the trigger, and
assemble the trial yourself.

``` r
# 1. Draw one stochastic enrollment and dropout plan.
direct_plan <- stochastic_schedule(
  sample_size = sample_size,
  arms = arms,
  allocation = allocation,
  enrollment = enrollment,
  dropout = dropout
)

# 2. Build the timer from that plan.
direct_timer <- Timer$new(name = "direct_timer")
add_timepoints(direct_timer, direct_plan)

# 3. Compute the planned arm sizes from the timer input.
direct_n_by_arm <- vapply(
  arms,
  function(arm_name) {
    as.integer(sum(direct_plan$enroll[direct_plan$arm == arm_name]))
  },
  integer(1)
)

# 4. Instantiate each arm population explicitly.
pop_pbo <- Population$new(
  name = "pbo",
  data = make_outcome_data(direct_n_by_arm[["pbo"]], mean_shift = 0)
)
pop_trt <- Population$new(
  name = "trt",
  data = make_outcome_data(direct_n_by_arm[["trt"]], mean_shift = delta)
)

# 5. Define the final analysis with a condition object.
direct_final <- make_final_condition()

# 6. Assemble the Trial object.
direct_trial <- Trial$new(
  name = "direct_trial",
  timer = direct_timer,
  population = list(pop_pbo, pop_trt),
  conditions = list(direct_final)
)
```

In this style, the timer is a first-class object. You can inspect the
time grid, modify individual timepoints, or replace the whole schedule
before running.

``` r
head(direct_plan)
#>        time arm enroll drop
#> 1 0.1983368 trt      1    0
#> 2 0.8592321 pbo      1    0
#> 3 1.1427231 trt      1    0
#> 4 1.1809150 pbo      1    0
#> 5 1.6540916 pbo      1    0
#> 6 3.1177188 trt      1    0
direct_timer$get_end_timepoint()
#> [1] 1110.694
```

``` r
# 7. Run the single direct trial.
direct_trial$run()
```

``` r
collect_results(direct_trial)
#>   replicate timepoint analysis  n   mean_pbo  mean_trt    p_value
#> 1         1  39.32605    final 40 -0.2365031 0.3138118 0.07738052
```

For many stochastic replicates, the direct pattern is usually wrapped in
a small constructor function.

``` r
build_direct_trial <- function(name) {
  # 1. Generate a fresh stochastic plan for this replicate.
  plan <- stochastic_schedule(
    sample_size = sample_size,
    arms = arms,
    allocation = allocation,
    enrollment = enrollment,
    dropout = dropout
  )

  # 2. Turn that plan into a timer.
  timer <- Timer$new(name = paste0(name, "_timer"))
  add_timepoints(timer, plan)

  # 3. Size each arm from the generated plan.
  n_by_arm <- vapply(
    arms,
    function(arm_name) {
      as.integer(sum(plan$enroll[plan$arm == arm_name], na.rm = TRUE))
    },
    integer(1)
  )

  # 4. Instantiate populations explicitly.
  populations <- list(
    Population$new(
      name = "pbo",
      data = make_outcome_data(n_by_arm[["pbo"]], mean_shift = 0)
    ),
    Population$new(
      name = "trt",
      data = make_outcome_data(n_by_arm[["trt"]], mean_shift = delta)
    )
  )

  # 5. Create a fresh condition object for this replicate.
  conditions <- list(make_final_condition())

  # 6. Return the assembled trial.
  Trial$new(
    name = name,
    timer = timer,
    population = populations,
    conditions = conditions
  )
}
```

``` r
direct_trials <- lapply(seq_len(n_reps), function(i) {
  build_direct_trial(paste0("direct_", i))
})

run_trials(direct_trials)
```

``` r
direct_results <- collect_results(direct_trials)
head(direct_results)
#>   replicate timepoint analysis  n    mean_pbo  mean_trt     p_value
#> 1         1  41.37433    final 40 -0.15072151 0.4098671 0.054209202
#> 2         2  46.89639    final 40 -0.24803270 0.6183037 0.006305818
#> 3         3  39.48357    final 40  0.09571239 0.6604110 0.027212451
#> 4         4  58.20288    final 40  0.27907034 0.5495985 0.538070676
#> 5         5  41.54251    final 40  0.32055414 0.1959181 0.687516418
#> 6         6  37.94218    final 40  0.31120700 0.4925153 0.491698998
```

## Style B: Generator API

Use the generator API when you want the same design with less setup. You
define functions for the arm populations and analyses, then let
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
build the timers, conditions, and trials internally.

``` r
# 1. Define one population generator per arm.
population_generators <- list(
  pbo = function(n) make_outcome_data(n, mean_shift = 0),
  trt = function(n) make_outcome_data(n, mean_shift = delta)
)

# 2. Define the analysis generator with a condition object.
analysis_generators <- list(
  final = list(
    trigger = enroll_trigger(1.0, sample_size),
    analysis = final_analysis
  )
)
```

The trigger must be a `trigger` object such as
`enroll_trigger(1.0, sample_size)`. No quotation helpers are needed.

``` r
# 3. Ask rxsim to build 100 independent trials.
generator_trials <- replicate_trial(
  trial_name = "generator",
  sample_size = sample_size,
  arms = arms,
  allocation = allocation,
  enrollment = enrollment,
  dropout = dropout,
  analysis_generators = analysis_generators,
  population_generators = population_generators,
  n = n_reps
)

# 4. Run all generated trials.
run_trials(generator_trials)
```

``` r
generator_results <- collect_results(generator_trials)
head(generator_results)
#>   replicate timepoint analysis  n    mean_pbo  mean_trt    p_value
#> 1         1  42.90242    final 40  0.17727332 0.3057526 0.76280511
#> 2         2  33.13166    final 40 -0.26273620 0.3530079 0.09358722
#> 3         3  41.63035    final 40  0.04993081 0.2030695 0.68585888
#> 4         4  55.59473    final 40 -0.33050838 0.4749362 0.01359331
#> 5         5  30.49985    final 40  0.18127567 0.6535889 0.18802339
#> 6         6  36.88387    final 40  0.01316165 0.5411758 0.06682232
```

This style is shorter because
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
handles the timer construction, population sizing, and condition
creation for every replicate.

## When to use Style A

Prefer direct instantiation when you want to work directly with the
underlying objects.

- Manually manipulate the `Timer` before running.
- Use deterministic schedules from
  [`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
  or a hand-built plan.
- Debug a single trial by inspecting `timer`, `population`, and
  `locked_data`.
- Build multi-stage designs where the schedule is easier to express
  directly, such as the seamless pattern in [Example
  5](https://boehringer-ingelheim.github.io/rxsim/articles/example-5.md).

## When to use Style B

Prefer the generator API when the design follows the standard workflow
that
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
already knows how to build.

- Run scenario grids generated with `expand_grid()`.
- Simulate standard stochastic enrollment and dropout with minimal
  setup.
- Write scripts that need many replicates and little boilerplate.
- Keep design changes local to a few small generator functions.

## Mixing the styles

A common workflow is to start with the direct style for one trial,
inspect it, then scale out with
[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md)
once you are happy with the setup.

``` r
# 1. Build one trial directly so you can inspect it.
mixed_trial <- build_direct_trial("mixed_1")
head(do.call(rbind, mixed_trial$timer$timelist))

# 2. Extend the same direct pattern to more replicates.
mixed_trials <- c(
  list(mixed_trial),
  lapply(2:5, function(i) build_direct_trial(paste0("mixed_", i)))
)

# 3. Run them together with the usual batch helper.
run_trials(mixed_trials)
```

If you intentionally want the same fixed timer and populations in each
replicate,
[`clone_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/clone_trial.md)
is also useful. That pattern is most natural for deterministic
schedules, not for stochastic enrollment.

## Next steps

- Add interims and multiple analyses: [Conditions and
  Triggers](https://boehringer-ingelheim.github.io/rxsim/articles/conditions.md)
- Learn more about trigger objects and analyses: [Conditions and
  Triggers](https://boehringer-ingelheim.github.io/rxsim/articles/conditions.md)
- Inspect the `Trial` class documentation: [Trial
  reference](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)

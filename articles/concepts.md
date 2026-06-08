# Core Concepts

## Overview

rxsim organises a clinical trial simulation around four collaborating
objects. A `Population` owns the subject-level data and tracks each
subject’s enrollment and dropout times. A `Timer` drives the trial
clock: it stores discrete timepoints per arm and defines when the
simulation clock advances. A `Condition` pairs a filter expression with
an optional analysis function and manages its own trigger state - it
fires when the snapshot data meets a criterion. A `Trial` orchestrates
the simulation by iterating over timepoints, updating populations,
snapshotting the enrolled cohort, and collecting results.

``` mermaid
graph LR
  P1(Control) --> TR(Trial)
  P2(Treatment) --> TR
  TI(Timer) --> TR
  CO(Condition) --> TR
  TR --> LD(Locked Data)
  TR --> RS(Results)
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
  [`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)
  cannot express
- use `Trial$new()` directly for a one-off single-run simulation (as in
  [Example
  5](https://boehringer-ingelheim.github.io/rxsim/articles/example-5.md))

The sections below give a short summary of each building block. For the
full reference, see the dedicated vignettes linked in the table.

## Building blocks at a glance

| Class        | Role                                                                        | Deep dive                                                                                      |
|--------------|-----------------------------------------------------------------------------|------------------------------------------------------------------------------------------------|
| `Population` | Holds subject-level endpoint data and enrollment/dropout state for one arm  | [Population](https://boehringer-ingelheim.github.io/rxsim/articles/class-population.md)        |
| `Timer`      | Stores the trial clock: when subjects enroll or drop in each arm            | [Enrollment and Dropout](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment.md)  |
| `Condition`  | Pairs a trigger expression with an analysis function; manages trigger state | [Conditions and Triggers](https://boehringer-ingelheim.github.io/rxsim/articles/conditions.md) |
| `Trial`      | Orchestrates the simulation loop; stores snapshots and results              | [Trial reference](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)             |

## How the pieces fit together

A typical two-arm simulation follows this pattern:

``` r
# 1. Populations: one per arm
pop_pbo <- Population$new("pbo", pbo_data)
pop_trt <- Population$new("trt", trt_data)

# 2. Timer: register enrollment/dropout events
tmr <- Timer$new("my_timer")
add_timepoints(tmr, stochastic_schedule(n, arms, alloc, enroll_fn, dropout_fn))

# 3. Condition: fire at full enrollment, run a t-test
cond <- Condition$new(
  where    = enroll_trigger(1.0, n),
  analysis = function(df, t) { ... },
  name     = "final"
)

# 4. Trial: assemble and run
trial <- Trial$new("trial", timer=tmr,
                   population=list(pop_pbo, pop_trt),
                   conditions=list(cond))
trial$run()
```

For a step-by-step walkthrough see [Enrollment and
Dropout](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment.md).
For the generator shortcut (`replicate_trial`) see [Two API
Styles](https://boehringer-ingelheim.github.io/rxsim/articles/api-styles.md).

## Next steps

- [Enrollment and
  Dropout](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment.md) -
  stochastic and deterministic enrollment
- [Population](https://boehringer-ingelheim.github.io/rxsim/articles/class-population.md) -
  subject data format and enrollment tracking
- [Conditions and
  Triggers](https://boehringer-ingelheim.github.io/rxsim/articles/conditions.md) -
  triggers and analysis functions
- [Trial
  reference](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md) -
  class documentation for run loop and outputs
- [Two API
  Styles](https://boehringer-ingelheim.github.io/rxsim/articles/api-styles.md) -
  direct vs. generator API

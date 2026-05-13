# Package index

## Core classes

R6 classes that form the simulation building blocks

- [`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
  : Condition: Stateful trigger and analysis unit
- [`Timer`](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md)
  : Timer: Track timed events across arms
- [`Trial`](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)
  : Trial: Simulate a multi‑arm clinical trial
- [`Population`](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md)
  : Population: Manage a patient population

## Trial construction helpers

Functions to build and replicate trials

- [`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
  : Create Multiple Trials with Generated Populations
- [`clone_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/clone_trial.md)
  : Clone a Trial Object Multiple Times
- [`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md)
  : Run Multiple Trial Objects

## Timer helpers

Helpers for building Timer timepoint lists

- [`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md)
  : Add Timepoints to a Timer
- [`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
  : Generate Piecewise-Linear Enrollment and Dropout Plan

## Condition helpers

Convenience constructors for Condition objects

- [`value_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`count_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`enroll_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`calendar_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`` `&`( ``*`<rxsim_trigger>`*`)`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`` `|`( ``*`<rxsim_trigger>`*`)`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  : Build Safe Trial Triggers
- [`trigger_by_calendar()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_calendar.md)
  : Trigger Analysis at a Calendar Time
- [`trigger_by_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_fraction.md)
  : Trigger Analysis at a Sample Fraction

## Population helpers

Helpers for building Population generators

- [`gen_population()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_population.md)
  : Generate a Population Object
- [`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
  : Generate Trial Enrollment and Dropout Plan

## Results helpers

Functions to collect and display simulation results

- [`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
  : Collect Trial Results Across Replicates
- [`prettify_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/prettify_results.md)
  : Format Trial Results as a Data Frame

## Utilities

- [`get_col_names()`](https://boehringer-ingelheim.github.io/rxsim/reference/get_col_names.md)
  : Extract Column Names from Populations
- [`vector_to_dataframe()`](https://boehringer-ingelheim.github.io/rxsim/reference/vector_to_dataframe.md)
  : Convert Vector to Population Data Frame

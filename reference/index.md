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
- [`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
  : Generate a Deterministic Enrollment and Dropout Schedule
- [`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)
  : Generate a Stochastic Enrollment and Dropout Schedule

## Condition helpers

Convenience constructors for Condition objects

- [`value_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`count_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`enroll_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`calendar_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`` `&`( ``*`<trigger>`*`)`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  [`` `|`( ``*`<trigger>`*`)`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  : Build Trial Triggers
- [`condition_calendar_time()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_calendar_time.md)
  : Build a Condition that Fires at a Calendar Time
- [`condition_enrollment_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_enrollment_fraction.md)
  : Build a Condition that Fires at an Enrollment Fraction

## Population helpers

Helpers for building Population generators

- [`gen_population()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_population.md)
  : Generate a Population Object
- [`as_population_data()`](https://boehringer-ingelheim.github.io/rxsim/reference/as_population_data.md)
  : Create a Population-Compatible Data Frame from a Vector

## Results helpers

Functions to collect and display simulation results

- [`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
  : Collect Trial Results Across Replicates

## Utilities

- [`get_col_names()`](https://boehringer-ingelheim.github.io/rxsim/reference/get_col_names.md)
  : Extract Column Names from Populations

# Create Multiple Trials with Generated Populations

Generates `n` independent `Trial` objects from a template, timer, and
population generators. Arms must match between plan and population
specifications.

## Usage

``` r
replicate_trial(
  trial_name = "name",
  sample_size,
  arms,
  allocation,
  enrollment,
  dropout,
  analysis_generators,
  population_generators,
  n
)
```

## Arguments

- trial_name:

  `character` Base name for generated trials (index appended). Defaults
  to `"name"`.

- sample_size:

  `integer` Total sample size across all arms.

- arms:

  `character` vector of arm identifiers.

- allocation:

  `numeric` vector of arm allocation ratios.

- enrollment:

  `function` that generates inter-enrollment times.

- dropout:

  `function` that generates inter-dropout times.

- analysis_generators:

  `list` (named) of analysis trigger specifications.

- population_generators:

  `list` (named) of population generator functions.

- n:

  `integer` Number of trials to create.

## Value

`list` of `n` `Trial` objects with indexed names, cloned timers, and
generated populations.

## See also

[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md),
[`gen_population()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_population.md),
[`clone_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/clone_trial.md),
[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md).

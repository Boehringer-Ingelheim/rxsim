# Clone a Trial Object Multiple Times

Creates deep copies of a `Trial` R6 object with independent timer and
population instances.

## Usage

``` r
clone_trial(trial, n = 1)
```

## Arguments

- trial:

  `Trial` R6 object to clone.

- n:

  `integer` Number of clones to create. Defaults to `1`.

## Value

`list` of `n` independently cloned `Trial` objects.

## See also

[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md),
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md),
[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md).

## Examples

``` r
# Create 5 independent copies of a trial:
# clones <- clone_trial(trial, n = 5)
```

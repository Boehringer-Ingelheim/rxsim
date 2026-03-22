# Run Multiple Trial Objects

Executes the `$run()` method for each trial in a list.

## Usage

``` r
run_trials(trials)
```

## Arguments

- trials:

  `list` of `Trial` R6 objects.

## Value

`list` of results from each trial's `$run()` method.

## See also

[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md),
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md).

## Examples

``` r
# results <- run_trials(trials)
```

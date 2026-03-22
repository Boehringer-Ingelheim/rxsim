# Generate a Population Object

Creates a `Population` R6 object by executing a generator function.

## Usage

``` r
gen_population(name, generator, sample_size = 1)
```

## Arguments

- name:

  `character` Population name (e.g., arm label).

- generator:

  `function` that returns a data.frame-like object.

- sample_size:

  `integer` Number of subjects to generate.

## Value

`Population` R6 object.

## See also

[Population](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md),
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md),
[`vector_to_dataframe()`](https://boehringer-ingelheim.github.io/rxsim/reference/vector_to_dataframe.md).

## Examples

``` r
gen_control <- function(n) data.frame(id = 1:n, value = rnorm(n), readout_time = 0)
pop <- gen_population("control", gen_control, sample_size = 50)
```

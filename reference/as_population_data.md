# Create a Population-Compatible Data Frame from a Vector

Converts a numeric vector to a data frame with the columns required by
[`Population`](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md):
`id`, `data`, and `readout_time`.

## Usage

``` r
as_population_data(data)
```

## Arguments

- data:

  `numeric` vector of endpoint values (one per subject).

## Value

`data.frame` with columns: `id` (integer), `data` (numeric),
`readout_time` (0 for all subjects).

## See also

[Population](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md).

## Examples

``` r
as_population_data(rnorm(5))
#>   id       data readout_time
#> 1  1  0.5038124            0
#> 2  2  2.5283366            0
#> 3  3  0.5490967            0
#> 4  4  0.2382129            0
#> 5  5 -1.0488931            0
```

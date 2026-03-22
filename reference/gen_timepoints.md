# Generate Piecewise-Linear Enrollment and Dropout Plan

Creates a time-indexed plan with piecewise constant enrollment and
dropout rates.

## Usage

``` r
gen_timepoints(sample_size, arms, allocation, enrollment, dropout)
```

## Arguments

- sample_size:

  `integer` Trial sample size.

- arms:

  `character` vector of arm identifiers.

- allocation:

  `numeric` vector of allocation ratios.

- enrollment:

  `list` with `end_time` (numeric endpoints) and `rate` (numeric
  subjects/unit time for each period).

- dropout:

  `list` with `end_time` and `rate` (same structure).

## Value

`data.frame` with columns: `time`, `arm`, `enroller`, `dropper`.

## See also

[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
for random inter-event times,
[`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md).

## Examples

``` r
gen_timepoints(
  sample_size = 100,
  arms = c("A", "B"),
  allocation = c(2, 1),
  enrollment = list(
    end_time = c(4, 8, 12),
    rate = c(6, 12, 18)
  ),
  dropout = list(
    end_time = c(5, 9, 13),
    rate = c(0, 3, 6)
  )
)
#> # A tibble: 20 × 4
#>     time arm   enroller dropper
#>    <dbl> <chr>    <int>   <int>
#>  1     1 A            4       0
#>  2     2 A            4       0
#>  3     3 A            4       0
#>  4     4 A            4       0
#>  5     5 A            8       0
#>  6     6 A            8       2
#>  7     7 A            8       2
#>  8     8 A            8       2
#>  9     9 A           12       2
#> 10    10 A            7       2
#> 11     1 B            2       0
#> 12     2 B            2       0
#> 13     3 B            2       0
#> 14     4 B            2       0
#> 15     5 B            4       0
#> 16     6 B            4       1
#> 17     7 B            4       1
#> 18     8 B            4       1
#> 19     9 B            6       1
#> 20    10 B            3       1
```

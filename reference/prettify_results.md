# Format Trial Results as a Data Frame

Converts trial results to a single data frame with all measurements.

## Usage

``` r
prettify_results(results)
```

## Arguments

- results:

  `list` Trial results (nested by time).

## Value

`data.frame` with columns: `time` and measurement columns.

## See also

[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)
for generating results.

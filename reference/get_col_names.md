# Extract Column Names from Populations

Collects all data frame column names from one or more populations.

## Usage

``` r
get_col_names(populations)
```

## Arguments

- populations:

  [`Population`](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md)
  object or `list` of
  [`Population`](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md)
  objects.

## Value

`character` vector of unique column names.

## See also

[Population](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md).

## Examples

``` r
pop1 <- Population$new(name = "P1", data = data.frame(
id = 1:10,
age = runif(10, 20, 60),
readout_time = 0
))
pop2 <- Population$new(name = "P2", data = data.frame(
id = 1:10,
weight = runif(10, 150, 250),
readout_time = 0
))
get_col_names(list(pop1, pop2))
#> [1] "id"           "age"          "readout_time" "arm"          "weight"      
#> [6] "time"         "enroll_time"  "drop_time"    "measure_time"
```

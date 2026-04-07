# Example 2: Two arm \| Interim \| Continuous \| t-test

``` r
library(rxsim)
```

This example extends [Example
1](https://boehringer-ingelheim.github.io/rxsim/articles/example-1.md)
by introducing two named analyses — `interim` and `final` — that fire at
different enrollment milestones. Understanding how analysis names
propagate through
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
is the foundation for multi-analysis designs shown in later examples.

## Scenario

Capture scenario parameters. We will assume piece-wise linear
enrollment.

``` r
sample_size <- 100
arms        <- c("pbo", "trt")
allocation  <- c(1, 1)
delta       <- 0.2  # true treatment - placebo mean difference
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn    <- function(n) rexp(n, rate = 0.01)
scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation  = list(allocation),
  delta       = delta
)
```

`allocation = c(1, 1)` specifies balanced randomisation.
[`tidyr::expand_grid()`](https://tidyr.tidyverse.org/reference/expand_grid.html)
embeds design parameters directly into each result row for traceability
across parameter sweeps.

## Populations

Define population generators.

``` r
population_generators <- list(
  pbo = function(n) data.frame(
    id = 1:n,
    value = rnorm(n),
    readout_time = 1
  ),
  trt = function(n) data.frame(
    id = 1:n,
    value = rnorm(n, mean = delta),
    readout_time = 1
  )
)
```

Each generator is a function of `n` returning a `data.frame` with one
row per subject. `readout_time = 1` means the endpoint is observed 1
time unit after enrollment. The treatment arm has a mean shift of δ =
0.2 SD units over placebo.

## Triggers & Analysis

Final analysis at full enrollment.

``` r
analysis_generators <- list(
  interim = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= (!!sample_size) / 2
    ),
    analysis = function(df, timer){
      df_enrolled <- df |> subset(!is.na(enroll_time))
      tt <- t.test(value ~ arm, data = df_enrolled)
      data.frame(
        scenario,
        n_total = nrow(df_enrolled),
        mean_pbo = mean(df_enrolled$value[df_enrolled$arm == "pbo"]),
        mean_trt = mean(df_enrolled$value[df_enrolled$arm == "trt"]),
        p_value = unname(tt$p.value),
        stringsAsFactors = FALSE
      )
    }
  ),
  
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer){
      df_enrolled <- df |> subset(!is.na(enroll_time))
      tt <- t.test(value ~ arm, data = df_enrolled)
      data.frame(
        scenario,
        n_total = nrow(df_enrolled),
        mean_pbo = mean(df_enrolled$value[df_enrolled$arm == "pbo"]),
        mean_trt = mean(df_enrolled$value[df_enrolled$arm == "trt"]),
        p_value = unname(tt$p.value),
        stringsAsFactors = FALSE
      )
    }
  )
)
```

Two analyses are defined. The `interim` analysis fires at half
enrollment and the `final` analysis fires when all subjects are
enrolled. These names appear in the `analysis` column of
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
output, making it straightforward to distinguish interim from final
results when stacking rows across timepoints.

## Trial

Make multiple trial replicates.

``` r
trials <- replicate_trial(
  trial_name = "test_trial",
  sample_size = sample_size,
  arms = arms,
  allocation = allocation,
  enrollment = enrollment_fn,
  dropout = dropout_fn,
  analysis_generators = analysis_generators,
  population_generators = population_generators,
  n = 3
)
```

## Simulate

To simulate all replicates:

``` r
run_trials(trials)
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: test_trial_1
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[2]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: test_trial_2
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[3]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: test_trial_3
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
```

Bind one row per replicate into one data frame.

``` r
replicate_results <- collect_results(trials)
replicate_results
#>   replicate timepoint analysis sample_size allocation delta n_total    mean_pbo
#> 1         1  48.42750  interim         100       1, 1   0.2      50 -0.15458732
#> 2         1  87.36847    final         100       1, 1   0.2     100 -0.02475451
#> 3         2  40.26049  interim         100       1, 1   0.2      50 -0.04934773
#> 4         2  83.92580    final         100       1, 1   0.2     100 -0.05831826
#> 5         3  43.51109  interim         100       1, 1   0.2      50 -0.14213989
#> 6         3  96.44952    final         100       1, 1   0.2     100 -0.25739017
#>      mean_trt    p_value
#> 1  0.42624025 0.01626231
#> 2  0.32833062 0.05877739
#> 3  0.02698932 0.80070480
#> 4  0.13052493 0.36397787
#> 5 -0.13708169 0.98578393
#> 6  0.08945477 0.09253380
```

[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
prepends `replicate`, `timepoint`, and `analysis` columns to your result
data. Each row shows either `"interim"` or `"final"` in the `analysis`
column, reflecting which named analysis produced it. This naming
convention becomes essential in later examples where multiple named
analyses fire per replicate.

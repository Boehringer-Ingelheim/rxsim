# Example 2: Two arm \| Interim \| Continuous \| t-test

``` r
library(rxsim)
```

## Scenario

Capture scenario parameters. We will assume piece-wise linear
enrollment.

``` r
sample_size <- 100
arms <- c("pbo", "trt")
allocation <- c(1,1)
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn <- function(n) rexp(n, rate = 0.01)
scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation = list(allocation)
)
```

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
    value = rnorm(n, 0.2),
    readout_time = 1
  )
)
```

## Triggers & Analysis

Final analysis at full enrollment.

``` r
analysis_generators <- list(
  interim_final = list(
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
replicate_results <- do.call(rbind, lapply(seq_along(trials), function(i) {
  data.frame(
    replicate = i,
    trials[[i]]$results[[1]]$interim_final,
    stringsAsFactors = FALSE
  )
}))

replicate_results
#>   replicate sample_size allocation n_total    mean_pbo   mean_trt    p_value
#> 1         1         100       1, 1     100 -0.02475451 0.32833062 0.05877739
#> 2         2         100       1, 1     100 -0.05831826 0.13052493 0.36397787
#> 3         3         100       1, 1     100 -0.25739017 0.08945477 0.09253380
```

# Example Multi Run: Two arm \| Interim \| Continuous \| t-test

``` r
library(rxsim)
```

## Scenario

Capture scenario parameters. We will assume piece-wise linear
enrollment.

``` r
set.seed(123)
sample_size <- 6
arms        <- c("pbo", "trt")
allocation  <- c(1, 1)
enrollment  <- function(n) rexp(n, rate = 1)
dropout     <- function(n) rexp(n, rate = 0.01)
scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation = list(allocation)
)
```

## Arms

Generate populations and endpoint data using generators (functions).

``` r
#update the generators to have sample size 50

population_generators <- list(
  pbo = function(n) vector_to_dataframe(stats::rnorm(n)),
  trt = function(n) vector_to_dataframe(stats::rnorm(n, mean = 2, sd = 3.14))
)
```

## Triggers & Analysis

Run one final analysis at full enrollment.

``` r
analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer) {
      df_enrolled <- df |> subset(!is.na(enroll_time))
      mean_pbo <- df_enrolled |> subset(arm == "pbo") |> with(mean(data))
      mean_trt <- df_enrolled |> subset(arm == "trt") |> with(mean(data))
      data.frame(
        scenario,
        mean_pbo = mean_pbo,
        mean_trt = mean_trt,
        n_total = nrow(df_enrolled),
        stringsAsFactors = FALSE
      )
    }
  )
)
```

## Trial

Make the trial objects for when data is using generators.

``` r
trials <- replicate_trial(
  trial_name    = "trial",
  sample_size   = sample_size,
  arms          = arms,
  allocation    = allocation,
  enrollment    = enrollment,
  dropout       = dropout,
  analysis_generators = analysis_generators,
  population_generators = population_generators,
  n             = 3
)
```

## Simulate

Run all replicates and bind one row per replicate.

``` r
run_trials(trials)
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: trial_1
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
#>     name: trial_2
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
#>     name: trial_3
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6

replicate_results <- do.call(rbind, lapply(seq_along(trials), function(i) {
  data.frame(
    replicate = i,
    trials[[i]]$results[[1]]$final,
    stringsAsFactors = FALSE
  )
}))

replicate_results
#>   replicate sample_size allocation   mean_pbo   mean_trt n_total
#> 1         1           6       1, 1  0.3666536 2.57810702       6
#> 2         2           6       1, 1  0.2398118 0.05936882       6
#> 3         3           6       1, 1 -0.1066843 5.16653968       6
```

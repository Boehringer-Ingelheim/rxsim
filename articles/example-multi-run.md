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
delta       <- 2    # true treatment - placebo mean difference
enrollment  <- function(n) rexp(n, rate = 1)
dropout     <- function(n) rexp(n, rate = 0.01)
scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation  = list(allocation),
  delta       = delta
)
```

## Arms

Generate populations and endpoint data using generators (functions).

``` r
population_generators <- list(
  pbo = function(n) vector_to_dataframe(stats::rnorm(n)),
  trt = function(n) vector_to_dataframe(stats::rnorm(n, mean = delta, sd = 3.14))
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
replicate_results <- collect_results(trials)
replicate_results
#>    replicate  timepoint analysis sample_size allocation delta   mean_pbo
#> 1          1   3.153412    final           6       1, 1     2  0.3666536
#> 2          1  31.422729    final           6       1, 1     2  0.3666536
#> 3          1  45.949410    final           6       1, 1     2  0.3666536
#> 4          1 318.573056    final           6       1, 1     2  0.3666536
#> 5          1 321.488401    final           6       1, 1     2  0.3666536
#> 6          1 421.971407    final           6       1, 1     2  0.3666536
#> 7          1 469.992879    final           6       1, 1     2  0.3666536
#> 8          2   9.851882    final           6       1, 1     2  0.2398118
#> 9          2 160.585234    final           6       1, 1     2  0.2398118
#> 10         2 310.259521    final           6       1, 1     2  0.2398118
#> 11         2 467.324776    final           6       1, 1     2  0.2398118
#> 12         2 470.501550    final           6       1, 1     2  0.2398118
#> 13         2 530.286519    final           6       1, 1     2  0.2398118
#> 14         2 747.070494    final           6       1, 1     2  0.2398118
#> 15         3  10.932265    final           6       1, 1     2 -0.1066843
#> 16         3 224.830569    final           6       1, 1     2 -0.1066843
#> 17         3 361.203999    final           6       1, 1     2 -0.1066843
#> 18         3 418.843166    final           6       1, 1     2 -0.1066843
#> 19         3 691.370751    final           6       1, 1     2 -0.1066843
#> 20         3 822.587055    final           6       1, 1     2 -0.1066843
#> 21         3 831.646190    final           6       1, 1     2 -0.1066843
#>      mean_trt n_total
#> 1  2.57810702       6
#> 2  2.57810702       6
#> 3  2.57810702       6
#> 4  2.57810702       6
#> 5  2.57810702       6
#> 6  2.57810702       6
#> 7  2.57810702       6
#> 8  0.05936882       6
#> 9  0.05936882       6
#> 10 0.05936882       6
#> 11 0.05936882       6
#> 12 0.05936882       6
#> 13 0.05936882       6
#> 14 0.05936882       6
#> 15 5.16653968       6
#> 16 5.16653968       6
#> 17 5.16653968       6
#> 18 5.16653968       6
#> 19 5.16653968       6
#> 20 5.16653968       6
#> 21 5.16653968       6
```

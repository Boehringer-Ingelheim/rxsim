# Getting Started

## Setup package

``` r
library(rxsim)
```

## Scenario

Capture essential trial parameters. These will form scenarios.

``` r
sample_size <- 15
allocation <- c(1)
scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation = list(allocation)
)
```

## Arm / Population / Cohort

Let’s assume a continuous endpoint that follows a normal distribution.

``` r
population_generators <- list(
  control = function(n) vector_to_dataframe(rnorm(n))
)
```

## Trial Parameters

Define enrollment and dropout profiles.

``` r
arms <- c("control")
enrollment <- function(n) rexp(n, rate = 1)
dropout <- function(n) rexp(n, rate = 0.1)
```

## Event Triggers

Define analysis conditions.

``` r
analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer) {
      data.frame(
        scenario,
        n_enrolled = sum(!is.na(df$enroll_time)),
        mean_data = mean(df$data),
        stringsAsFactors = FALSE
      )
    }
  )
)
```

## Collect & Run

Create multiple replicates of the trial:

``` r
trials <- replicate_trial(
  trial_name = "test_trial",
  sample_size = sample_size,
  arms = arms,
  allocation = allocation,
  enrollment = enrollment,
  dropout = dropout,
  analysis_generators = analysis_generators,
  population_generators = population_generators,
  n = 3
)
```

Run all replicates:

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

replicate_results <- do.call(rbind, lapply(seq_along(trials), function(i) {
  data.frame(
    replicate = i,
    trials[[i]]$results[[1]]$final,
    stringsAsFactors = FALSE
  )
}))

replicate_results
#>   replicate sample_size allocation n_enrolled   mean_data
#> 1         1          15          1         15  0.24991547
#> 2         2          15          1         15 -0.12779334
#> 3         3          15          1         15 -0.03437337
```

You can also inspect locked data from any replicate:

``` r
trials[[1]]$locked_data |> lapply(head)
#> $time_20.8297407415279
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811        NA          2
#> 3  3 -0.16267634            0 control   10.060903        NA          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 20.82974
#> 2        14.224811 20.82974
#> 3        10.060903 20.82974
#> 4        16.176464 20.82974
#> 5         3.040110 20.82974
#> 6         4.885593 20.82974
#> 
#> $time_21.5488625294756
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811        NA          2
#> 3  3 -0.16267634            0 control   10.060903        NA          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 21.54886
#> 2        14.224811 21.54886
#> 3        10.060903 21.54886
#> 4        16.176464 21.54886
#> 5         3.040110 21.54886
#> 6         4.885593 21.54886
#> 
#> $time_21.9054609428141
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811        NA          2
#> 3  3 -0.16267634            0 control   10.060903        NA          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 21.90546
#> 2        14.224811 21.90546
#> 3        10.060903 21.90546
#> 4        16.176464 21.90546
#> 5         3.040110 21.90546
#> 6         4.885593 21.90546
#> 
#> $time_25.8219385224375
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903        NA          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 25.82194
#> 2        14.224811 25.82194
#> 3        10.060903 25.82194
#> 4        16.176464 25.82194
#> 5         3.040110 25.82194
#> 6         4.885593 25.82194
#> 
#> $time_29.5930585901234
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903        NA          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 29.59306
#> 2        14.224811 29.59306
#> 3        10.060903 29.59306
#> 4        16.176464 29.59306
#> 5         3.040110 29.59306
#> 6         4.885593 29.59306
#> 
#> $time_65.8138979666
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903        NA          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time    time
#> 1        17.483763 65.8139
#> 2        14.224811 65.8139
#> 3        10.060903 65.8139
#> 4        16.176464 65.8139
#> 5         3.040110 65.8139
#> 6         4.885593 65.8139
#> 
#> $time_69.9425747219538
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903        NA          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 69.94257
#> 2        14.224811 69.94257
#> 3        10.060903 69.94257
#> 4        16.176464 69.94257
#> 5         3.040110 69.94257
#> 6         4.885593 69.94257
#> 
#> $time_71.1925878081882
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903  71.19259          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 71.19259
#> 2        14.224811 71.19259
#> 3        10.060903 71.19259
#> 4        16.176464 71.19259
#> 5         3.040110 71.19259
#> 6         4.885593 71.19259
#> 
#> $time_74.7958461599836
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903  71.19259          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 74.79585
#> 2        14.224811 74.79585
#> 3        10.060903 74.79585
#> 4        16.176464 74.79585
#> 5         3.040110 74.79585
#> 6         4.885593 74.79585
#> 
#> $time_82.3929067106495
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763        NA          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903  71.19259          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 82.39291
#> 2        14.224811 82.39291
#> 3        10.060903 82.39291
#> 4        16.176464 82.39291
#> 5         3.040110 82.39291
#> 6         4.885593 82.39291
#> 
#> $time_87.7864027965392
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763  87.78640          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903  71.19259          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time    time
#> 1        17.483763 87.7864
#> 2        14.224811 87.7864
#> 3        10.060903 87.7864
#> 4        16.176464 87.7864
#> 5         3.040110 87.7864
#> 6         4.885593 87.7864
#> 
#> $time_90.4840121289538
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763  87.78640          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903  71.19259          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593        NA          6
#>   measurement_time     time
#> 1        17.483763 90.48401
#> 2        14.224811 90.48401
#> 3        10.060903 90.48401
#> 4        16.176464 90.48401
#> 5         3.040110 90.48401
#> 6         4.885593 90.48401
#> 
#> $time_100.025154031562
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763  87.78640          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903  71.19259          3
#> 4  4 -0.82731017            0 control   16.176464        NA          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593 100.02515          6
#>   measurement_time     time
#> 1        17.483763 100.0252
#> 2        14.224811 100.0252
#> 3        10.060903 100.0252
#> 4        16.176464 100.0252
#> 5         3.040110 100.0252
#> 6         4.885593 100.0252
#> 
#> $time_117.862169160347
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763  87.78640          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903  71.19259          3
#> 4  4 -0.82731017            0 control   16.176464 117.86217          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593 100.02515          6
#>   measurement_time     time
#> 1        17.483763 117.8622
#> 2        14.224811 117.8622
#> 3        10.060903 117.8622
#> 4        16.176464 117.8622
#> 5         3.040110 117.8622
#> 6         4.885593 117.8622
#> 
#> $time_118.493639975626
#>   id        data readout_time     arm enroll_time drop_time subject_id
#> 1  1  0.89496163            0 control   17.483763  87.78640          1
#> 2  2  0.06730444            0 control   14.224811  25.82194          2
#> 3  3 -0.16267634            0 control   10.060903  71.19259          3
#> 4  4 -0.82731017            0 control   16.176464 117.86217          4
#> 5  5  1.87650562            0 control    3.040110        NA          5
#> 6  6  0.76644020            0 control    4.885593 100.02515          6
#>   measurement_time     time
#> 1        17.483763 118.4936
#> 2        14.224811 118.4936
#> 3        10.060903 118.4936
#> 4        16.176464 118.4936
#> 5         3.040110 118.4936
#> 6         4.885593 118.4936
```

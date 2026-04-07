# Example 1: Two arm \| Fixed design \| Continuous

``` r
library(rxsim)
```

This example simulates a Phase II/III parallel-group trial evaluating a
new treatment against placebo on a continuous primary endpoint (e.g., a
biomarker score or symptom scale). With 100 subjects randomised 1:1, we
run a two-sample t-test at full enrollment and track the p-value and arm
means across simulation replicates. This is the simplest complete rxsim
workflow and a good starting point before moving to more complex
designs. See [Getting
Started](https://boehringer-ingelheim.github.io/rxsim/articles/getting-started.md)
for an overview of the package.

## Scenario

Capture scenario parameters. We will assume piece-wise linear
enrollment.

``` r
sample_size <- 100
arms <- c("pbo", "trt")
allocation <- c(1,1)
delta <- 0.2
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn <- function(n) rexp(n, rate = 0.01)
scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation = list(allocation),
  delta = delta
)
```

`allocation = c(1, 1)` specifies balanced randomisation, equal expected
numbers per arm. Using
[`tidyr::expand_grid()`](https://tidyr.tidyverse.org/reference/expand_grid.html)
to build the `scenario` data.frame embeds the design parameters directly
into each analysis result row, making results self-documenting and easy
to compare across parameter sweeps.

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
    value = rnorm(n, delta),
    readout_time = 1
  )
)
```

Each population generator is a function of `n` that returns a
`data.frame` with one row per subject. `readout_time = 1` means the
endpoint is observed exactly 1 time unit after a subject enrolls. The
placebo arm has a mean of 0 and the treatment arm a mean of `delta`, a
small-to-moderate standardised effect.

## Triggers & Analysis

We want to do a t-test when `sample_size` subjects have been enrolled.

``` r
analysis_generators <- list(
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

The trigger fires when the cumulative count of enrolled subjects
(`sum(!is.na(enroll_time))`) reaches `sample_size`. Inside the analysis
function, `subset(!is.na(enroll_time))` drops subjects who have been
allocated but not yet enrolled. rxsim pre-generates the full allocation
list so filtering is essential. The two-sample t-test then compares
endpoint values between the two enrolled arms.

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
#> 1         1  87.36847    final         100       1, 1   0.2     100 -0.02475451
#> 2         2  83.92580    final         100       1, 1   0.2     100 -0.05831826
#> 3         3  96.44952    final         100       1, 1   0.2     100 -0.25739017
#>     mean_trt    p_value
#> 1 0.32833062 0.05877739
#> 2 0.13052493 0.36397787
#> 3 0.08945477 0.09253380
```

[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
returns a data.frame with `replicate`, `timepoint`, and `analysis`
columns prepended to your analysis columns. The `p_value` column varies
across replicates because each replicate draws fresh random data.

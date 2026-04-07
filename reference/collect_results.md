# Collect Trial Results Across Replicates

Gathers analysis outputs from one or more `Trial` objects into a single
tidy data frame. Works with any number of named analyses (e.g., both an
interim and a final) and any number of replicates.

## Usage

``` r
collect_results(trials, analysis = NULL)
```

## Arguments

- trials:

  A `Trial` R6 object **or** a `list` of `Trial` objects (as returned by
  [`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)).

- analysis:

  `character` or `NULL`. When supplied, only analyses whose name matches
  one of these values are included. Defaults to `NULL` (all analyses).

## Value

A `data.frame` with columns:

- `replicate` `integer` Index of the trial replicate (1-based).

- `timepoint` `numeric` Calendar time at which the analysis fired.

- `analysis` `character` Name of the analysis (as given in
  `analysis_generators` or via a trigger helper).

- Additional columns from the value returned by each analysis function.

## Details

Each analysis function may return either a `data.frame` (the standard
pattern) or a named `list`; both are coerced to a single-row data frame
and stacked. Analyses that return `NULL` or `NA` are silently skipped.

When `trials` is a single `Trial` object (e.g., from a one-off
`Trial$new()` + `$run()` call), the `replicate` column is always `1`.

## See also

[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md),
[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md),
[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md).

## Examples

``` r
# --- replicated trial ---
pop_gens <- list(
  control   = function(n) vector_to_dataframe(rnorm(n)),
  treatment = function(n) vector_to_dataframe(rnorm(n, 0.5))
)
an_gens <- list(
  final = list(
    trigger  = rlang::exprs(sum(!is.na(enroll_time)) >= 20L),
    analysis = function(df, timer) {
      data.frame(mean_ctrl = mean(df$data[df$arm == "control"]))
    }
  )
)
trials <- replicate_trial(
  trial_name = "ex", sample_size = 20L,
  arms = c("control", "treatment"), allocation = c(1, 1),
  enrollment = function(n) rexp(n, 1), dropout = function(n) rexp(n, 0.01),
  analysis_generators = an_gens, population_generators = pop_gens, n = 3
)
run_trials(trials)
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: ex_1
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
#>     name: ex_2
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
#>     name: ex_3
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
collect_results(trials)
#>   replicate timepoint analysis  mean_ctrl
#> 1         1  22.03729    final  0.1219045
#> 2         2  20.51980    final -0.2771686
#> 3         3  17.30206    final -0.2397427

# --- filter to a specific analysis name ---
collect_results(trials, analysis = "final")
#>   replicate timepoint analysis  mean_ctrl
#> 1         1  22.03729    final  0.1219045
#> 2         2  20.51980    final -0.2771686
#> 3         3  17.30206    final -0.2397427
```

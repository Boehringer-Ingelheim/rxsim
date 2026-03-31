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
#>    replicate  timepoint analysis  mean_ctrl
#> 1          1   22.03729    final  0.1219045
#> 2          1   50.66157    final  0.1219045
#> 3          1   76.61735    final  0.1219045
#> 4          1  336.30657    final  0.1219045
#> 5          1  459.20914    final  0.1219045
#> 6          1  538.27732    final  0.1219045
#> 7          1  601.20532    final  0.1219045
#> 8          1  726.66942    final  0.1219045
#> 9          1  785.53789    final  0.1219045
#> 10         1  898.46689    final  0.1219045
#> 11         1  940.50337    final  0.1219045
#> 12         1 1661.60413    final  0.1219045
#> 13         1 1746.17633    final  0.1219045
#> 14         1 1768.73053    final  0.1219045
#> 15         1 1878.76441    final  0.1219045
#> 16         1 2103.59498    final  0.1219045
#> 17         1 2239.96841    final  0.1219045
#> 18         1 2297.60757    final  0.1219045
#> 19         1 2570.13516    final  0.1219045
#> 20         1 2701.35146    final  0.1219045
#> 21         1 2710.41060    final  0.1219045
#> 22         2   20.51980    final -0.2771686
#> 23         2  157.04543    final -0.2771686
#> 24         2  183.04004    final -0.2771686
#> 25         2  368.73227    final -0.2771686
#> 26         2  415.05424    final -0.2771686
#> 27         2  438.65781    final -0.2771686
#> 28         2  556.86775    final -0.2771686
#> 29         2  562.83489    final -0.2771686
#> 30         2  603.15873    final -0.2771686
#> 31         2  697.45272    final -0.2771686
#> 32         2  739.11078    final -0.2771686
#> 33         2  814.43262    final -0.2771686
#> 34         2  833.30126    final -0.2771686
#> 35         2  920.98665    final -0.2771686
#> 36         2  939.99043    final -0.2771686
#> 37         2 1037.89825    final -0.2771686
#> 38         2 1070.23583    final -0.2771686
#> 39         2 1202.28359    final -0.2771686
#> 40         2 1234.12965    final -0.2771686
#> 41         2 1394.63638    final -0.2771686
#> 42         2 1409.20980    final -0.2771686
#> 43         3   17.30206    final -0.2397427
#> 44         3   51.71863    final -0.2397427
#> 45         3  200.03293    final -0.2397427
#> 46         3  327.98149    final -0.2397427
#> 47         3  470.62831    final -0.2397427
#> 48         3  682.45774    final -0.2397427
#> 49         3  824.87378    final -0.2397427
#> 50         3  913.28521    final -0.2397427
#> 51         3 1182.38390    final -0.2397427
#> 52         3 1405.37461    final -0.2397427
#> 53         3 1457.65505    final -0.2397427
#> 54         3 1464.52011    final -0.2397427
#> 55         3 1521.43516    final -0.2397427
#> 56         3 1958.02626    final -0.2397427
#> 57         3 2013.83943    final -0.2397427
#> 58         3 2059.71756    final -0.2397427
#> 59         3 2085.74393    final -0.2397427
#> 60         3 2155.62811    final -0.2397427
#> 61         3 2251.94695    final -0.2397427
#> 62         3 2325.09490    final -0.2397427
#> 63         3 2436.32286    final -0.2397427

# --- filter to a specific analysis name ---
collect_results(trials, analysis = "final")
#>    replicate  timepoint analysis  mean_ctrl
#> 1          1   22.03729    final  0.1219045
#> 2          1   50.66157    final  0.1219045
#> 3          1   76.61735    final  0.1219045
#> 4          1  336.30657    final  0.1219045
#> 5          1  459.20914    final  0.1219045
#> 6          1  538.27732    final  0.1219045
#> 7          1  601.20532    final  0.1219045
#> 8          1  726.66942    final  0.1219045
#> 9          1  785.53789    final  0.1219045
#> 10         1  898.46689    final  0.1219045
#> 11         1  940.50337    final  0.1219045
#> 12         1 1661.60413    final  0.1219045
#> 13         1 1746.17633    final  0.1219045
#> 14         1 1768.73053    final  0.1219045
#> 15         1 1878.76441    final  0.1219045
#> 16         1 2103.59498    final  0.1219045
#> 17         1 2239.96841    final  0.1219045
#> 18         1 2297.60757    final  0.1219045
#> 19         1 2570.13516    final  0.1219045
#> 20         1 2701.35146    final  0.1219045
#> 21         1 2710.41060    final  0.1219045
#> 22         2   20.51980    final -0.2771686
#> 23         2  157.04543    final -0.2771686
#> 24         2  183.04004    final -0.2771686
#> 25         2  368.73227    final -0.2771686
#> 26         2  415.05424    final -0.2771686
#> 27         2  438.65781    final -0.2771686
#> 28         2  556.86775    final -0.2771686
#> 29         2  562.83489    final -0.2771686
#> 30         2  603.15873    final -0.2771686
#> 31         2  697.45272    final -0.2771686
#> 32         2  739.11078    final -0.2771686
#> 33         2  814.43262    final -0.2771686
#> 34         2  833.30126    final -0.2771686
#> 35         2  920.98665    final -0.2771686
#> 36         2  939.99043    final -0.2771686
#> 37         2 1037.89825    final -0.2771686
#> 38         2 1070.23583    final -0.2771686
#> 39         2 1202.28359    final -0.2771686
#> 40         2 1234.12965    final -0.2771686
#> 41         2 1394.63638    final -0.2771686
#> 42         2 1409.20980    final -0.2771686
#> 43         3   17.30206    final -0.2397427
#> 44         3   51.71863    final -0.2397427
#> 45         3  200.03293    final -0.2397427
#> 46         3  327.98149    final -0.2397427
#> 47         3  470.62831    final -0.2397427
#> 48         3  682.45774    final -0.2397427
#> 49         3  824.87378    final -0.2397427
#> 50         3  913.28521    final -0.2397427
#> 51         3 1182.38390    final -0.2397427
#> 52         3 1405.37461    final -0.2397427
#> 53         3 1457.65505    final -0.2397427
#> 54         3 1464.52011    final -0.2397427
#> 55         3 1521.43516    final -0.2397427
#> 56         3 1958.02626    final -0.2397427
#> 57         3 2013.83943    final -0.2397427
#> 58         3 2059.71756    final -0.2397427
#> 59         3 2085.74393    final -0.2397427
#> 60         3 2155.62811    final -0.2397427
#> 61         3 2251.94695    final -0.2397427
#> 62         3 2325.09490    final -0.2397427
#> 63         3 2436.32286    final -0.2397427
```

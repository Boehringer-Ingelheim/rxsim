# Example 3: Two arm \| Fixed design \| Two correlated continuous endpoints \| t-test

``` r
# Core package that provides Timer, Population, Trial, gen_timepoints, add_timepoints, ...
library(rxsim)

# Utilities
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(MASS)   # for mvrnorm (simulate correlated endpoints)
#> 
#> Attaching package: 'MASS'
#> The following object is masked from 'package:dplyr':
#> 
#>     select
```

## Scenario

We consider a **two-arm, fixed design** with **two correlated continuous
endpoints** per subject.  
Enrollment is piece-wise linear; analysis happens **once at the final
time**.

``` r
# Total target sample size
sample_size <- 100

# Arms and allocation (balanced)
arms       <- c("pbo", "trt")
allocation <- c(1, 1)

# Piece-wise linear enrollment and dropout rates (per unit time)
# (You can tweak these to match your program's cadence.)
enrollment <- list(
  end_time = c(4, 8, 12, 16),
  rate     = c(5, 10, 15, 20)
)
dropout <- list(
  end_time = c(4, 8, 12, 16),
  rate     = c(1, 2, 4, 8)
)

# Endpoint model parameters
# Endpoint correlation structure (common across arms for simplicity)
rho  <- 0.60            # correlation between endpoint1 and endpoint2
sd1  <- 1.00            # SD of endpoint1
sd2  <- 1.00            # SD of endpoint2

# Mean structure per arm
mu1_pbo <- 0.00         # Control mean for endpoint1
mu2_pbo <- 0.00         # Control mean for endpoint2

# Treatment effects (mean shifts)
delta1  <- 0.30         # Treatment - control difference for endpoint1
delta2  <- 0.20         # Treatment - control difference for endpoint2

mu1_trt <- mu1_pbo + delta1
mu2_trt <- mu2_pbo + delta2

# Construct covariance matrix
Sigma <- matrix(
  c(sd1^2, rho*sd1*sd2,
    rho*sd1*sd2, sd2^2),
  nrow = 2, byrow = TRUE
)

scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation = list(allocation),
  rho = rho,
  delta1 = delta1,
  delta2 = delta2
)
```

## Time points

Generate timeline of discrete timepoints and add them to a `Timer`.

``` r
timepoints <- gen_timepoints(
  sample_size = sample_size,
  arms        = arms,
  allocation  = allocation,
  enrollment  = enrollment,
  dropout     = dropout
)

t <- Timer$new(name = "trial_timer_two_endpoints")
add_timepoints(t, timepoints)

final_time <- t$get_end_timepoint()
final_time
#> [1] 12
```

## Arms (Populations)

For each arm we simulate **two correlated endpoints** per subject using
[`MASS::mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html).  
We keep **50** subjects per arm (consistent with `sample_size = 100` and
equal allocation).

``` r
# Create closure to capture parameters
mk_population_generator <- function(mu_y1, mu_y2) {
  function(n) {
    xy <- MASS::mvrnorm(
      n     = n,
      mu    = c(mu_y1, mu_y2),
      Sigma = Sigma
    )
    data.frame(
      id = seq_len(n),
      y1 = xy[, 1],
      y2 = xy[, 2],
      readout_time = 1
    )
  }
}

population_generators <- list(
  pbo = mk_population_generator(mu1_pbo, mu2_pbo),
  trt = mk_population_generator(mu1_trt, mu2_trt)
)
```

## Trial Parameters

Define enrollment/dropout profiles and parameters.

``` r
arms <- c("pbo", "trt")
allocation <- c(1, 1)
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn <- function(n) rexp(n, rate = 0.01)
```

## Trigger & Analysis

This is a **fixed design**: perform analysis **once**.  
We run two independent two-sample t-tests (one per endpoint).

``` r
analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer) {
      df_e <- df |> subset(!is.na(enroll_time))

      p1 <- t.test(y1 ~ arm, data = df_e)$p.value
      p2 <- t.test(y2 ~ arm, data = df_e)$p.value
      padj <- p.adjust(c(p1, p2), method = "holm")

      data.frame(
        scenario,
        p_y1      = unname(p1),
        p_y2      = unname(p2),
        p_holm_y1 = unname(padj[1]),
        p_holm_y2 = unname(padj[2]),
        n_total   = nrow(df_e),
        n_pbo     = sum(df_e$arm == "pbo"),
        n_trt     = sum(df_e$arm == "trt"),
        stringsAsFactors = FALSE
      )
    }
  )
)
```

## Trial

Create and run the trial.  
(We supply a `seed` for reproducibility inside the trial engine—e.g.,
enrollment/dropout sampling.)

``` r
trials <- replicate_trial(
  trial_name = "two_endpoints_fixed_design",
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

``` r
run_trials(trials)
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: two_endpoints_fixed_design_1
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
#>     name: two_endpoints_fixed_design_2
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
#>     name: two_endpoints_fixed_design_3
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
```

## Results

One row per replicate, row-bound across all replicates:

``` r
replicate_results <- do.call(rbind, lapply(seq_along(trials), function(i) {
  data.frame(
    replicate = i,
    trials[[i]]$results[[1]]$final,
    stringsAsFactors = FALSE
  )
}))

replicate_results
#>   replicate sample_size allocation rho delta1 delta2       p_y1      p_y2
#> 1         1         100       1, 1 0.6    0.3    0.2 0.07643188 0.6506456
#> 2         2         100       1, 1 0.6    0.3    0.2 0.01064244 0.2144255
#> 3         3         100       1, 1 0.6    0.3    0.2 0.17952451 0.3090734
#>    p_holm_y1 p_holm_y2 n_total n_pbo n_trt
#> 1 0.15286377 0.6506456     100    50    50
#> 2 0.02128487 0.2144255     100    50    50
#> 3 0.35904902 0.3590490     100    50    50
```

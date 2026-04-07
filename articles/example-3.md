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

Many trials evaluate co-primary or key secondary endpoints
simultaneously. Ignoring endpoint correlation when planning a study
leads to incorrect power estimates. This example shows how to simulate
two correlated continuous endpoints (e.g., a primary efficacy score and
a key secondary biomarker) and apply multiplicity adjustment via Holm’s
procedure. See [Example
1](https://boehringer-ingelheim.github.io/rxsim/articles/example-1.md)
for the simpler single-endpoint baseline.

## Scenario

We consider a **two-arm, fixed design** with **two correlated continuous
endpoints** per subject. Enrollment is piece-wise linear; analysis
happens **once at the final time**.

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

`rho = 0.60` represents a moderate positive correlation between the two
endpoints — plausible when both reflect the same underlying biological
process. The covariance matrix `Sigma` is constructed from the marginal
SDs and correlation using the standard identity cov(y1, y2) = rho × sd1
× sd2. Treatment shifts are delta1 = 0.30 SD for the primary endpoint
and delta2 = 0.20 SD for the secondary, reflecting a scenario where the
primary endpoint is more sensitive to treatment.

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

[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
converts piecewise-constant enrollment and dropout rate lists into a
deterministic set of calendar timepoints — in contrast to the stochastic
`rexp` calls used in [Examples
1](https://boehringer-ingelheim.github.io/rxsim/articles/example-1.md)
and
[2](https://boehringer-ingelheim.github.io/rxsim/articles/example-2.md).
`t$get_end_timepoint()` returns the latest timepoint at which all
expected enrollment and follow-up will have completed, providing the
calendar time for the single final analysis.

## Arms (Populations)

For each arm we simulate **two correlated endpoints** per subject using
[`MASS::mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html). We keep
**50** subjects per arm (consistent with `sample_size = 100` and equal
allocation).

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

[`MASS::mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) draws `n`
samples from a multivariate normal with the given mean vector and
covariance matrix `Sigma`. The closure pattern
(`mk_population_generator`) captures the arm-specific means at
definition time, so each call to `population_generators$pbo(n)` or
`population_generators$trt(n)` produces data under the correct arm
parameters. Endpoint correlation flows entirely through the shared
`Sigma` matrix — changing `rho` in the scenario automatically updates
both arms.

## Trial Parameters

Define enrollment/dropout profiles and parameters.

``` r
arms <- c("pbo", "trt")
allocation <- c(1, 1)
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn <- function(n) rexp(n, rate = 0.01)
```

Note that `enrollment_fn` and `dropout_fn` here are stochastic functions
passed to
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md),
which uses them in
[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
to generate per-replicate enrollment and dropout times. They are
separate from the piecewise-constant `enrollment` and `dropout` lists
used by
[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
above — those lists drive the deterministic timer, while these functions
drive per-subject event times within each simulated replicate.

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

The two t-tests are run independently on the enrolled subset, producing
unadjusted p-values `p_y1` and `p_y2`. `p.adjust(..., method = "holm")`
applies Holm’s stepwise procedure: `padj[1]` corresponds to the endpoint
with the smaller raw p-value (most significant) and `padj[2]` to the
larger. The adjusted values account for the familywise error rate across
the two endpoints, so `p_holm_y1 >= p_y1` always holds.

## Trial

Create and run the trial. (We supply a `seed` for reproducibility inside
the trial engine—e.g., enrollment/dropout sampling.)

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
replicate_results <- collect_results(trials)
replicate_results
#>   replicate timepoint analysis sample_size allocation rho delta1 delta2
#> 1         1  87.36847    final         100       1, 1 0.6    0.3    0.2
#> 2         2  83.92580    final         100       1, 1 0.6    0.3    0.2
#> 3         3  96.44952    final         100       1, 1 0.6    0.3    0.2
#>         p_y1      p_y2  p_holm_y1 p_holm_y2 n_total n_pbo n_trt
#> 1 0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 2 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 3 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
```

[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
prepends `replicate`, `timepoint`, and `analysis` columns to the four
p-value columns. `p_holm_y1` and `p_holm_y2` will always be ≥ their
unadjusted counterparts because Holm correction is more conservative.
Because the two endpoints are positively correlated (rho = 0.6),
replicates that reject for y1 tend to also reject for y2 — joint
rejection probability exceeds what you would expect for independent
endpoints under the same effect sizes.

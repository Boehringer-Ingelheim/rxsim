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
#>     replicate  timepoint analysis sample_size allocation rho delta1 delta2
#> 1           1   87.36847    final         100       1, 1 0.6    0.3    0.2
#> 2           1  347.10044    final         100       1, 1 0.6    0.3    0.2
#> 3           1  364.07755    final         100       1, 1 0.6    0.3    0.2
#> 4           1  448.68714    final         100       1, 1 0.6    0.3    0.2
#> 5           1  469.18697    final         100       1, 1 0.6    0.3    0.2
#> 6           1  487.15448    final         100       1, 1 0.6    0.3    0.2
#> 7           1  538.98616    final         100       1, 1 0.6    0.3    0.2
#> 8           1  606.20122    final         100       1, 1 0.6    0.3    0.2
#> 9           1  658.76511    final         100       1, 1 0.6    0.3    0.2
#> 10          1  794.98780    final         100       1, 1 0.6    0.3    0.2
#> 11          1  944.07711    final         100       1, 1 0.6    0.3    0.2
#> 12          1 1181.37320    final         100       1, 1 0.6    0.3    0.2
#> 13          1 1212.56973    final         100       1, 1 0.6    0.3    0.2
#> 14          1 1232.97050    final         100       1, 1 0.6    0.3    0.2
#> 15          1 1264.36967    final         100       1, 1 0.6    0.3    0.2
#> 16          1 1365.41125    final         100       1, 1 0.6    0.3    0.2
#> 17          1 1414.98418    final         100       1, 1 0.6    0.3    0.2
#> 18          1 1580.61613    final         100       1, 1 0.6    0.3    0.2
#> 19          1 1905.27865    final         100       1, 1 0.6    0.3    0.2
#> 20          1 1928.66156    final         100       1, 1 0.6    0.3    0.2
#> 21          1 1948.35656    final         100       1, 1 0.6    0.3    0.2
#> 22          1 2080.41273    final         100       1, 1 0.6    0.3    0.2
#> 23          1 2116.29566    final         100       1, 1 0.6    0.3    0.2
#> 24          1 2471.76965    final         100       1, 1 0.6    0.3    0.2
#> 25          1 2506.45921    final         100       1, 1 0.6    0.3    0.2
#> 26          1 2624.28188    final         100       1, 1 0.6    0.3    0.2
#> 27          1 2670.38674    final         100       1, 1 0.6    0.3    0.2
#> 28          1 2721.59496    final         100       1, 1 0.6    0.3    0.2
#> 29          1 2730.16239    final         100       1, 1 0.6    0.3    0.2
#> 30          1 2897.71439    final         100       1, 1 0.6    0.3    0.2
#> 31          1 2943.10416    final         100       1, 1 0.6    0.3    0.2
#> 32          1 2970.64127    final         100       1, 1 0.6    0.3    0.2
#> 33          1 3098.51997    final         100       1, 1 0.6    0.3    0.2
#> 34          1 3101.99965    final         100       1, 1 0.6    0.3    0.2
#> 35          1 3117.23713    final         100       1, 1 0.6    0.3    0.2
#> 36          1 3264.59026    final         100       1, 1 0.6    0.3    0.2
#> 37          1 3324.84528    final         100       1, 1 0.6    0.3    0.2
#> 38          1 3376.45601    final         100       1, 1 0.6    0.3    0.2
#> 39          1 3382.96904    final         100       1, 1 0.6    0.3    0.2
#> 40          1 3392.33000    final         100       1, 1 0.6    0.3    0.2
#> 41          1 3653.75656    final         100       1, 1 0.6    0.3    0.2
#> 42          1 3778.41117    final         100       1, 1 0.6    0.3    0.2
#> 43          1 3954.92212    final         100       1, 1 0.6    0.3    0.2
#> 44          1 3993.06729    final         100       1, 1 0.6    0.3    0.2
#> 45          1 4028.10899    final         100       1, 1 0.6    0.3    0.2
#> 46          1 4041.70989    final         100       1, 1 0.6    0.3    0.2
#> 47          1 4165.86312    final         100       1, 1 0.6    0.3    0.2
#> 48          1 4195.93190    final         100       1, 1 0.6    0.3    0.2
#> 49          1 4258.85114    final         100       1, 1 0.6    0.3    0.2
#> 50          1 4480.33893    final         100       1, 1 0.6    0.3    0.2
#> 51          1 4485.70500    final         100       1, 1 0.6    0.3    0.2
#> 52          1 4538.39996    final         100       1, 1 0.6    0.3    0.2
#> 53          1 4621.85701    final         100       1, 1 0.6    0.3    0.2
#> 54          1 4980.81111    final         100       1, 1 0.6    0.3    0.2
#> 55          1 5022.77206    final         100       1, 1 0.6    0.3    0.2
#> 56          1 5025.26114    final         100       1, 1 0.6    0.3    0.2
#> 57          1 5283.65388    final         100       1, 1 0.6    0.3    0.2
#> 58          1 5322.97061    final         100       1, 1 0.6    0.3    0.2
#> 59          1 5454.57682    final         100       1, 1 0.6    0.3    0.2
#> 60          1 5667.53959    final         100       1, 1 0.6    0.3    0.2
#> 61          1 5807.36313    final         100       1, 1 0.6    0.3    0.2
#> 62          1 5822.36146    final         100       1, 1 0.6    0.3    0.2
#> 63          1 6024.11144    final         100       1, 1 0.6    0.3    0.2
#> 64          1 6033.82973    final         100       1, 1 0.6    0.3    0.2
#> 65          1 6099.33001    final         100       1, 1 0.6    0.3    0.2
#> 66          1 6240.00557    final         100       1, 1 0.6    0.3    0.2
#> 67          1 6313.79473    final         100       1, 1 0.6    0.3    0.2
#> 68          1 6398.47626    final         100       1, 1 0.6    0.3    0.2
#> 69          1 6843.06418    final         100       1, 1 0.6    0.3    0.2
#> 70          1 6893.26695    final         100       1, 1 0.6    0.3    0.2
#> 71          1 6937.06974    final         100       1, 1 0.6    0.3    0.2
#> 72          1 6940.14080    final         100       1, 1 0.6    0.3    0.2
#> 73          1 6951.95172    final         100       1, 1 0.6    0.3    0.2
#> 74          1 7300.51551    final         100       1, 1 0.6    0.3    0.2
#> 75          1 7357.70574    final         100       1, 1 0.6    0.3    0.2
#> 76          1 7422.76732    final         100       1, 1 0.6    0.3    0.2
#> 77          1 7434.85486    final         100       1, 1 0.6    0.3    0.2
#> 78          1 7455.60494    final         100       1, 1 0.6    0.3    0.2
#> 79          1 7516.69768    final         100       1, 1 0.6    0.3    0.2
#> 80          1 7811.33067    final         100       1, 1 0.6    0.3    0.2
#> 81          1 7857.91280    final         100       1, 1 0.6    0.3    0.2
#> 82          1 7997.65149    final         100       1, 1 0.6    0.3    0.2
#> 83          1 8034.86411    final         100       1, 1 0.6    0.3    0.2
#> 84          1 8063.25236    final         100       1, 1 0.6    0.3    0.2
#> 85          1 8163.50249    final         100       1, 1 0.6    0.3    0.2
#> 86          1 8287.81000    final         100       1, 1 0.6    0.3    0.2
#> 87          1 8329.31008    final         100       1, 1 0.6    0.3    0.2
#> 88          1 8368.84116    final         100       1, 1 0.6    0.3    0.2
#> 89          1 8400.65785    final         100       1, 1 0.6    0.3    0.2
#> 90          1 8571.29608    final         100       1, 1 0.6    0.3    0.2
#> 91          1 8584.67803    final         100       1, 1 0.6    0.3    0.2
#> 92          1 8625.88371    final         100       1, 1 0.6    0.3    0.2
#> 93          1 8644.56538    final         100       1, 1 0.6    0.3    0.2
#> 94          1 8700.35951    final         100       1, 1 0.6    0.3    0.2
#> 95          1 8828.78289    final         100       1, 1 0.6    0.3    0.2
#> 96          1 8850.81043    final         100       1, 1 0.6    0.3    0.2
#> 97          1 9065.76198    final         100       1, 1 0.6    0.3    0.2
#> 98          1 9187.05539    final         100       1, 1 0.6    0.3    0.2
#> 99          2   83.92580    final         100       1, 1 0.6    0.3    0.2
#> 100         2  162.86506    final         100       1, 1 0.6    0.3    0.2
#> 101         2  164.52703    final         100       1, 1 0.6    0.3    0.2
#> 102         2  409.86390    final         100       1, 1 0.6    0.3    0.2
#> 103         2  410.46923    final         100       1, 1 0.6    0.3    0.2
#> 104         2  711.15813    final         100       1, 1 0.6    0.3    0.2
#> 105         2  775.96471    final         100       1, 1 0.6    0.3    0.2
#> 106         2  859.67365    final         100       1, 1 0.6    0.3    0.2
#> 107         2 1025.88990    final         100       1, 1 0.6    0.3    0.2
#> 108         2 1131.25232    final         100       1, 1 0.6    0.3    0.2
#> 109         2 1156.72552    final         100       1, 1 0.6    0.3    0.2
#> 110         2 1221.70859    final         100       1, 1 0.6    0.3    0.2
#> 111         2 1376.43033    final         100       1, 1 0.6    0.3    0.2
#> 112         2 1382.39410    final         100       1, 1 0.6    0.3    0.2
#> 113         2 1383.07312    final         100       1, 1 0.6    0.3    0.2
#> 114         2 1452.59714    final         100       1, 1 0.6    0.3    0.2
#> 115         2 1458.82499    final         100       1, 1 0.6    0.3    0.2
#> 116         2 1487.90153    final         100       1, 1 0.6    0.3    0.2
#> 117         2 1574.45724    final         100       1, 1 0.6    0.3    0.2
#> 118         2 1687.75601    final         100       1, 1 0.6    0.3    0.2
#> 119         2 1885.57444    final         100       1, 1 0.6    0.3    0.2
#> 120         2 1930.83150    final         100       1, 1 0.6    0.3    0.2
#> 121         2 2103.96690    final         100       1, 1 0.6    0.3    0.2
#> 122         2 2266.56908    final         100       1, 1 0.6    0.3    0.2
#> 123         2 2482.22496    final         100       1, 1 0.6    0.3    0.2
#> 124         2 2829.74784    final         100       1, 1 0.6    0.3    0.2
#> 125         2 2876.89371    final         100       1, 1 0.6    0.3    0.2
#> 126         2 2904.98490    final         100       1, 1 0.6    0.3    0.2
#> 127         2 2982.17747    final         100       1, 1 0.6    0.3    0.2
#> 128         2 3051.66387    final         100       1, 1 0.6    0.3    0.2
#> 129         2 3207.73134    final         100       1, 1 0.6    0.3    0.2
#> 130         2 3233.42137    final         100       1, 1 0.6    0.3    0.2
#> 131         2 3329.18229    final         100       1, 1 0.6    0.3    0.2
#> 132         2 3335.55715    final         100       1, 1 0.6    0.3    0.2
#> 133         2 3383.11618    final         100       1, 1 0.6    0.3    0.2
#> 134         2 3397.06429    final         100       1, 1 0.6    0.3    0.2
#> 135         2 3436.36201    final         100       1, 1 0.6    0.3    0.2
#> 136         2 3476.51551    final         100       1, 1 0.6    0.3    0.2
#> 137         2 3781.58448    final         100       1, 1 0.6    0.3    0.2
#> 138         2 3890.17196    final         100       1, 1 0.6    0.3    0.2
#> 139         2 3907.82140    final         100       1, 1 0.6    0.3    0.2
#> 140         2 3957.68072    final         100       1, 1 0.6    0.3    0.2
#> 141         2 4046.83243    final         100       1, 1 0.6    0.3    0.2
#> 142         2 4070.13216    final         100       1, 1 0.6    0.3    0.2
#> 143         2 4236.18139    final         100       1, 1 0.6    0.3    0.2
#> 144         2 4366.13531    final         100       1, 1 0.6    0.3    0.2
#> 145         2 5058.56822    final         100       1, 1 0.6    0.3    0.2
#> 146         2 5195.06935    final         100       1, 1 0.6    0.3    0.2
#> 147         2 5247.85622    final         100       1, 1 0.6    0.3    0.2
#> 148         2 5306.75175    final         100       1, 1 0.6    0.3    0.2
#> 149         2 5447.39854    final         100       1, 1 0.6    0.3    0.2
#> 150         2 5453.48043    final         100       1, 1 0.6    0.3    0.2
#> 151         2 5578.89134    final         100       1, 1 0.6    0.3    0.2
#> 152         2 5613.89604    final         100       1, 1 0.6    0.3    0.2
#> 153         2 5872.18930    final         100       1, 1 0.6    0.3    0.2
#> 154         2 6068.78939    final         100       1, 1 0.6    0.3    0.2
#> 155         2 6175.40092    final         100       1, 1 0.6    0.3    0.2
#> 156         2 6190.98047    final         100       1, 1 0.6    0.3    0.2
#> 157         2 6229.16404    final         100       1, 1 0.6    0.3    0.2
#> 158         2 6323.41052    final         100       1, 1 0.6    0.3    0.2
#> 159         2 6401.16754    final         100       1, 1 0.6    0.3    0.2
#> 160         2 6496.77482    final         100       1, 1 0.6    0.3    0.2
#> 161         2 6560.01235    final         100       1, 1 0.6    0.3    0.2
#> 162         2 6569.19544    final         100       1, 1 0.6    0.3    0.2
#> 163         2 6571.86589    final         100       1, 1 0.6    0.3    0.2
#> 164         2 6807.20046    final         100       1, 1 0.6    0.3    0.2
#> 165         2 6896.87699    final         100       1, 1 0.6    0.3    0.2
#> 166         2 7109.88073    final         100       1, 1 0.6    0.3    0.2
#> 167         2 7124.19290    final         100       1, 1 0.6    0.3    0.2
#> 168         2 7159.79051    final         100       1, 1 0.6    0.3    0.2
#> 169         2 7197.93433    final         100       1, 1 0.6    0.3    0.2
#> 170         2 7316.61338    final         100       1, 1 0.6    0.3    0.2
#> 171         2 7581.32242    final         100       1, 1 0.6    0.3    0.2
#> 172         2 7847.68211    final         100       1, 1 0.6    0.3    0.2
#> 173         2 7902.64242    final         100       1, 1 0.6    0.3    0.2
#> 174         2 8062.46831    final         100       1, 1 0.6    0.3    0.2
#> 175         2 8189.62566    final         100       1, 1 0.6    0.3    0.2
#> 176         2 8233.68681    final         100       1, 1 0.6    0.3    0.2
#> 177         2 8338.28947    final         100       1, 1 0.6    0.3    0.2
#> 178         2 8346.03816    final         100       1, 1 0.6    0.3    0.2
#> 179         2 8349.95850    final         100       1, 1 0.6    0.3    0.2
#> 180         2 8359.23202    final         100       1, 1 0.6    0.3    0.2
#> 181         2 8429.83124    final         100       1, 1 0.6    0.3    0.2
#> 182         2 8441.34870    final         100       1, 1 0.6    0.3    0.2
#> 183         2 8555.06586    final         100       1, 1 0.6    0.3    0.2
#> 184         2 8585.46590    final         100       1, 1 0.6    0.3    0.2
#> 185         2 8614.05634    final         100       1, 1 0.6    0.3    0.2
#> 186         2 8780.06920    final         100       1, 1 0.6    0.3    0.2
#> 187         2 8849.20991    final         100       1, 1 0.6    0.3    0.2
#> 188         2 8851.85363    final         100       1, 1 0.6    0.3    0.2
#> 189         2 8867.06403    final         100       1, 1 0.6    0.3    0.2
#> 190         2 8968.75154    final         100       1, 1 0.6    0.3    0.2
#> 191         2 9149.58844    final         100       1, 1 0.6    0.3    0.2
#> 192         2 9263.90968    final         100       1, 1 0.6    0.3    0.2
#> 193         2 9389.80570    final         100       1, 1 0.6    0.3    0.2
#> 194         2 9399.46535    final         100       1, 1 0.6    0.3    0.2
#> 195         2 9476.42755    final         100       1, 1 0.6    0.3    0.2
#> 196         2 9483.19609    final         100       1, 1 0.6    0.3    0.2
#> 197         2 9531.04625    final         100       1, 1 0.6    0.3    0.2
#> 198         2 9574.95081    final         100       1, 1 0.6    0.3    0.2
#> 199         2 9650.05829    final         100       1, 1 0.6    0.3    0.2
#> 200         3   96.44952    final         100       1, 1 0.6    0.3    0.2
#> 201         3  641.30535    final         100       1, 1 0.6    0.3    0.2
#> 202         3  642.31945    final         100       1, 1 0.6    0.3    0.2
#> 203         3  674.93478    final         100       1, 1 0.6    0.3    0.2
#> 204         3  769.22693    final         100       1, 1 0.6    0.3    0.2
#> 205         3  842.08596    final         100       1, 1 0.6    0.3    0.2
#> 206         3  875.35735    final         100       1, 1 0.6    0.3    0.2
#> 207         3  907.68533    final         100       1, 1 0.6    0.3    0.2
#> 208         3  997.34185    final         100       1, 1 0.6    0.3    0.2
#> 209         3 1072.05718    final         100       1, 1 0.6    0.3    0.2
#> 210         3 1163.39035    final         100       1, 1 0.6    0.3    0.2
#> 211         3 1191.59587    final         100       1, 1 0.6    0.3    0.2
#> 212         3 1562.07009    final         100       1, 1 0.6    0.3    0.2
#> 213         3 1899.90611    final         100       1, 1 0.6    0.3    0.2
#> 214         3 1940.35632    final         100       1, 1 0.6    0.3    0.2
#> 215         3 2154.49290    final         100       1, 1 0.6    0.3    0.2
#> 216         3 2447.37535    final         100       1, 1 0.6    0.3    0.2
#> 217         3 2822.93260    final         100       1, 1 0.6    0.3    0.2
#> 218         3 2941.35975    final         100       1, 1 0.6    0.3    0.2
#> 219         3 3110.94827    final         100       1, 1 0.6    0.3    0.2
#> 220         3 3148.98925    final         100       1, 1 0.6    0.3    0.2
#> 221         3 3217.53725    final         100       1, 1 0.6    0.3    0.2
#> 222         3 3337.18543    final         100       1, 1 0.6    0.3    0.2
#> 223         3 3367.90075    final         100       1, 1 0.6    0.3    0.2
#> 224         3 3529.54127    final         100       1, 1 0.6    0.3    0.2
#> 225         3 3561.12345    final         100       1, 1 0.6    0.3    0.2
#> 226         3 3709.40028    final         100       1, 1 0.6    0.3    0.2
#> 227         3 3840.98646    final         100       1, 1 0.6    0.3    0.2
#> 228         3 3981.35337    final         100       1, 1 0.6    0.3    0.2
#> 229         3 4027.76788    final         100       1, 1 0.6    0.3    0.2
#> 230         3 4160.72530    final         100       1, 1 0.6    0.3    0.2
#> 231         3 4186.89686    final         100       1, 1 0.6    0.3    0.2
#> 232         3 4399.55299    final         100       1, 1 0.6    0.3    0.2
#> 233         3 4557.13615    final         100       1, 1 0.6    0.3    0.2
#> 234         3 4638.33157    final         100       1, 1 0.6    0.3    0.2
#> 235         3 4639.82467    final         100       1, 1 0.6    0.3    0.2
#> 236         3 4671.16697    final         100       1, 1 0.6    0.3    0.2
#> 237         3 4728.99685    final         100       1, 1 0.6    0.3    0.2
#> 238         3 4770.76598    final         100       1, 1 0.6    0.3    0.2
#> 239         3 4982.05892    final         100       1, 1 0.6    0.3    0.2
#> 240         3 4984.32119    final         100       1, 1 0.6    0.3    0.2
#> 241         3 5183.69330    final         100       1, 1 0.6    0.3    0.2
#> 242         3 5246.51642    final         100       1, 1 0.6    0.3    0.2
#> 243         3 5277.27244    final         100       1, 1 0.6    0.3    0.2
#> 244         3 5387.88795    final         100       1, 1 0.6    0.3    0.2
#> 245         3 5407.20693    final         100       1, 1 0.6    0.3    0.2
#> 246         3 5452.74780    final         100       1, 1 0.6    0.3    0.2
#> 247         3 5468.71062    final         100       1, 1 0.6    0.3    0.2
#> 248         3 5542.58750    final         100       1, 1 0.6    0.3    0.2
#> 249         3 5763.47315    final         100       1, 1 0.6    0.3    0.2
#> 250         3 5771.96148    final         100       1, 1 0.6    0.3    0.2
#> 251         3 5984.73519    final         100       1, 1 0.6    0.3    0.2
#> 252         3 6004.08873    final         100       1, 1 0.6    0.3    0.2
#> 253         3 6300.91345    final         100       1, 1 0.6    0.3    0.2
#> 254         3 6311.73592    final         100       1, 1 0.6    0.3    0.2
#> 255         3 6313.31644    final         100       1, 1 0.6    0.3    0.2
#> 256         3 6417.45507    final         100       1, 1 0.6    0.3    0.2
#> 257         3 6569.04139    final         100       1, 1 0.6    0.3    0.2
#> 258         3 6679.38170    final         100       1, 1 0.6    0.3    0.2
#> 259         3 6809.27436    final         100       1, 1 0.6    0.3    0.2
#> 260         3 6813.54794    final         100       1, 1 0.6    0.3    0.2
#> 261         3 6849.44519    final         100       1, 1 0.6    0.3    0.2
#> 262         3 6970.02494    final         100       1, 1 0.6    0.3    0.2
#> 263         3 6987.69901    final         100       1, 1 0.6    0.3    0.2
#> 264         3 7057.59733    final         100       1, 1 0.6    0.3    0.2
#> 265         3 7120.24351    final         100       1, 1 0.6    0.3    0.2
#> 266         3 7148.95843    final         100       1, 1 0.6    0.3    0.2
#> 267         3 7367.40178    final         100       1, 1 0.6    0.3    0.2
#> 268         3 7413.14504    final         100       1, 1 0.6    0.3    0.2
#> 269         3 7447.56728    final         100       1, 1 0.6    0.3    0.2
#> 270         3 7511.52748    final         100       1, 1 0.6    0.3    0.2
#> 271         3 7545.62878    final         100       1, 1 0.6    0.3    0.2
#> 272         3 7575.01451    final         100       1, 1 0.6    0.3    0.2
#> 273         3 7653.02453    final         100       1, 1 0.6    0.3    0.2
#> 274         3 7718.38987    final         100       1, 1 0.6    0.3    0.2
#> 275         3 7742.26552    final         100       1, 1 0.6    0.3    0.2
#> 276         3 7810.47626    final         100       1, 1 0.6    0.3    0.2
#> 277         3 7929.15282    final         100       1, 1 0.6    0.3    0.2
#> 278         3 8138.48451    final         100       1, 1 0.6    0.3    0.2
#> 279         3 8213.91817    final         100       1, 1 0.6    0.3    0.2
#> 280         3 8356.04399    final         100       1, 1 0.6    0.3    0.2
#> 281         3 8394.89321    final         100       1, 1 0.6    0.3    0.2
#> 282         3 8567.68352    final         100       1, 1 0.6    0.3    0.2
#> 283         3 8612.32526    final         100       1, 1 0.6    0.3    0.2
#> 284         3 8683.38045    final         100       1, 1 0.6    0.3    0.2
#> 285         3 8693.11929    final         100       1, 1 0.6    0.3    0.2
#> 286         3 8832.87024    final         100       1, 1 0.6    0.3    0.2
#> 287         3 8900.00150    final         100       1, 1 0.6    0.3    0.2
#> 288         3 8988.58802    final         100       1, 1 0.6    0.3    0.2
#> 289         3 9110.43423    final         100       1, 1 0.6    0.3    0.2
#> 290         3 9272.97598    final         100       1, 1 0.6    0.3    0.2
#> 291         3 9283.41623    final         100       1, 1 0.6    0.3    0.2
#> 292         3 9313.77786    final         100       1, 1 0.6    0.3    0.2
#> 293         3 9394.88550    final         100       1, 1 0.6    0.3    0.2
#> 294         3 9562.00349    final         100       1, 1 0.6    0.3    0.2
#> 295         3 9594.28157    final         100       1, 1 0.6    0.3    0.2
#> 296         3 9701.40004    final         100       1, 1 0.6    0.3    0.2
#> 297         3 9706.49654    final         100       1, 1 0.6    0.3    0.2
#> 298         3 9730.28383    final         100       1, 1 0.6    0.3    0.2
#>           p_y1      p_y2  p_holm_y1 p_holm_y2 n_total n_pbo n_trt
#> 1   0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 2   0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 3   0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 4   0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 5   0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 6   0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 7   0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 8   0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 9   0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 10  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 11  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 12  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 13  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 14  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 15  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 16  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 17  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 18  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 19  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 20  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 21  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 22  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 23  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 24  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 25  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 26  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 27  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 28  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 29  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 30  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 31  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 32  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 33  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 34  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 35  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 36  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 37  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 38  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 39  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 40  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 41  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 42  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 43  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 44  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 45  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 46  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 47  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 48  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 49  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 50  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 51  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 52  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 53  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 54  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 55  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 56  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 57  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 58  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 59  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 60  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 61  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 62  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 63  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 64  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 65  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 66  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 67  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 68  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 69  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 70  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 71  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 72  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 73  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 74  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 75  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 76  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 77  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 78  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 79  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 80  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 81  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 82  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 83  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 84  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 85  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 86  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 87  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 88  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 89  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 90  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 91  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 92  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 93  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 94  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 95  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 96  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 97  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 98  0.07643188 0.6506456 0.15286377 0.6506456     100    50    50
#> 99  0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 100 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 101 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 102 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 103 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 104 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 105 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 106 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 107 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 108 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 109 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 110 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 111 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 112 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 113 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 114 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 115 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 116 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 117 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 118 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 119 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 120 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 121 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 122 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 123 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 124 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 125 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 126 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 127 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 128 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 129 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 130 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 131 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 132 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 133 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 134 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 135 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 136 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 137 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 138 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 139 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 140 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 141 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 142 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 143 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 144 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 145 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 146 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 147 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 148 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 149 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 150 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 151 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 152 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 153 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 154 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 155 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 156 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 157 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 158 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 159 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 160 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 161 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 162 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 163 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 164 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 165 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 166 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 167 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 168 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 169 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 170 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 171 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 172 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 173 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 174 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 175 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 176 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 177 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 178 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 179 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 180 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 181 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 182 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 183 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 184 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 185 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 186 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 187 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 188 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 189 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 190 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 191 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 192 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 193 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 194 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 195 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 196 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 197 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 198 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 199 0.01064244 0.2144255 0.02128487 0.2144255     100    50    50
#> 200 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 201 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 202 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 203 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 204 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 205 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 206 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 207 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 208 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 209 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 210 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 211 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 212 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 213 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 214 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 215 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 216 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 217 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 218 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 219 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 220 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 221 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 222 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 223 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 224 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 225 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 226 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 227 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 228 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 229 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 230 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 231 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 232 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 233 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 234 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 235 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 236 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 237 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 238 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 239 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 240 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 241 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 242 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 243 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 244 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 245 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 246 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 247 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 248 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 249 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 250 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 251 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 252 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 253 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 254 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 255 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 256 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 257 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 258 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 259 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 260 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 261 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 262 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 263 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 264 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 265 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 266 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 267 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 268 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 269 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 270 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 271 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 272 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 273 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 274 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 275 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 276 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 277 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 278 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 279 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 280 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 281 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 282 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 283 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 284 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 285 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 286 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 287 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 288 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 289 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 290 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 291 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 292 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 293 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 294 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 295 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 296 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 297 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
#> 298 0.17952451 0.3090734 0.35904902 0.3590490     100    50    50
```

[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
prepends `replicate`, `timepoint`, and `analysis` columns to the four
p-value columns. `p_holm_y1` and `p_holm_y2` will always be ≥ their
unadjusted counterparts because Holm correction is more conservative.
Because the two endpoints are positively correlated (rho = 0.6),
replicates that reject for y1 tend to also reject for y2 — joint
rejection probability exceeds what you would expect for independent
endpoints under the same effect sizes.

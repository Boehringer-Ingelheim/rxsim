# Example 4: Two arm \| Fixed design \| Correlated continuous & time-to-event endpoints \| t-test + survival analysis

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
library(MASS)      # for mvrnorm (simulate correlated latent variables)
#> 
#> Attaching package: 'MASS'
#> The following object is masked from 'package:dplyr':
#> 
#>     select
library(survival)  # for Surv, survdiff, coxph
```

Oncology and cardiovascular trials often co-primary a continuous
biomarker alongside a time-to-event (TTE) endpoint such as
progression-free survival or overall survival. Correlation between a
biomarker and TTE arises naturally through shared biological mechanisms.
This example demonstrates simulating such a joint endpoint structure and
applying both a t-test (for the continuous component) and a log-rank
test plus Cox proportional hazards model (for TTE) at a single final
analysis. See [Example
3](https://boehringer-ingelheim.github.io/rxsim/articles/example-3.md)
for the analogous two-correlated-continuous setup.

## Scenario

We consider a **two-arm fixed design** with a **continuous endpoint
(Y)** and a **time-to-event (TTE)** endpoint per subject. The endpoints
are **correlated** via a shared Gaussian latent structure. Enrollment is
piece-wise linear; analysis happens **once at the final time**.

``` r
# Total target N (across both arms)
sample_size <- 1200

# Arms and allocation (balanced)
arms       <- c("pbo", "trt")
allocation <- c(1, 1)

enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn <- function(n) rexp(n, rate = 0.01)

# Correlation between continuous endpoint driver (Z1) and frailty for TTE (Z2)
rho <- 0.50

# Continuous endpoint parameters
mu_pbo <- 0.0
delta  <- 0.3      # true treatment - placebo mean difference
sd_y   <- 1.0

# TTE model parameters (Exponential PH with subject frailty exp(eta))
lambda0 <- 0.08    # baseline hazard (per time unit)
HR_trt  <- 0.70    # treatment hazard ratio (<1 favors treatment)

# Frailty scale (eta = sigma_frailty * Z2)
sigma_frailty <- 0.6

# Covariance of (Z1, Z2)
SigmaZ <- matrix(c(1, rho, rho, 1), nrow = 2)

scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation  = list(allocation),
  delta       = delta,
  rho         = rho,
  hr_trt_true = HR_trt
)
```

`HR_trt = 0.70` encodes a 30% hazard reduction in the treatment arm — a
clinically meaningful effect in many oncology settings. `lambda0 = 0.08`
is the baseline event rate per time unit; combined with n=1200, this
ensures a substantial number of TTE events even with only 3 replicates.
`sigma_frailty = 0.6` controls the strength of the subject-level frailty
term: larger values introduce more between-subject heterogeneity in TTE,
making the frailty the shared biological driver linking the continuous
and TTE endpoints through the bivariate latent (Z1, Z2) with
`rho = 0.50`.

## Time points

Generate the discrete timepoints and add them to a `Timer`.

``` r
plan <- gen_plan(
  sample_size = sample_size,
  arms        = arms,
  allocation  = allocation,
  enrollment  = enrollment_fn,
  dropout     = dropout_fn
)
```

[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)
is called with `enrollment_fn` as a stochastic function, so each
invocation by
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
generates a fresh per-replicate enrollment schedule. This differs from
the piecewise-constant list approach in [Example
3](https://boehringer-ingelheim.github.io/rxsim/articles/example-3.md):
here there is no deterministic timer pre-built — enrollment timing
varies across replicates, which is the more realistic simulation regime
for large trials.

## Arms (Populations)

Per arm we simulate subject-level latent bivariate normals `(Z1, Z2)`
with correlation `rho`. - The **continuous** endpoint is
`Y = mu_arm + sd_y * Z1`. - The **TTE** endpoint uses a
proportional-hazards exponential model with subject-specific frailty:
`hazard_i = lambda_arm * exp(eta_i)`, where `eta_i = sigma_frailty * Z2`
and `lambda_trt = lambda0 * HR_trt`. Event time
`T_i ~ Exponential(rate = hazard_i)`.

``` r
# Create generator function for arm-specific data
mk_population_generator <- function(mu_y, lambda_arm) {
  function(n) {
    Z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = SigmaZ)
    Z1 <- Z[, 1]
    Z2 <- Z[, 2]
    y  <- mu_y + sd_y * Z1
    eta <- sigma_frailty * Z2
    rate_i <- lambda_arm * exp(eta)
    tte_true <- rexp(n, rate = rate_i)

    data.frame(
      id = seq_len(n),
      y  = y,
      tte_true = tte_true,
      readout_time = 1
    )
  }
}

lambda_pbo <- lambda0
lambda_trt <- lambda0 * HR_trt

population_generators <- list(
  pbo = mk_population_generator(mu_y = mu_pbo, lambda_arm = lambda_pbo),
  trt = mk_population_generator(mu_y = mu_pbo + delta, lambda_arm = lambda_trt)
)
```

Each subject gets a latent pair (Z1, Z2) drawn from a bivariate standard
normal with correlation `rho`. Z1 drives the continuous endpoint:
`Y = mu_arm + sd_y * Z1`. Z2 drives TTE through a multiplicative
frailty: `eta = sigma_frailty * Z2`, so
`hazard_i = lambda_arm * exp(eta_i)`, and `T_i ~ Exponential(hazard_i)`.
The shared latent structure means subjects with higher Z2 (higher
frailty, shorter TTE) also tend to have higher Z1 (higher Y), modelling
a biological scenario where the biomarker and TTE share a common
underlying driver.

## Trigger & Analysis

This is a **fixed design**: analyze **once at the final timepoint**. We
run a **two-sample t-test** for the continuous endpoint and a **log-rank
test** plus **Cox PH** for the TTE endpoint.

``` r
analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer) {
      df_e <- df[!is.na(df$enroll_time), , drop = FALSE]
      if (nrow(df_e) == 0) {
        return(data.frame(scenario, note = "no enrolled subjects", stringsAsFactors = FALSE))
      }

      # Build observed TTE with censoring
      censor_drop  <- ifelse(is.na(df_e$drop_time), Inf, pmax(0, df_e$drop_time - df_e$enroll_time))
      censor_admin <- pmax(0, max(df_e$enroll_time) + 100 - df_e$enroll_time)

      tte_true <- df_e$tte_true
      t_obs    <- pmin(tte_true, censor_drop, censor_admin)
      status   <- as.integer(tte_true <= pmin(censor_drop, censor_admin))

      # Continuous endpoint t-test
      p_t <- tryCatch(
        t.test(y ~ arm, data = df_e)$p.value,
        error = function(e) NA_real_
      )

      # Survival analysis: log-rank and Cox
      surv_p <- NA_real_
      hr <- hr_lo <- hr_hi <- NA_real_
      events <- sum(status)

      if (length(unique(df_e$arm)) == 2 && events > 0 && all(t_obs >= 0)) {
        S <- Surv(time = t_obs, event = status)

        surv_p <- tryCatch({
          sd <- survdiff(S ~ arm, data = df_e)
          1 - pchisq(sd$chisq, df = length(sd$n) - 1)
        }, error = function(e) NA_real_)

        cox <- tryCatch(coxph(S ~ arm, data = df_e), error = function(e) NULL)
        if (!is.null(cox)) {
          s <- summary(cox)
          hr    <- unname(s$coef[1, "exp(coef)"])
          hr_lo <- unname(s$conf.int[1, "lower .95"])
          hr_hi <- unname(s$conf.int[1, "upper .95"])
        }
      }

      data.frame(
        scenario,
        n_total    = nrow(df_e),
        n_pbo      = sum(df_e$arm == "pbo"),
        n_trt      = sum(df_e$arm == "trt"),
        p_ttest_y  = p_t,
        logrank_p  = surv_p,
        hr_cox     = hr,
        hr_ci_lo   = hr_lo,
        hr_ci_hi   = hr_hi,
        events     = events,
        stringsAsFactors = FALSE
      )
    }
  )
)
```

`censor_drop` censors a subject at the time elapsed between enrollment
and dropout (or infinity if no dropout). `censor_admin` provides
administrative censoring: each subject is followed for 100 time units
after the last enrollment, a common end-of-study rule. `t_obs` and
`status` capture what would actually be observed. The
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object is then
passed to [`survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html)
for the log-rank p-value and
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) for the Cox HR
and 95% CI. `tryCatch` guards against edge cases (e.g., one arm with no
events) so the simulation continues rather than stopping on an error.

## Trial

Create and run multiple trial replicates.

``` r
trials <- replicate_trial(
  trial_name = "two_endpoints_continuous_tte_fixed",
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
#>     name: two_endpoints_continuous_tte_fixed_1
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
#>     name: two_endpoints_continuous_tte_fixed_2
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
#>     name: two_endpoints_continuous_tte_fixed_3
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
#>      replicate  timepoint analysis sample_size allocation delta rho hr_trt_true
#> 1            1   1141.764    final        1200       1, 1   0.3 0.5         0.7
#> 2            1   1159.349    final        1200       1, 1   0.3 0.5         0.7
#> 3            1   1220.721    final        1200       1, 1   0.3 0.5         0.7
#> 4            1   1220.729    final        1200       1, 1   0.3 0.5         0.7
#> 5            1   1314.392    final        1200       1, 1   0.3 0.5         0.7
#> 6            1   1433.419    final        1200       1, 1   0.3 0.5         0.7
#> 7            1   1661.305    final        1200       1, 1   0.3 0.5         0.7
#> 8            1   1832.187    final        1200       1, 1   0.3 0.5         0.7
#> 9            1   1833.612    final        1200       1, 1   0.3 0.5         0.7
#> 10           1   1894.058    final        1200       1, 1   0.3 0.5         0.7
#> 11           1   1902.107    final        1200       1, 1   0.3 0.5         0.7
#> 12           1   2007.626    final        1200       1, 1   0.3 0.5         0.7
#> 13           1   2104.565    final        1200       1, 1   0.3 0.5         0.7
#> 14           1   2110.249    final        1200       1, 1   0.3 0.5         0.7
#> 15           1   2383.996    final        1200       1, 1   0.3 0.5         0.7
#> 16           1   2386.712    final        1200       1, 1   0.3 0.5         0.7
#> 17           1   2401.376    final        1200       1, 1   0.3 0.5         0.7
#> 18           1   2441.727    final        1200       1, 1   0.3 0.5         0.7
#> 19           1   2526.420    final        1200       1, 1   0.3 0.5         0.7
#> 20           1   2617.872    final        1200       1, 1   0.3 0.5         0.7
#> 21           1   2700.665    final        1200       1, 1   0.3 0.5         0.7
#> 22           1   2714.308    final        1200       1, 1   0.3 0.5         0.7
#> 23           1   2772.291    final        1200       1, 1   0.3 0.5         0.7
#> 24           1   2942.121    final        1200       1, 1   0.3 0.5         0.7
#> 25           1   3025.929    final        1200       1, 1   0.3 0.5         0.7
#> 26           1   3220.106    final        1200       1, 1   0.3 0.5         0.7
#> 27           1   3285.424    final        1200       1, 1   0.3 0.5         0.7
#> 28           1   3331.494    final        1200       1, 1   0.3 0.5         0.7
#> 29           1   3348.092    final        1200       1, 1   0.3 0.5         0.7
#> 30           1   3418.150    final        1200       1, 1   0.3 0.5         0.7
#> 31           1   3677.117    final        1200       1, 1   0.3 0.5         0.7
#> 32           1   3884.464    final        1200       1, 1   0.3 0.5         0.7
#> 33           1   4069.567    final        1200       1, 1   0.3 0.5         0.7
#> 34           1   4318.380    final        1200       1, 1   0.3 0.5         0.7
#> 35           1   4473.983    final        1200       1, 1   0.3 0.5         0.7
#> 36           1   4532.176    final        1200       1, 1   0.3 0.5         0.7
#> 37           1   5443.813    final        1200       1, 1   0.3 0.5         0.7
#> 38           1   5465.190    final        1200       1, 1   0.3 0.5         0.7
#> 39           1   5663.733    final        1200       1, 1   0.3 0.5         0.7
#> 40           1   5692.756    final        1200       1, 1   0.3 0.5         0.7
#> 41           1   5735.559    final        1200       1, 1   0.3 0.5         0.7
#> 42           1   5955.931    final        1200       1, 1   0.3 0.5         0.7
#> 43           1   6010.345    final        1200       1, 1   0.3 0.5         0.7
#> 44           1   6242.124    final        1200       1, 1   0.3 0.5         0.7
#> 45           1   6249.703    final        1200       1, 1   0.3 0.5         0.7
#> 46           1   6300.414    final        1200       1, 1   0.3 0.5         0.7
#> 47           1   6373.729    final        1200       1, 1   0.3 0.5         0.7
#> 48           1   6394.622    final        1200       1, 1   0.3 0.5         0.7
#> 49           1   6401.580    final        1200       1, 1   0.3 0.5         0.7
#> 50           1   6556.954    final        1200       1, 1   0.3 0.5         0.7
#> 51           1   6858.245    final        1200       1, 1   0.3 0.5         0.7
#> 52           1   6858.821    final        1200       1, 1   0.3 0.5         0.7
#> 53           1   6908.014    final        1200       1, 1   0.3 0.5         0.7
#> 54           1   6968.204    final        1200       1, 1   0.3 0.5         0.7
#> 55           1   6990.195    final        1200       1, 1   0.3 0.5         0.7
#> 56           1   7073.929    final        1200       1, 1   0.3 0.5         0.7
#> 57           1   7474.671    final        1200       1, 1   0.3 0.5         0.7
#> 58           1   7880.245    final        1200       1, 1   0.3 0.5         0.7
#> 59           1   7918.291    final        1200       1, 1   0.3 0.5         0.7
#> 60           1   8024.672    final        1200       1, 1   0.3 0.5         0.7
#> 61           1   8076.234    final        1200       1, 1   0.3 0.5         0.7
#> 62           1   8116.622    final        1200       1, 1   0.3 0.5         0.7
#> 63           1   8153.450    final        1200       1, 1   0.3 0.5         0.7
#> 64           1   8175.058    final        1200       1, 1   0.3 0.5         0.7
#> 65           1   8182.032    final        1200       1, 1   0.3 0.5         0.7
#> 66           1   8224.690    final        1200       1, 1   0.3 0.5         0.7
#> 67           1   8293.122    final        1200       1, 1   0.3 0.5         0.7
#> 68           1   8321.264    final        1200       1, 1   0.3 0.5         0.7
#> 69           1   8436.449    final        1200       1, 1   0.3 0.5         0.7
#> 70           1   8495.532    final        1200       1, 1   0.3 0.5         0.7
#> 71           1   8583.731    final        1200       1, 1   0.3 0.5         0.7
#> 72           1   8872.671    final        1200       1, 1   0.3 0.5         0.7
#> 73           1   8968.424    final        1200       1, 1   0.3 0.5         0.7
#> 74           1   8972.019    final        1200       1, 1   0.3 0.5         0.7
#> 75           1   8973.745    final        1200       1, 1   0.3 0.5         0.7
#> 76           1   9049.253    final        1200       1, 1   0.3 0.5         0.7
#> 77           1   9103.892    final        1200       1, 1   0.3 0.5         0.7
#> 78           1   9127.207    final        1200       1, 1   0.3 0.5         0.7
#> 79           1   9165.994    final        1200       1, 1   0.3 0.5         0.7
#> 80           1   9172.373    final        1200       1, 1   0.3 0.5         0.7
#> 81           1   9184.235    final        1200       1, 1   0.3 0.5         0.7
#> 82           1   9350.593    final        1200       1, 1   0.3 0.5         0.7
#> 83           1   9405.683    final        1200       1, 1   0.3 0.5         0.7
#> 84           1   9898.925    final        1200       1, 1   0.3 0.5         0.7
#> 85           1  10002.083    final        1200       1, 1   0.3 0.5         0.7
#> 86           1  10008.242    final        1200       1, 1   0.3 0.5         0.7
#> 87           1  10014.825    final        1200       1, 1   0.3 0.5         0.7
#> 88           1  10052.109    final        1200       1, 1   0.3 0.5         0.7
#> 89           1  10078.489    final        1200       1, 1   0.3 0.5         0.7
#> 90           1  10202.586    final        1200       1, 1   0.3 0.5         0.7
#> 91           1  10540.998    final        1200       1, 1   0.3 0.5         0.7
#> 92           1  10571.371    final        1200       1, 1   0.3 0.5         0.7
#> 93           1  10586.093    final        1200       1, 1   0.3 0.5         0.7
#> 94           1  10672.135    final        1200       1, 1   0.3 0.5         0.7
#> 95           1  10707.780    final        1200       1, 1   0.3 0.5         0.7
#> 96           1  10895.881    final        1200       1, 1   0.3 0.5         0.7
#> 97           1  10967.526    final        1200       1, 1   0.3 0.5         0.7
#> 98           1  11098.440    final        1200       1, 1   0.3 0.5         0.7
#> 99           1  11129.826    final        1200       1, 1   0.3 0.5         0.7
#> 100          1  11133.979    final        1200       1, 1   0.3 0.5         0.7
#> 101          1  11235.128    final        1200       1, 1   0.3 0.5         0.7
#> 102          1  11251.172    final        1200       1, 1   0.3 0.5         0.7
#> 103          1  11263.565    final        1200       1, 1   0.3 0.5         0.7
#> 104          1  11321.386    final        1200       1, 1   0.3 0.5         0.7
#> 105          1  11340.679    final        1200       1, 1   0.3 0.5         0.7
#> 106          1  11376.291    final        1200       1, 1   0.3 0.5         0.7
#> 107          1  11399.005    final        1200       1, 1   0.3 0.5         0.7
#> 108          1  11515.493    final        1200       1, 1   0.3 0.5         0.7
#> 109          1  11717.161    final        1200       1, 1   0.3 0.5         0.7
#> 110          1  11758.908    final        1200       1, 1   0.3 0.5         0.7
#> 111          1  11852.310    final        1200       1, 1   0.3 0.5         0.7
#> 112          1  11973.350    final        1200       1, 1   0.3 0.5         0.7
#> 113          1  12097.104    final        1200       1, 1   0.3 0.5         0.7
#> 114          1  12179.243    final        1200       1, 1   0.3 0.5         0.7
#> 115          1  12227.070    final        1200       1, 1   0.3 0.5         0.7
#> 116          1  12284.994    final        1200       1, 1   0.3 0.5         0.7
#> 117          1  12366.745    final        1200       1, 1   0.3 0.5         0.7
#> 118          1  12442.277    final        1200       1, 1   0.3 0.5         0.7
#> 119          1  12625.515    final        1200       1, 1   0.3 0.5         0.7
#> 120          1  12680.159    final        1200       1, 1   0.3 0.5         0.7
#> 121          1  12687.728    final        1200       1, 1   0.3 0.5         0.7
#> 122          1  12722.359    final        1200       1, 1   0.3 0.5         0.7
#> 123          1  12803.744    final        1200       1, 1   0.3 0.5         0.7
#> 124          1  12869.078    final        1200       1, 1   0.3 0.5         0.7
#> 125          1  12967.420    final        1200       1, 1   0.3 0.5         0.7
#> 126          1  12988.375    final        1200       1, 1   0.3 0.5         0.7
#> 127          1  13020.182    final        1200       1, 1   0.3 0.5         0.7
#> 128          1  13094.580    final        1200       1, 1   0.3 0.5         0.7
#> 129          1  13149.683    final        1200       1, 1   0.3 0.5         0.7
#> 130          1  13256.810    final        1200       1, 1   0.3 0.5         0.7
#> 131          1  13387.492    final        1200       1, 1   0.3 0.5         0.7
#> 132          1  13432.191    final        1200       1, 1   0.3 0.5         0.7
#> 133          1  13444.212    final        1200       1, 1   0.3 0.5         0.7
#> 134          1  13801.240    final        1200       1, 1   0.3 0.5         0.7
#> 135          1  13939.000    final        1200       1, 1   0.3 0.5         0.7
#> 136          1  14047.070    final        1200       1, 1   0.3 0.5         0.7
#> 137          1  14280.284    final        1200       1, 1   0.3 0.5         0.7
#> 138          1  14324.670    final        1200       1, 1   0.3 0.5         0.7
#> 139          1  14334.699    final        1200       1, 1   0.3 0.5         0.7
#> 140          1  14386.436    final        1200       1, 1   0.3 0.5         0.7
#> 141          1  14467.906    final        1200       1, 1   0.3 0.5         0.7
#> 142          1  14577.942    final        1200       1, 1   0.3 0.5         0.7
#> 143          1  14685.528    final        1200       1, 1   0.3 0.5         0.7
#> 144          1  14738.502    final        1200       1, 1   0.3 0.5         0.7
#> 145          1  15071.289    final        1200       1, 1   0.3 0.5         0.7
#> 146          1  15102.486    final        1200       1, 1   0.3 0.5         0.7
#> 147          1  15118.890    final        1200       1, 1   0.3 0.5         0.7
#> 148          1  15133.889    final        1200       1, 1   0.3 0.5         0.7
#> 149          1  15152.922    final        1200       1, 1   0.3 0.5         0.7
#> 150          1  15159.369    final        1200       1, 1   0.3 0.5         0.7
#> 151          1  15367.230    final        1200       1, 1   0.3 0.5         0.7
#> 152          1  15412.689    final        1200       1, 1   0.3 0.5         0.7
#> 153          1  15428.561    final        1200       1, 1   0.3 0.5         0.7
#> 154          1  15438.171    final        1200       1, 1   0.3 0.5         0.7
#> 155          1  15446.924    final        1200       1, 1   0.3 0.5         0.7
#> 156          1  15627.782    final        1200       1, 1   0.3 0.5         0.7
#> 157          1  15798.180    final        1200       1, 1   0.3 0.5         0.7
#> 158          1  15968.107    final        1200       1, 1   0.3 0.5         0.7
#> 159          1  16055.772    final        1200       1, 1   0.3 0.5         0.7
#> 160          1  16111.390    final        1200       1, 1   0.3 0.5         0.7
#> 161          1  16236.221    final        1200       1, 1   0.3 0.5         0.7
#> 162          1  16257.719    final        1200       1, 1   0.3 0.5         0.7
#> 163          1  16287.846    final        1200       1, 1   0.3 0.5         0.7
#> 164          1  16295.555    final        1200       1, 1   0.3 0.5         0.7
#> 165          1  16387.554    final        1200       1, 1   0.3 0.5         0.7
#> 166          1  16507.444    final        1200       1, 1   0.3 0.5         0.7
#> 167          1  16614.234    final        1200       1, 1   0.3 0.5         0.7
#> 168          1  16762.540    final        1200       1, 1   0.3 0.5         0.7
#> 169          1  16869.888    final        1200       1, 1   0.3 0.5         0.7
#> 170          1  17172.714    final        1200       1, 1   0.3 0.5         0.7
#> 171          1  17310.203    final        1200       1, 1   0.3 0.5         0.7
#> 172          1  17326.327    final        1200       1, 1   0.3 0.5         0.7
#> 173          1  17361.833    final        1200       1, 1   0.3 0.5         0.7
#> 174          1  17653.367    final        1200       1, 1   0.3 0.5         0.7
#> 175          1  17675.920    final        1200       1, 1   0.3 0.5         0.7
#> 176          1  17954.931    final        1200       1, 1   0.3 0.5         0.7
#> 177          1  18025.847    final        1200       1, 1   0.3 0.5         0.7
#> 178          1  18063.846    final        1200       1, 1   0.3 0.5         0.7
#> 179          1  18097.533    final        1200       1, 1   0.3 0.5         0.7
#> 180          1  18100.902    final        1200       1, 1   0.3 0.5         0.7
#> 181          1  18111.950    final        1200       1, 1   0.3 0.5         0.7
#> 182          1  18274.788    final        1200       1, 1   0.3 0.5         0.7
#> 183          1  18342.202    final        1200       1, 1   0.3 0.5         0.7
#> 184          1  18372.418    final        1200       1, 1   0.3 0.5         0.7
#> 185          1  18496.648    final        1200       1, 1   0.3 0.5         0.7
#> 186          1  18863.870    final        1200       1, 1   0.3 0.5         0.7
#> 187          1  18920.451    final        1200       1, 1   0.3 0.5         0.7
#> 188          1  19017.208    final        1200       1, 1   0.3 0.5         0.7
#> 189          1  19041.432    final        1200       1, 1   0.3 0.5         0.7
#> 190          1  19045.308    final        1200       1, 1   0.3 0.5         0.7
#> 191          1  19386.979    final        1200       1, 1   0.3 0.5         0.7
#> 192          1  19417.963    final        1200       1, 1   0.3 0.5         0.7
#> 193          1  19494.432    final        1200       1, 1   0.3 0.5         0.7
#> 194          1  19685.205    final        1200       1, 1   0.3 0.5         0.7
#> 195          1  19788.168    final        1200       1, 1   0.3 0.5         0.7
#> 196          1  20041.960    final        1200       1, 1   0.3 0.5         0.7
#> 197          1  20087.217    final        1200       1, 1   0.3 0.5         0.7
#> 198          1  20256.321    final        1200       1, 1   0.3 0.5         0.7
#> 199          1  20264.305    final        1200       1, 1   0.3 0.5         0.7
#> 200          1  20297.981    final        1200       1, 1   0.3 0.5         0.7
#> 201          1  20440.675    final        1200       1, 1   0.3 0.5         0.7
#> 202          1  20642.540    final        1200       1, 1   0.3 0.5         0.7
#> 203          1  20693.636    final        1200       1, 1   0.3 0.5         0.7
#> 204          1  20713.538    final        1200       1, 1   0.3 0.5         0.7
#> 205          1  20795.977    final        1200       1, 1   0.3 0.5         0.7
#> 206          1  20831.353    final        1200       1, 1   0.3 0.5         0.7
#> 207          1  20981.705    final        1200       1, 1   0.3 0.5         0.7
#> 208          1  21017.238    final        1200       1, 1   0.3 0.5         0.7
#> 209          1  21039.173    final        1200       1, 1   0.3 0.5         0.7
#> 210          1  21107.484    final        1200       1, 1   0.3 0.5         0.7
#> 211          1  21126.264    final        1200       1, 1   0.3 0.5         0.7
#> 212          1  21202.265    final        1200       1, 1   0.3 0.5         0.7
#> 213          1  21208.462    final        1200       1, 1   0.3 0.5         0.7
#> 214          1  21239.913    final        1200       1, 1   0.3 0.5         0.7
#> 215          1  21316.560    final        1200       1, 1   0.3 0.5         0.7
#> 216          1  21384.638    final        1200       1, 1   0.3 0.5         0.7
#> 217          1  21413.640    final        1200       1, 1   0.3 0.5         0.7
#> 218          1  21471.589    final        1200       1, 1   0.3 0.5         0.7
#> 219          1  21708.894    final        1200       1, 1   0.3 0.5         0.7
#> 220          1  21806.261    final        1200       1, 1   0.3 0.5         0.7
#> 221          1  21991.868    final        1200       1, 1   0.3 0.5         0.7
#> 222          1  22120.815    final        1200       1, 1   0.3 0.5         0.7
#> 223          1  22141.774    final        1200       1, 1   0.3 0.5         0.7
#> 224          1  22278.012    final        1200       1, 1   0.3 0.5         0.7
#> 225          1  22293.376    final        1200       1, 1   0.3 0.5         0.7
#> 226          1  22295.110    final        1200       1, 1   0.3 0.5         0.7
#> 227          1  22311.741    final        1200       1, 1   0.3 0.5         0.7
#> 228          1  22476.197    final        1200       1, 1   0.3 0.5         0.7
#> 229          1  22543.368    final        1200       1, 1   0.3 0.5         0.7
#> 230          1  22597.397    final        1200       1, 1   0.3 0.5         0.7
#> 231          1  22643.584    final        1200       1, 1   0.3 0.5         0.7
#> 232          1  22710.164    final        1200       1, 1   0.3 0.5         0.7
#> 233          1  23225.571    final        1200       1, 1   0.3 0.5         0.7
#> 234          1  23261.688    final        1200       1, 1   0.3 0.5         0.7
#> 235          1  23273.250    final        1200       1, 1   0.3 0.5         0.7
#> 236          1  23326.889    final        1200       1, 1   0.3 0.5         0.7
#> 237          1  23344.951    final        1200       1, 1   0.3 0.5         0.7
#> 238          1  23346.680    final        1200       1, 1   0.3 0.5         0.7
#> 239          1  23527.619    final        1200       1, 1   0.3 0.5         0.7
#> 240          1  23534.196    final        1200       1, 1   0.3 0.5         0.7
#> 241          1  23636.095    final        1200       1, 1   0.3 0.5         0.7
#> 242          1  23824.169    final        1200       1, 1   0.3 0.5         0.7
#> 243          1  23855.533    final        1200       1, 1   0.3 0.5         0.7
#> 244          1  23904.279    final        1200       1, 1   0.3 0.5         0.7
#> 245          1  23904.551    final        1200       1, 1   0.3 0.5         0.7
#> 246          1  23967.282    final        1200       1, 1   0.3 0.5         0.7
#> 247          1  24091.318    final        1200       1, 1   0.3 0.5         0.7
#> 248          1  24182.514    final        1200       1, 1   0.3 0.5         0.7
#> 249          1  24251.020    final        1200       1, 1   0.3 0.5         0.7
#> 250          1  24395.650    final        1200       1, 1   0.3 0.5         0.7
#> 251          1  24606.467    final        1200       1, 1   0.3 0.5         0.7
#> 252          1  24927.951    final        1200       1, 1   0.3 0.5         0.7
#> 253          1  24961.622    final        1200       1, 1   0.3 0.5         0.7
#> 254          1  25036.648    final        1200       1, 1   0.3 0.5         0.7
#> 255          1  25043.968    final        1200       1, 1   0.3 0.5         0.7
#> 256          1  25163.413    final        1200       1, 1   0.3 0.5         0.7
#> 257          1  25528.405    final        1200       1, 1   0.3 0.5         0.7
#> 258          1  26029.116    final        1200       1, 1   0.3 0.5         0.7
#> 259          1  26182.552    final        1200       1, 1   0.3 0.5         0.7
#> 260          1  26289.099    final        1200       1, 1   0.3 0.5         0.7
#> 261          1  26344.565    final        1200       1, 1   0.3 0.5         0.7
#> 262          1  26349.095    final        1200       1, 1   0.3 0.5         0.7
#> 263          1  26367.363    final        1200       1, 1   0.3 0.5         0.7
#> 264          1  26539.886    final        1200       1, 1   0.3 0.5         0.7
#> 265          1  26701.027    final        1200       1, 1   0.3 0.5         0.7
#> 266          1  26789.604    final        1200       1, 1   0.3 0.5         0.7
#> 267          1  26820.696    final        1200       1, 1   0.3 0.5         0.7
#> 268          1  26885.655    final        1200       1, 1   0.3 0.5         0.7
#> 269          1  26958.855    final        1200       1, 1   0.3 0.5         0.7
#> 270          1  27004.581    final        1200       1, 1   0.3 0.5         0.7
#> 271          1  27049.778    final        1200       1, 1   0.3 0.5         0.7
#> 272          1  27091.890    final        1200       1, 1   0.3 0.5         0.7
#> 273          1  27165.233    final        1200       1, 1   0.3 0.5         0.7
#> 274          1  27468.416    final        1200       1, 1   0.3 0.5         0.7
#> 275          1  27494.768    final        1200       1, 1   0.3 0.5         0.7
#> 276          1  27805.969    final        1200       1, 1   0.3 0.5         0.7
#> 277          1  28006.737    final        1200       1, 1   0.3 0.5         0.7
#> 278          1  28156.302    final        1200       1, 1   0.3 0.5         0.7
#> 279          1  28370.043    final        1200       1, 1   0.3 0.5         0.7
#> 280          1  28685.086    final        1200       1, 1   0.3 0.5         0.7
#> 281          1  28745.306    final        1200       1, 1   0.3 0.5         0.7
#> 282          1  28859.105    final        1200       1, 1   0.3 0.5         0.7
#> 283          1  28881.034    final        1200       1, 1   0.3 0.5         0.7
#> 284          1  29112.953    final        1200       1, 1   0.3 0.5         0.7
#> 285          1  29193.973    final        1200       1, 1   0.3 0.5         0.7
#> 286          1  29196.284    final        1200       1, 1   0.3 0.5         0.7
#> 287          1  29554.454    final        1200       1, 1   0.3 0.5         0.7
#> 288          1  29775.538    final        1200       1, 1   0.3 0.5         0.7
#> 289          1  29872.505    final        1200       1, 1   0.3 0.5         0.7
#> 290          1  29939.059    final        1200       1, 1   0.3 0.5         0.7
#> 291          1  30039.428    final        1200       1, 1   0.3 0.5         0.7
#> 292          1  30230.542    final        1200       1, 1   0.3 0.5         0.7
#> 293          1  30428.689    final        1200       1, 1   0.3 0.5         0.7
#> 294          1  30805.606    final        1200       1, 1   0.3 0.5         0.7
#> 295          1  30991.173    final        1200       1, 1   0.3 0.5         0.7
#> 296          1  31236.322    final        1200       1, 1   0.3 0.5         0.7
#> 297          1  31301.369    final        1200       1, 1   0.3 0.5         0.7
#> 298          1  31330.378    final        1200       1, 1   0.3 0.5         0.7
#> 299          1  31365.765    final        1200       1, 1   0.3 0.5         0.7
#> 300          1  31420.590    final        1200       1, 1   0.3 0.5         0.7
#> 301          1  31448.841    final        1200       1, 1   0.3 0.5         0.7
#> 302          1  31646.915    final        1200       1, 1   0.3 0.5         0.7
#> 303          1  31673.588    final        1200       1, 1   0.3 0.5         0.7
#> 304          1  31893.766    final        1200       1, 1   0.3 0.5         0.7
#> 305          1  32026.482    final        1200       1, 1   0.3 0.5         0.7
#> 306          1  32104.665    final        1200       1, 1   0.3 0.5         0.7
#> 307          1  32153.901    final        1200       1, 1   0.3 0.5         0.7
#> 308          1  32166.360    final        1200       1, 1   0.3 0.5         0.7
#> 309          1  32194.182    final        1200       1, 1   0.3 0.5         0.7
#> 310          1  32244.915    final        1200       1, 1   0.3 0.5         0.7
#> 311          1  32550.163    final        1200       1, 1   0.3 0.5         0.7
#> 312          1  32574.062    final        1200       1, 1   0.3 0.5         0.7
#> 313          1  32594.590    final        1200       1, 1   0.3 0.5         0.7
#> 314          1  32644.075    final        1200       1, 1   0.3 0.5         0.7
#> 315          1  32827.507    final        1200       1, 1   0.3 0.5         0.7
#> 316          1  33176.836    final        1200       1, 1   0.3 0.5         0.7
#> 317          1  33311.552    final        1200       1, 1   0.3 0.5         0.7
#> 318          1  33505.917    final        1200       1, 1   0.3 0.5         0.7
#> 319          1  33530.967    final        1200       1, 1   0.3 0.5         0.7
#> 320          1  33545.731    final        1200       1, 1   0.3 0.5         0.7
#> 321          1  33590.801    final        1200       1, 1   0.3 0.5         0.7
#> 322          1  33636.316    final        1200       1, 1   0.3 0.5         0.7
#> 323          1  33892.726    final        1200       1, 1   0.3 0.5         0.7
#> 324          1  33954.893    final        1200       1, 1   0.3 0.5         0.7
#> 325          1  34221.668    final        1200       1, 1   0.3 0.5         0.7
#> 326          1  34244.162    final        1200       1, 1   0.3 0.5         0.7
#> 327          1  34261.618    final        1200       1, 1   0.3 0.5         0.7
#> 328          1  34264.943    final        1200       1, 1   0.3 0.5         0.7
#> 329          1  34332.236    final        1200       1, 1   0.3 0.5         0.7
#> 330          1  34407.674    final        1200       1, 1   0.3 0.5         0.7
#> 331          1  34475.395    final        1200       1, 1   0.3 0.5         0.7
#> 332          1  34650.215    final        1200       1, 1   0.3 0.5         0.7
#> 333          1  34657.134    final        1200       1, 1   0.3 0.5         0.7
#> 334          1  34694.207    final        1200       1, 1   0.3 0.5         0.7
#> 335          1  34835.846    final        1200       1, 1   0.3 0.5         0.7
#> 336          1  34919.336    final        1200       1, 1   0.3 0.5         0.7
#> 337          1  34982.134    final        1200       1, 1   0.3 0.5         0.7
#> 338          1  35025.535    final        1200       1, 1   0.3 0.5         0.7
#> 339          1  35053.623    final        1200       1, 1   0.3 0.5         0.7
#> 340          1  35117.701    final        1200       1, 1   0.3 0.5         0.7
#> 341          1  35525.376    final        1200       1, 1   0.3 0.5         0.7
#> 342          1  35532.541    final        1200       1, 1   0.3 0.5         0.7
#> 343          1  35607.209    final        1200       1, 1   0.3 0.5         0.7
#> 344          1  35621.952    final        1200       1, 1   0.3 0.5         0.7
#> 345          1  35696.975    final        1200       1, 1   0.3 0.5         0.7
#> 346          1  35745.234    final        1200       1, 1   0.3 0.5         0.7
#> 347          1  35916.731    final        1200       1, 1   0.3 0.5         0.7
#> 348          1  35926.550    final        1200       1, 1   0.3 0.5         0.7
#> 349          1  36086.320    final        1200       1, 1   0.3 0.5         0.7
#> 350          1  36171.135    final        1200       1, 1   0.3 0.5         0.7
#> 351          1  36389.911    final        1200       1, 1   0.3 0.5         0.7
#> 352          1  36404.378    final        1200       1, 1   0.3 0.5         0.7
#> 353          1  36444.340    final        1200       1, 1   0.3 0.5         0.7
#> 354          1  36595.439    final        1200       1, 1   0.3 0.5         0.7
#> 355          1  36686.895    final        1200       1, 1   0.3 0.5         0.7
#> 356          1  36765.356    final        1200       1, 1   0.3 0.5         0.7
#> 357          1  36786.637    final        1200       1, 1   0.3 0.5         0.7
#> 358          1  36902.934    final        1200       1, 1   0.3 0.5         0.7
#> 359          1  36907.462    final        1200       1, 1   0.3 0.5         0.7
#> 360          1  36984.549    final        1200       1, 1   0.3 0.5         0.7
#> 361          1  37267.682    final        1200       1, 1   0.3 0.5         0.7
#> 362          1  37458.804    final        1200       1, 1   0.3 0.5         0.7
#> 363          1  37468.075    final        1200       1, 1   0.3 0.5         0.7
#> 364          1  37610.970    final        1200       1, 1   0.3 0.5         0.7
#> 365          1  37901.127    final        1200       1, 1   0.3 0.5         0.7
#> 366          1  38053.978    final        1200       1, 1   0.3 0.5         0.7
#> 367          1  38062.651    final        1200       1, 1   0.3 0.5         0.7
#> 368          1  38137.227    final        1200       1, 1   0.3 0.5         0.7
#> 369          1  38312.712    final        1200       1, 1   0.3 0.5         0.7
#> 370          1  38331.573    final        1200       1, 1   0.3 0.5         0.7
#> 371          1  38602.015    final        1200       1, 1   0.3 0.5         0.7
#> 372          1  38693.123    final        1200       1, 1   0.3 0.5         0.7
#> 373          1  38915.341    final        1200       1, 1   0.3 0.5         0.7
#> 374          1  39130.686    final        1200       1, 1   0.3 0.5         0.7
#> 375          1  39519.125    final        1200       1, 1   0.3 0.5         0.7
#> 376          1  39636.383    final        1200       1, 1   0.3 0.5         0.7
#> 377          1  39768.890    final        1200       1, 1   0.3 0.5         0.7
#> 378          1  39914.132    final        1200       1, 1   0.3 0.5         0.7
#> 379          1  39962.040    final        1200       1, 1   0.3 0.5         0.7
#> 380          1  39972.391    final        1200       1, 1   0.3 0.5         0.7
#> 381          1  40013.623    final        1200       1, 1   0.3 0.5         0.7
#> 382          1  40060.732    final        1200       1, 1   0.3 0.5         0.7
#> 383          1  40338.020    final        1200       1, 1   0.3 0.5         0.7
#> 384          1  40391.727    final        1200       1, 1   0.3 0.5         0.7
#> 385          1  40439.832    final        1200       1, 1   0.3 0.5         0.7
#> 386          1  40448.523    final        1200       1, 1   0.3 0.5         0.7
#> 387          1  40526.618    final        1200       1, 1   0.3 0.5         0.7
#> 388          1  40596.063    final        1200       1, 1   0.3 0.5         0.7
#> 389          1  40635.798    final        1200       1, 1   0.3 0.5         0.7
#> 390          1  40852.067    final        1200       1, 1   0.3 0.5         0.7
#> 391          1  40933.802    final        1200       1, 1   0.3 0.5         0.7
#> 392          1  41006.498    final        1200       1, 1   0.3 0.5         0.7
#> 393          1  41117.780    final        1200       1, 1   0.3 0.5         0.7
#> 394          1  41131.775    final        1200       1, 1   0.3 0.5         0.7
#> 395          1  41373.642    final        1200       1, 1   0.3 0.5         0.7
#> 396          1  41388.414    final        1200       1, 1   0.3 0.5         0.7
#> 397          1  41469.474    final        1200       1, 1   0.3 0.5         0.7
#> 398          1  41522.367    final        1200       1, 1   0.3 0.5         0.7
#> 399          1  41563.089    final        1200       1, 1   0.3 0.5         0.7
#> 400          1  41646.846    final        1200       1, 1   0.3 0.5         0.7
#> 401          1  41735.287    final        1200       1, 1   0.3 0.5         0.7
#> 402          1  41970.701    final        1200       1, 1   0.3 0.5         0.7
#> 403          1  41981.099    final        1200       1, 1   0.3 0.5         0.7
#> 404          1  42048.116    final        1200       1, 1   0.3 0.5         0.7
#> 405          1  42110.071    final        1200       1, 1   0.3 0.5         0.7
#> 406          1  42152.938    final        1200       1, 1   0.3 0.5         0.7
#> 407          1  42195.246    final        1200       1, 1   0.3 0.5         0.7
#> 408          1  42619.383    final        1200       1, 1   0.3 0.5         0.7
#> 409          1  42646.336    final        1200       1, 1   0.3 0.5         0.7
#> 410          1  42680.399    final        1200       1, 1   0.3 0.5         0.7
#> 411          1  42987.895    final        1200       1, 1   0.3 0.5         0.7
#> 412          1  43379.933    final        1200       1, 1   0.3 0.5         0.7
#> 413          1  43428.303    final        1200       1, 1   0.3 0.5         0.7
#> 414          1  43491.680    final        1200       1, 1   0.3 0.5         0.7
#> 415          1  43529.397    final        1200       1, 1   0.3 0.5         0.7
#> 416          1  43533.388    final        1200       1, 1   0.3 0.5         0.7
#> 417          1  43746.064    final        1200       1, 1   0.3 0.5         0.7
#> 418          1  43790.583    final        1200       1, 1   0.3 0.5         0.7
#> 419          1  43799.675    final        1200       1, 1   0.3 0.5         0.7
#> 420          1  43892.256    final        1200       1, 1   0.3 0.5         0.7
#> 421          1  43934.172    final        1200       1, 1   0.3 0.5         0.7
#> 422          1  44036.412    final        1200       1, 1   0.3 0.5         0.7
#> 423          1  44082.095    final        1200       1, 1   0.3 0.5         0.7
#> 424          1  44177.254    final        1200       1, 1   0.3 0.5         0.7
#> 425          1  44240.899    final        1200       1, 1   0.3 0.5         0.7
#> 426          1  44247.099    final        1200       1, 1   0.3 0.5         0.7
#> 427          1  44287.383    final        1200       1, 1   0.3 0.5         0.7
#> 428          1  44330.913    final        1200       1, 1   0.3 0.5         0.7
#> 429          1  44561.863    final        1200       1, 1   0.3 0.5         0.7
#> 430          1  44698.340    final        1200       1, 1   0.3 0.5         0.7
#> 431          1  44853.588    final        1200       1, 1   0.3 0.5         0.7
#> 432          1  44890.240    final        1200       1, 1   0.3 0.5         0.7
#> 433          1  45001.360    final        1200       1, 1   0.3 0.5         0.7
#> 434          1  45032.996    final        1200       1, 1   0.3 0.5         0.7
#> 435          1  45062.968    final        1200       1, 1   0.3 0.5         0.7
#> 436          1  45192.810    final        1200       1, 1   0.3 0.5         0.7
#> 437          1  45261.788    final        1200       1, 1   0.3 0.5         0.7
#> 438          1  45557.194    final        1200       1, 1   0.3 0.5         0.7
#> 439          1  45559.225    final        1200       1, 1   0.3 0.5         0.7
#> 440          1  45576.466    final        1200       1, 1   0.3 0.5         0.7
#> 441          1  45687.750    final        1200       1, 1   0.3 0.5         0.7
#> 442          1  45959.921    final        1200       1, 1   0.3 0.5         0.7
#> 443          1  45987.708    final        1200       1, 1   0.3 0.5         0.7
#> 444          1  46030.919    final        1200       1, 1   0.3 0.5         0.7
#> 445          1  46117.872    final        1200       1, 1   0.3 0.5         0.7
#> 446          1  46415.341    final        1200       1, 1   0.3 0.5         0.7
#> 447          1  46505.043    final        1200       1, 1   0.3 0.5         0.7
#> 448          1  46605.528    final        1200       1, 1   0.3 0.5         0.7
#> 449          1  46611.231    final        1200       1, 1   0.3 0.5         0.7
#> 450          1  46890.859    final        1200       1, 1   0.3 0.5         0.7
#> 451          1  46929.706    final        1200       1, 1   0.3 0.5         0.7
#> 452          1  46932.298    final        1200       1, 1   0.3 0.5         0.7
#> 453          1  47041.678    final        1200       1, 1   0.3 0.5         0.7
#> 454          1  47178.845    final        1200       1, 1   0.3 0.5         0.7
#> 455          1  47247.920    final        1200       1, 1   0.3 0.5         0.7
#> 456          1  47279.711    final        1200       1, 1   0.3 0.5         0.7
#> 457          1  47369.133    final        1200       1, 1   0.3 0.5         0.7
#> 458          1  47546.927    final        1200       1, 1   0.3 0.5         0.7
#> 459          1  47549.117    final        1200       1, 1   0.3 0.5         0.7
#> 460          1  47564.685    final        1200       1, 1   0.3 0.5         0.7
#> 461          1  47608.483    final        1200       1, 1   0.3 0.5         0.7
#> 462          1  47717.113    final        1200       1, 1   0.3 0.5         0.7
#> 463          1  47972.294    final        1200       1, 1   0.3 0.5         0.7
#> 464          1  47979.302    final        1200       1, 1   0.3 0.5         0.7
#> 465          1  48084.630    final        1200       1, 1   0.3 0.5         0.7
#> 466          1  48142.177    final        1200       1, 1   0.3 0.5         0.7
#> 467          1  48311.346    final        1200       1, 1   0.3 0.5         0.7
#> 468          1  48363.582    final        1200       1, 1   0.3 0.5         0.7
#> 469          1  48388.435    final        1200       1, 1   0.3 0.5         0.7
#> 470          1  48400.590    final        1200       1, 1   0.3 0.5         0.7
#> 471          1  48445.419    final        1200       1, 1   0.3 0.5         0.7
#> 472          1  48561.420    final        1200       1, 1   0.3 0.5         0.7
#> 473          1  48563.109    final        1200       1, 1   0.3 0.5         0.7
#> 474          1  48601.374    final        1200       1, 1   0.3 0.5         0.7
#> 475          1  48662.434    final        1200       1, 1   0.3 0.5         0.7
#> 476          1  48763.331    final        1200       1, 1   0.3 0.5         0.7
#> 477          1  48832.940    final        1200       1, 1   0.3 0.5         0.7
#> 478          1  49430.510    final        1200       1, 1   0.3 0.5         0.7
#> 479          1  49453.770    final        1200       1, 1   0.3 0.5         0.7
#> 480          1  49615.760    final        1200       1, 1   0.3 0.5         0.7
#> 481          1  49718.535    final        1200       1, 1   0.3 0.5         0.7
#> 482          1  49736.138    final        1200       1, 1   0.3 0.5         0.7
#> 483          1  49736.751    final        1200       1, 1   0.3 0.5         0.7
#> 484          1  49858.777    final        1200       1, 1   0.3 0.5         0.7
#> 485          1  49952.051    final        1200       1, 1   0.3 0.5         0.7
#> 486          1  49976.492    final        1200       1, 1   0.3 0.5         0.7
#> 487          1  50048.585    final        1200       1, 1   0.3 0.5         0.7
#> 488          1  50149.038    final        1200       1, 1   0.3 0.5         0.7
#> 489          1  50245.789    final        1200       1, 1   0.3 0.5         0.7
#> 490          1  50333.040    final        1200       1, 1   0.3 0.5         0.7
#> 491          1  50610.353    final        1200       1, 1   0.3 0.5         0.7
#> 492          1  50988.999    final        1200       1, 1   0.3 0.5         0.7
#> 493          1  51020.720    final        1200       1, 1   0.3 0.5         0.7
#> 494          1  51053.109    final        1200       1, 1   0.3 0.5         0.7
#> 495          1  51057.411    final        1200       1, 1   0.3 0.5         0.7
#> 496          1  51180.144    final        1200       1, 1   0.3 0.5         0.7
#> 497          1  51208.725    final        1200       1, 1   0.3 0.5         0.7
#> 498          1  51238.712    final        1200       1, 1   0.3 0.5         0.7
#> 499          1  51321.273    final        1200       1, 1   0.3 0.5         0.7
#> 500          1  51340.191    final        1200       1, 1   0.3 0.5         0.7
#> 501          1  51374.169    final        1200       1, 1   0.3 0.5         0.7
#> 502          1  51435.111    final        1200       1, 1   0.3 0.5         0.7
#> 503          1  51461.163    final        1200       1, 1   0.3 0.5         0.7
#> 504          1  51485.337    final        1200       1, 1   0.3 0.5         0.7
#> 505          1  51603.520    final        1200       1, 1   0.3 0.5         0.7
#> 506          1  51682.263    final        1200       1, 1   0.3 0.5         0.7
#> 507          1  51689.392    final        1200       1, 1   0.3 0.5         0.7
#> 508          1  51704.533    final        1200       1, 1   0.3 0.5         0.7
#> 509          1  51771.779    final        1200       1, 1   0.3 0.5         0.7
#> 510          1  51812.632    final        1200       1, 1   0.3 0.5         0.7
#> 511          1  51833.897    final        1200       1, 1   0.3 0.5         0.7
#> 512          1  51992.646    final        1200       1, 1   0.3 0.5         0.7
#> 513          1  52066.445    final        1200       1, 1   0.3 0.5         0.7
#> 514          1  52117.713    final        1200       1, 1   0.3 0.5         0.7
#> 515          1  52170.349    final        1200       1, 1   0.3 0.5         0.7
#> 516          1  52361.249    final        1200       1, 1   0.3 0.5         0.7
#> 517          1  52361.665    final        1200       1, 1   0.3 0.5         0.7
#> 518          1  52613.944    final        1200       1, 1   0.3 0.5         0.7
#> 519          1  52704.409    final        1200       1, 1   0.3 0.5         0.7
#> 520          1  52709.902    final        1200       1, 1   0.3 0.5         0.7
#> 521          1  52718.864    final        1200       1, 1   0.3 0.5         0.7
#> 522          1  52771.500    final        1200       1, 1   0.3 0.5         0.7
#> 523          1  52782.576    final        1200       1, 1   0.3 0.5         0.7
#> 524          1  52890.898    final        1200       1, 1   0.3 0.5         0.7
#> 525          1  52974.192    final        1200       1, 1   0.3 0.5         0.7
#> 526          1  53243.267    final        1200       1, 1   0.3 0.5         0.7
#> 527          1  53279.556    final        1200       1, 1   0.3 0.5         0.7
#> 528          1  53279.800    final        1200       1, 1   0.3 0.5         0.7
#> 529          1  53457.806    final        1200       1, 1   0.3 0.5         0.7
#> 530          1  53575.049    final        1200       1, 1   0.3 0.5         0.7
#> 531          1  53602.529    final        1200       1, 1   0.3 0.5         0.7
#> 532          1  53629.415    final        1200       1, 1   0.3 0.5         0.7
#> 533          1  53729.401    final        1200       1, 1   0.3 0.5         0.7
#> 534          1  53750.967    final        1200       1, 1   0.3 0.5         0.7
#> 535          1  53781.960    final        1200       1, 1   0.3 0.5         0.7
#> 536          1  53805.050    final        1200       1, 1   0.3 0.5         0.7
#> 537          1  54039.224    final        1200       1, 1   0.3 0.5         0.7
#> 538          1  54204.468    final        1200       1, 1   0.3 0.5         0.7
#> 539          1  54421.538    final        1200       1, 1   0.3 0.5         0.7
#> 540          1  54699.286    final        1200       1, 1   0.3 0.5         0.7
#> 541          1  54714.678    final        1200       1, 1   0.3 0.5         0.7
#> 542          1  54725.833    final        1200       1, 1   0.3 0.5         0.7
#> 543          1  54862.872    final        1200       1, 1   0.3 0.5         0.7
#> 544          1  55067.128    final        1200       1, 1   0.3 0.5         0.7
#> 545          1  55284.304    final        1200       1, 1   0.3 0.5         0.7
#> 546          1  55478.773    final        1200       1, 1   0.3 0.5         0.7
#> 547          1  55515.848    final        1200       1, 1   0.3 0.5         0.7
#> 548          1  55530.254    final        1200       1, 1   0.3 0.5         0.7
#> 549          1  55706.289    final        1200       1, 1   0.3 0.5         0.7
#> 550          1  55945.125    final        1200       1, 1   0.3 0.5         0.7
#> 551          1  55989.679    final        1200       1, 1   0.3 0.5         0.7
#> 552          1  56075.544    final        1200       1, 1   0.3 0.5         0.7
#> 553          1  56138.844    final        1200       1, 1   0.3 0.5         0.7
#> 554          1  56525.814    final        1200       1, 1   0.3 0.5         0.7
#> 555          1  56553.832    final        1200       1, 1   0.3 0.5         0.7
#> 556          1  56567.575    final        1200       1, 1   0.3 0.5         0.7
#> 557          1  56613.866    final        1200       1, 1   0.3 0.5         0.7
#> 558          1  56865.899    final        1200       1, 1   0.3 0.5         0.7
#> 559          1  56966.219    final        1200       1, 1   0.3 0.5         0.7
#> 560          1  57028.772    final        1200       1, 1   0.3 0.5         0.7
#> 561          1  57045.993    final        1200       1, 1   0.3 0.5         0.7
#> 562          1  57124.053    final        1200       1, 1   0.3 0.5         0.7
#> 563          1  57283.258    final        1200       1, 1   0.3 0.5         0.7
#> 564          1  57285.737    final        1200       1, 1   0.3 0.5         0.7
#> 565          1  57304.281    final        1200       1, 1   0.3 0.5         0.7
#> 566          1  57514.570    final        1200       1, 1   0.3 0.5         0.7
#> 567          1  57704.657    final        1200       1, 1   0.3 0.5         0.7
#> 568          1  57793.150    final        1200       1, 1   0.3 0.5         0.7
#> 569          1  57819.829    final        1200       1, 1   0.3 0.5         0.7
#> 570          1  57820.648    final        1200       1, 1   0.3 0.5         0.7
#> 571          1  58043.399    final        1200       1, 1   0.3 0.5         0.7
#> 572          1  58335.286    final        1200       1, 1   0.3 0.5         0.7
#> 573          1  58365.937    final        1200       1, 1   0.3 0.5         0.7
#> 574          1  58379.242    final        1200       1, 1   0.3 0.5         0.7
#> 575          1  58390.011    final        1200       1, 1   0.3 0.5         0.7
#> 576          1  58406.542    final        1200       1, 1   0.3 0.5         0.7
#> 577          1  58440.588    final        1200       1, 1   0.3 0.5         0.7
#> 578          1  58519.503    final        1200       1, 1   0.3 0.5         0.7
#> 579          1  58540.316    final        1200       1, 1   0.3 0.5         0.7
#> 580          1  58565.771    final        1200       1, 1   0.3 0.5         0.7
#> 581          1  58758.564    final        1200       1, 1   0.3 0.5         0.7
#> 582          1  59143.721    final        1200       1, 1   0.3 0.5         0.7
#> 583          1  59355.617    final        1200       1, 1   0.3 0.5         0.7
#> 584          1  59493.158    final        1200       1, 1   0.3 0.5         0.7
#> 585          1  59505.528    final        1200       1, 1   0.3 0.5         0.7
#> 586          1  59647.357    final        1200       1, 1   0.3 0.5         0.7
#> 587          1  59693.002    final        1200       1, 1   0.3 0.5         0.7
#> 588          1  59844.366    final        1200       1, 1   0.3 0.5         0.7
#> 589          1  60051.141    final        1200       1, 1   0.3 0.5         0.7
#> 590          1  60170.834    final        1200       1, 1   0.3 0.5         0.7
#> 591          1  60259.325    final        1200       1, 1   0.3 0.5         0.7
#> 592          1  60309.027    final        1200       1, 1   0.3 0.5         0.7
#> 593          1  60549.985    final        1200       1, 1   0.3 0.5         0.7
#> 594          1  60568.952    final        1200       1, 1   0.3 0.5         0.7
#> 595          1  60597.351    final        1200       1, 1   0.3 0.5         0.7
#> 596          1  60621.245    final        1200       1, 1   0.3 0.5         0.7
#> 597          1  60680.862    final        1200       1, 1   0.3 0.5         0.7
#> 598          1  60961.538    final        1200       1, 1   0.3 0.5         0.7
#> 599          1  61036.935    final        1200       1, 1   0.3 0.5         0.7
#> 600          1  61134.133    final        1200       1, 1   0.3 0.5         0.7
#> 601          1  61218.277    final        1200       1, 1   0.3 0.5         0.7
#> 602          1  61317.698    final        1200       1, 1   0.3 0.5         0.7
#> 603          1  61692.751    final        1200       1, 1   0.3 0.5         0.7
#> 604          1  61706.673    final        1200       1, 1   0.3 0.5         0.7
#> 605          1  61900.935    final        1200       1, 1   0.3 0.5         0.7
#> 606          1  61939.962    final        1200       1, 1   0.3 0.5         0.7
#> 607          1  61951.216    final        1200       1, 1   0.3 0.5         0.7
#> 608          1  62097.969    final        1200       1, 1   0.3 0.5         0.7
#> 609          1  62139.682    final        1200       1, 1   0.3 0.5         0.7
#> 610          1  62353.886    final        1200       1, 1   0.3 0.5         0.7
#> 611          1  62629.503    final        1200       1, 1   0.3 0.5         0.7
#> 612          1  62630.932    final        1200       1, 1   0.3 0.5         0.7
#> 613          1  62679.548    final        1200       1, 1   0.3 0.5         0.7
#> 614          1  62765.430    final        1200       1, 1   0.3 0.5         0.7
#> 615          1  63081.481    final        1200       1, 1   0.3 0.5         0.7
#> 616          1  63252.240    final        1200       1, 1   0.3 0.5         0.7
#> 617          1  63429.838    final        1200       1, 1   0.3 0.5         0.7
#> 618          1  63478.441    final        1200       1, 1   0.3 0.5         0.7
#> 619          1  63537.895    final        1200       1, 1   0.3 0.5         0.7
#> 620          1  63738.552    final        1200       1, 1   0.3 0.5         0.7
#> 621          1  63766.193    final        1200       1, 1   0.3 0.5         0.7
#> 622          1  63849.852    final        1200       1, 1   0.3 0.5         0.7
#> 623          1  63908.178    final        1200       1, 1   0.3 0.5         0.7
#> 624          1  63958.746    final        1200       1, 1   0.3 0.5         0.7
#> 625          1  64042.874    final        1200       1, 1   0.3 0.5         0.7
#> 626          1  64043.964    final        1200       1, 1   0.3 0.5         0.7
#> 627          1  64092.586    final        1200       1, 1   0.3 0.5         0.7
#> 628          1  64384.392    final        1200       1, 1   0.3 0.5         0.7
#> 629          1  64622.427    final        1200       1, 1   0.3 0.5         0.7
#> 630          1  64726.149    final        1200       1, 1   0.3 0.5         0.7
#> 631          1  64896.634    final        1200       1, 1   0.3 0.5         0.7
#> 632          1  64941.244    final        1200       1, 1   0.3 0.5         0.7
#> 633          1  65104.476    final        1200       1, 1   0.3 0.5         0.7
#> 634          1  65145.380    final        1200       1, 1   0.3 0.5         0.7
#> 635          1  65190.911    final        1200       1, 1   0.3 0.5         0.7
#> 636          1  65206.203    final        1200       1, 1   0.3 0.5         0.7
#> 637          1  65261.575    final        1200       1, 1   0.3 0.5         0.7
#> 638          1  65313.911    final        1200       1, 1   0.3 0.5         0.7
#> 639          1  65480.112    final        1200       1, 1   0.3 0.5         0.7
#> 640          1  65615.082    final        1200       1, 1   0.3 0.5         0.7
#> 641          1  65661.938    final        1200       1, 1   0.3 0.5         0.7
#> 642          1  65662.522    final        1200       1, 1   0.3 0.5         0.7
#> 643          1  65706.474    final        1200       1, 1   0.3 0.5         0.7
#> 644          1  65892.242    final        1200       1, 1   0.3 0.5         0.7
#> 645          1  65974.778    final        1200       1, 1   0.3 0.5         0.7
#> 646          1  66061.141    final        1200       1, 1   0.3 0.5         0.7
#> 647          1  66120.557    final        1200       1, 1   0.3 0.5         0.7
#> 648          1  66155.779    final        1200       1, 1   0.3 0.5         0.7
#> 649          1  66179.610    final        1200       1, 1   0.3 0.5         0.7
#> 650          1  66232.036    final        1200       1, 1   0.3 0.5         0.7
#> 651          1  66251.738    final        1200       1, 1   0.3 0.5         0.7
#> 652          1  66300.402    final        1200       1, 1   0.3 0.5         0.7
#> 653          1  66409.281    final        1200       1, 1   0.3 0.5         0.7
#> 654          1  66500.354    final        1200       1, 1   0.3 0.5         0.7
#> 655          1  66531.357    final        1200       1, 1   0.3 0.5         0.7
#> 656          1  66549.506    final        1200       1, 1   0.3 0.5         0.7
#> 657          1  66607.863    final        1200       1, 1   0.3 0.5         0.7
#> 658          1  66711.053    final        1200       1, 1   0.3 0.5         0.7
#> 659          1  66789.055    final        1200       1, 1   0.3 0.5         0.7
#> 660          1  66864.624    final        1200       1, 1   0.3 0.5         0.7
#> 661          1  66987.837    final        1200       1, 1   0.3 0.5         0.7
#> 662          1  67045.495    final        1200       1, 1   0.3 0.5         0.7
#> 663          1  67085.442    final        1200       1, 1   0.3 0.5         0.7
#> 664          1  67256.160    final        1200       1, 1   0.3 0.5         0.7
#> 665          1  67272.713    final        1200       1, 1   0.3 0.5         0.7
#> 666          1  67274.156    final        1200       1, 1   0.3 0.5         0.7
#> 667          1  67302.717    final        1200       1, 1   0.3 0.5         0.7
#> 668          1  67529.413    final        1200       1, 1   0.3 0.5         0.7
#> 669          1  67531.311    final        1200       1, 1   0.3 0.5         0.7
#> 670          1  67596.524    final        1200       1, 1   0.3 0.5         0.7
#> 671          1  67779.716    final        1200       1, 1   0.3 0.5         0.7
#> 672          1  67782.172    final        1200       1, 1   0.3 0.5         0.7
#> 673          1  67791.783    final        1200       1, 1   0.3 0.5         0.7
#> 674          1  67853.264    final        1200       1, 1   0.3 0.5         0.7
#> 675          1  67870.560    final        1200       1, 1   0.3 0.5         0.7
#> 676          1  67951.335    final        1200       1, 1   0.3 0.5         0.7
#> 677          1  68104.969    final        1200       1, 1   0.3 0.5         0.7
#> 678          1  68239.093    final        1200       1, 1   0.3 0.5         0.7
#> 679          1  68267.837    final        1200       1, 1   0.3 0.5         0.7
#> 680          1  68317.975    final        1200       1, 1   0.3 0.5         0.7
#> 681          1  68365.241    final        1200       1, 1   0.3 0.5         0.7
#> 682          1  68436.908    final        1200       1, 1   0.3 0.5         0.7
#> 683          1  68526.630    final        1200       1, 1   0.3 0.5         0.7
#> 684          1  68593.085    final        1200       1, 1   0.3 0.5         0.7
#> 685          1  68678.206    final        1200       1, 1   0.3 0.5         0.7
#> 686          1  69203.430    final        1200       1, 1   0.3 0.5         0.7
#> 687          1  69300.440    final        1200       1, 1   0.3 0.5         0.7
#> 688          1  69421.129    final        1200       1, 1   0.3 0.5         0.7
#> 689          1  69445.937    final        1200       1, 1   0.3 0.5         0.7
#> 690          1  69615.857    final        1200       1, 1   0.3 0.5         0.7
#> 691          1  69862.915    final        1200       1, 1   0.3 0.5         0.7
#> 692          1  69911.816    final        1200       1, 1   0.3 0.5         0.7
#> 693          1  69926.664    final        1200       1, 1   0.3 0.5         0.7
#> 694          1  69968.406    final        1200       1, 1   0.3 0.5         0.7
#> 695          1  70066.329    final        1200       1, 1   0.3 0.5         0.7
#> 696          1  70141.532    final        1200       1, 1   0.3 0.5         0.7
#> 697          1  70242.060    final        1200       1, 1   0.3 0.5         0.7
#> 698          1  70283.679    final        1200       1, 1   0.3 0.5         0.7
#> 699          1  70471.452    final        1200       1, 1   0.3 0.5         0.7
#> 700          1  70621.840    final        1200       1, 1   0.3 0.5         0.7
#> 701          1  70842.453    final        1200       1, 1   0.3 0.5         0.7
#> 702          1  70847.229    final        1200       1, 1   0.3 0.5         0.7
#> 703          1  70954.739    final        1200       1, 1   0.3 0.5         0.7
#> 704          1  70993.322    final        1200       1, 1   0.3 0.5         0.7
#> 705          1  71216.207    final        1200       1, 1   0.3 0.5         0.7
#> 706          1  71549.422    final        1200       1, 1   0.3 0.5         0.7
#> 707          1  71604.842    final        1200       1, 1   0.3 0.5         0.7
#> 708          1  71610.793    final        1200       1, 1   0.3 0.5         0.7
#> 709          1  71610.969    final        1200       1, 1   0.3 0.5         0.7
#> 710          1  71930.270    final        1200       1, 1   0.3 0.5         0.7
#> 711          1  71986.941    final        1200       1, 1   0.3 0.5         0.7
#> 712          1  72016.981    final        1200       1, 1   0.3 0.5         0.7
#> 713          1  72027.707    final        1200       1, 1   0.3 0.5         0.7
#> 714          1  72237.131    final        1200       1, 1   0.3 0.5         0.7
#> 715          1  72399.348    final        1200       1, 1   0.3 0.5         0.7
#> 716          1  72411.934    final        1200       1, 1   0.3 0.5         0.7
#> 717          1  72597.935    final        1200       1, 1   0.3 0.5         0.7
#> 718          1  72683.866    final        1200       1, 1   0.3 0.5         0.7
#> 719          1  72710.482    final        1200       1, 1   0.3 0.5         0.7
#> 720          1  72786.144    final        1200       1, 1   0.3 0.5         0.7
#> 721          1  72899.158    final        1200       1, 1   0.3 0.5         0.7
#> 722          1  73049.240    final        1200       1, 1   0.3 0.5         0.7
#> 723          1  73056.043    final        1200       1, 1   0.3 0.5         0.7
#> 724          1  73145.138    final        1200       1, 1   0.3 0.5         0.7
#> 725          1  73167.030    final        1200       1, 1   0.3 0.5         0.7
#> 726          1  73349.820    final        1200       1, 1   0.3 0.5         0.7
#> 727          1  73454.442    final        1200       1, 1   0.3 0.5         0.7
#> 728          1  73518.224    final        1200       1, 1   0.3 0.5         0.7
#> 729          1  73675.313    final        1200       1, 1   0.3 0.5         0.7
#> 730          1  73781.804    final        1200       1, 1   0.3 0.5         0.7
#> 731          1  73974.138    final        1200       1, 1   0.3 0.5         0.7
#> 732          1  74027.510    final        1200       1, 1   0.3 0.5         0.7
#> 733          1  74099.947    final        1200       1, 1   0.3 0.5         0.7
#> 734          1  74106.016    final        1200       1, 1   0.3 0.5         0.7
#> 735          1  74185.641    final        1200       1, 1   0.3 0.5         0.7
#> 736          1  74473.979    final        1200       1, 1   0.3 0.5         0.7
#> 737          1  74559.185    final        1200       1, 1   0.3 0.5         0.7
#> 738          1  74789.498    final        1200       1, 1   0.3 0.5         0.7
#> 739          1  74992.125    final        1200       1, 1   0.3 0.5         0.7
#> 740          1  75007.126    final        1200       1, 1   0.3 0.5         0.7
#> 741          1  75015.477    final        1200       1, 1   0.3 0.5         0.7
#> 742          1  75227.327    final        1200       1, 1   0.3 0.5         0.7
#> 743          1  75255.453    final        1200       1, 1   0.3 0.5         0.7
#> 744          1  75258.165    final        1200       1, 1   0.3 0.5         0.7
#> 745          1  75423.625    final        1200       1, 1   0.3 0.5         0.7
#> 746          1  75551.646    final        1200       1, 1   0.3 0.5         0.7
#> 747          1  75579.496    final        1200       1, 1   0.3 0.5         0.7
#> 748          1  76000.296    final        1200       1, 1   0.3 0.5         0.7
#> 749          1  76163.305    final        1200       1, 1   0.3 0.5         0.7
#> 750          1  76278.043    final        1200       1, 1   0.3 0.5         0.7
#> 751          1  76292.819    final        1200       1, 1   0.3 0.5         0.7
#> 752          1  76311.980    final        1200       1, 1   0.3 0.5         0.7
#> 753          1  76330.693    final        1200       1, 1   0.3 0.5         0.7
#> 754          1  76461.071    final        1200       1, 1   0.3 0.5         0.7
#> 755          1  76479.902    final        1200       1, 1   0.3 0.5         0.7
#> 756          1  76503.941    final        1200       1, 1   0.3 0.5         0.7
#> 757          1  76742.629    final        1200       1, 1   0.3 0.5         0.7
#> 758          1  76753.105    final        1200       1, 1   0.3 0.5         0.7
#> 759          1  76756.966    final        1200       1, 1   0.3 0.5         0.7
#> 760          1  76770.393    final        1200       1, 1   0.3 0.5         0.7
#> 761          1  76793.008    final        1200       1, 1   0.3 0.5         0.7
#> 762          1  76947.625    final        1200       1, 1   0.3 0.5         0.7
#> 763          1  77018.630    final        1200       1, 1   0.3 0.5         0.7
#> 764          1  77238.789    final        1200       1, 1   0.3 0.5         0.7
#> 765          1  77396.528    final        1200       1, 1   0.3 0.5         0.7
#> 766          1  77704.362    final        1200       1, 1   0.3 0.5         0.7
#> 767          1  77770.197    final        1200       1, 1   0.3 0.5         0.7
#> 768          1  77870.671    final        1200       1, 1   0.3 0.5         0.7
#> 769          1  77935.400    final        1200       1, 1   0.3 0.5         0.7
#> 770          1  78210.653    final        1200       1, 1   0.3 0.5         0.7
#> 771          1  78257.086    final        1200       1, 1   0.3 0.5         0.7
#> 772          1  78282.631    final        1200       1, 1   0.3 0.5         0.7
#> 773          1  78397.149    final        1200       1, 1   0.3 0.5         0.7
#> 774          1  78427.711    final        1200       1, 1   0.3 0.5         0.7
#> 775          1  78614.937    final        1200       1, 1   0.3 0.5         0.7
#> 776          1  78668.688    final        1200       1, 1   0.3 0.5         0.7
#> 777          1  78817.994    final        1200       1, 1   0.3 0.5         0.7
#> 778          1  78997.262    final        1200       1, 1   0.3 0.5         0.7
#> 779          1  79138.482    final        1200       1, 1   0.3 0.5         0.7
#> 780          1  79279.415    final        1200       1, 1   0.3 0.5         0.7
#> 781          1  79485.731    final        1200       1, 1   0.3 0.5         0.7
#> 782          1  79491.140    final        1200       1, 1   0.3 0.5         0.7
#> 783          1  79591.150    final        1200       1, 1   0.3 0.5         0.7
#> 784          1  79638.770    final        1200       1, 1   0.3 0.5         0.7
#> 785          1  79661.359    final        1200       1, 1   0.3 0.5         0.7
#> 786          1  79829.832    final        1200       1, 1   0.3 0.5         0.7
#> 787          1  79883.513    final        1200       1, 1   0.3 0.5         0.7
#> 788          1  80078.860    final        1200       1, 1   0.3 0.5         0.7
#> 789          1  80121.267    final        1200       1, 1   0.3 0.5         0.7
#> 790          1  80239.443    final        1200       1, 1   0.3 0.5         0.7
#> 791          1  80282.073    final        1200       1, 1   0.3 0.5         0.7
#> 792          1  80655.436    final        1200       1, 1   0.3 0.5         0.7
#> 793          1  80762.118    final        1200       1, 1   0.3 0.5         0.7
#> 794          1  80836.961    final        1200       1, 1   0.3 0.5         0.7
#> 795          1  80874.957    final        1200       1, 1   0.3 0.5         0.7
#> 796          1  80913.703    final        1200       1, 1   0.3 0.5         0.7
#> 797          1  80986.217    final        1200       1, 1   0.3 0.5         0.7
#> 798          1  81021.858    final        1200       1, 1   0.3 0.5         0.7
#> 799          1  81246.601    final        1200       1, 1   0.3 0.5         0.7
#> 800          1  81255.458    final        1200       1, 1   0.3 0.5         0.7
#> 801          1  81306.628    final        1200       1, 1   0.3 0.5         0.7
#> 802          1  81418.790    final        1200       1, 1   0.3 0.5         0.7
#> 803          1  81528.742    final        1200       1, 1   0.3 0.5         0.7
#> 804          1  81549.923    final        1200       1, 1   0.3 0.5         0.7
#> 805          1  81940.714    final        1200       1, 1   0.3 0.5         0.7
#> 806          1  81955.120    final        1200       1, 1   0.3 0.5         0.7
#> 807          1  82101.796    final        1200       1, 1   0.3 0.5         0.7
#> 808          1  82599.821    final        1200       1, 1   0.3 0.5         0.7
#> 809          1  82635.760    final        1200       1, 1   0.3 0.5         0.7
#> 810          1  82980.967    final        1200       1, 1   0.3 0.5         0.7
#> 811          1  83004.201    final        1200       1, 1   0.3 0.5         0.7
#> 812          1  83039.889    final        1200       1, 1   0.3 0.5         0.7
#> 813          1  83425.486    final        1200       1, 1   0.3 0.5         0.7
#> 814          1  83588.937    final        1200       1, 1   0.3 0.5         0.7
#> 815          1  83591.303    final        1200       1, 1   0.3 0.5         0.7
#> 816          1  83758.347    final        1200       1, 1   0.3 0.5         0.7
#> 817          1  83794.264    final        1200       1, 1   0.3 0.5         0.7
#> 818          1  83867.759    final        1200       1, 1   0.3 0.5         0.7
#> 819          1  84035.378    final        1200       1, 1   0.3 0.5         0.7
#> 820          1  84268.144    final        1200       1, 1   0.3 0.5         0.7
#> 821          1  84302.503    final        1200       1, 1   0.3 0.5         0.7
#> 822          1  84327.655    final        1200       1, 1   0.3 0.5         0.7
#> 823          1  84397.802    final        1200       1, 1   0.3 0.5         0.7
#> 824          1  84468.712    final        1200       1, 1   0.3 0.5         0.7
#> 825          1  84582.175    final        1200       1, 1   0.3 0.5         0.7
#> 826          1  84685.948    final        1200       1, 1   0.3 0.5         0.7
#> 827          1  84689.767    final        1200       1, 1   0.3 0.5         0.7
#> 828          1  84713.605    final        1200       1, 1   0.3 0.5         0.7
#> 829          1  84721.573    final        1200       1, 1   0.3 0.5         0.7
#> 830          1  84926.175    final        1200       1, 1   0.3 0.5         0.7
#> 831          1  84950.045    final        1200       1, 1   0.3 0.5         0.7
#> 832          1  84969.403    final        1200       1, 1   0.3 0.5         0.7
#> 833          1  85270.051    final        1200       1, 1   0.3 0.5         0.7
#> 834          1  85311.921    final        1200       1, 1   0.3 0.5         0.7
#> 835          1  85356.314    final        1200       1, 1   0.3 0.5         0.7
#> 836          1  85551.853    final        1200       1, 1   0.3 0.5         0.7
#> 837          1  85615.035    final        1200       1, 1   0.3 0.5         0.7
#> 838          1  85882.407    final        1200       1, 1   0.3 0.5         0.7
#> 839          1  86023.044    final        1200       1, 1   0.3 0.5         0.7
#> 840          1  86075.732    final        1200       1, 1   0.3 0.5         0.7
#> 841          1  86124.843    final        1200       1, 1   0.3 0.5         0.7
#> 842          1  86198.431    final        1200       1, 1   0.3 0.5         0.7
#> 843          1  86223.655    final        1200       1, 1   0.3 0.5         0.7
#> 844          1  86254.809    final        1200       1, 1   0.3 0.5         0.7
#> 845          1  86288.413    final        1200       1, 1   0.3 0.5         0.7
#> 846          1  86397.548    final        1200       1, 1   0.3 0.5         0.7
#> 847          1  86581.357    final        1200       1, 1   0.3 0.5         0.7
#> 848          1  86586.439    final        1200       1, 1   0.3 0.5         0.7
#> 849          1  86823.738    final        1200       1, 1   0.3 0.5         0.7
#> 850          1  86858.752    final        1200       1, 1   0.3 0.5         0.7
#> 851          1  86884.466    final        1200       1, 1   0.3 0.5         0.7
#> 852          1  87042.675    final        1200       1, 1   0.3 0.5         0.7
#> 853          1  87106.100    final        1200       1, 1   0.3 0.5         0.7
#> 854          1  87150.812    final        1200       1, 1   0.3 0.5         0.7
#> 855          1  87213.201    final        1200       1, 1   0.3 0.5         0.7
#> 856          1  87263.153    final        1200       1, 1   0.3 0.5         0.7
#> 857          1  87298.762    final        1200       1, 1   0.3 0.5         0.7
#> 858          1  87343.501    final        1200       1, 1   0.3 0.5         0.7
#> 859          1  87427.267    final        1200       1, 1   0.3 0.5         0.7
#> 860          1  87496.294    final        1200       1, 1   0.3 0.5         0.7
#> 861          1  87740.070    final        1200       1, 1   0.3 0.5         0.7
#> 862          1  87745.807    final        1200       1, 1   0.3 0.5         0.7
#> 863          1  87801.460    final        1200       1, 1   0.3 0.5         0.7
#> 864          1  87824.240    final        1200       1, 1   0.3 0.5         0.7
#> 865          1  87897.916    final        1200       1, 1   0.3 0.5         0.7
#> 866          1  87923.901    final        1200       1, 1   0.3 0.5         0.7
#> 867          1  87951.096    final        1200       1, 1   0.3 0.5         0.7
#> 868          1  88004.800    final        1200       1, 1   0.3 0.5         0.7
#> 869          1  88043.950    final        1200       1, 1   0.3 0.5         0.7
#> 870          1  88172.973    final        1200       1, 1   0.3 0.5         0.7
#> 871          1  88229.416    final        1200       1, 1   0.3 0.5         0.7
#> 872          1  88381.354    final        1200       1, 1   0.3 0.5         0.7
#> 873          1  88385.075    final        1200       1, 1   0.3 0.5         0.7
#> 874          1  88388.975    final        1200       1, 1   0.3 0.5         0.7
#> 875          1  88417.211    final        1200       1, 1   0.3 0.5         0.7
#> 876          1  88546.580    final        1200       1, 1   0.3 0.5         0.7
#> 877          1  88995.384    final        1200       1, 1   0.3 0.5         0.7
#> 878          1  89210.229    final        1200       1, 1   0.3 0.5         0.7
#> 879          1  89237.870    final        1200       1, 1   0.3 0.5         0.7
#> 880          1  89569.505    final        1200       1, 1   0.3 0.5         0.7
#> 881          1  89697.908    final        1200       1, 1   0.3 0.5         0.7
#> 882          1  89740.026    final        1200       1, 1   0.3 0.5         0.7
#> 883          1  89989.268    final        1200       1, 1   0.3 0.5         0.7
#> 884          1  90016.250    final        1200       1, 1   0.3 0.5         0.7
#> 885          1  90070.440    final        1200       1, 1   0.3 0.5         0.7
#> 886          1  90206.014    final        1200       1, 1   0.3 0.5         0.7
#> 887          1  90244.202    final        1200       1, 1   0.3 0.5         0.7
#> 888          1  90259.377    final        1200       1, 1   0.3 0.5         0.7
#> 889          1  90328.579    final        1200       1, 1   0.3 0.5         0.7
#> 890          1  90465.499    final        1200       1, 1   0.3 0.5         0.7
#> 891          1  90912.708    final        1200       1, 1   0.3 0.5         0.7
#> 892          1  90941.565    final        1200       1, 1   0.3 0.5         0.7
#> 893          1  91065.252    final        1200       1, 1   0.3 0.5         0.7
#> 894          1  91150.260    final        1200       1, 1   0.3 0.5         0.7
#> 895          1  91153.674    final        1200       1, 1   0.3 0.5         0.7
#> 896          1  91439.745    final        1200       1, 1   0.3 0.5         0.7
#> 897          1  91573.882    final        1200       1, 1   0.3 0.5         0.7
#> 898          1  91726.079    final        1200       1, 1   0.3 0.5         0.7
#> 899          1  91748.031    final        1200       1, 1   0.3 0.5         0.7
#> 900          1  91757.563    final        1200       1, 1   0.3 0.5         0.7
#> 901          1  91869.645    final        1200       1, 1   0.3 0.5         0.7
#> 902          1  91891.327    final        1200       1, 1   0.3 0.5         0.7
#> 903          1  91907.225    final        1200       1, 1   0.3 0.5         0.7
#> 904          1  92028.979    final        1200       1, 1   0.3 0.5         0.7
#> 905          1  92092.554    final        1200       1, 1   0.3 0.5         0.7
#> 906          1  92168.529    final        1200       1, 1   0.3 0.5         0.7
#> 907          1  92202.058    final        1200       1, 1   0.3 0.5         0.7
#> 908          1  92322.598    final        1200       1, 1   0.3 0.5         0.7
#> 909          1  92338.463    final        1200       1, 1   0.3 0.5         0.7
#> 910          1  92852.540    final        1200       1, 1   0.3 0.5         0.7
#> 911          1  92997.370    final        1200       1, 1   0.3 0.5         0.7
#> 912          1  93152.625    final        1200       1, 1   0.3 0.5         0.7
#> 913          1  93200.881    final        1200       1, 1   0.3 0.5         0.7
#> 914          1  93285.270    final        1200       1, 1   0.3 0.5         0.7
#> 915          1  93289.545    final        1200       1, 1   0.3 0.5         0.7
#> 916          1  93470.072    final        1200       1, 1   0.3 0.5         0.7
#> 917          1  93704.261    final        1200       1, 1   0.3 0.5         0.7
#> 918          1  93744.784    final        1200       1, 1   0.3 0.5         0.7
#> 919          1  94502.942    final        1200       1, 1   0.3 0.5         0.7
#> 920          1  94568.659    final        1200       1, 1   0.3 0.5         0.7
#> 921          1  94765.767    final        1200       1, 1   0.3 0.5         0.7
#> 922          1  94776.774    final        1200       1, 1   0.3 0.5         0.7
#> 923          1  94798.720    final        1200       1, 1   0.3 0.5         0.7
#> 924          1  94918.257    final        1200       1, 1   0.3 0.5         0.7
#> 925          1  95077.490    final        1200       1, 1   0.3 0.5         0.7
#> 926          1  95117.002    final        1200       1, 1   0.3 0.5         0.7
#> 927          1  95187.321    final        1200       1, 1   0.3 0.5         0.7
#> 928          1  95417.783    final        1200       1, 1   0.3 0.5         0.7
#> 929          1  95562.265    final        1200       1, 1   0.3 0.5         0.7
#> 930          1  95655.485    final        1200       1, 1   0.3 0.5         0.7
#> 931          1  95717.316    final        1200       1, 1   0.3 0.5         0.7
#> 932          1  95823.550    final        1200       1, 1   0.3 0.5         0.7
#> 933          1  95839.724    final        1200       1, 1   0.3 0.5         0.7
#> 934          1  95950.162    final        1200       1, 1   0.3 0.5         0.7
#> 935          1  96010.223    final        1200       1, 1   0.3 0.5         0.7
#> 936          1  96069.026    final        1200       1, 1   0.3 0.5         0.7
#> 937          1  96076.565    final        1200       1, 1   0.3 0.5         0.7
#> 938          1  96137.216    final        1200       1, 1   0.3 0.5         0.7
#> 939          1  96272.500    final        1200       1, 1   0.3 0.5         0.7
#> 940          1  96509.099    final        1200       1, 1   0.3 0.5         0.7
#> 941          1  96572.704    final        1200       1, 1   0.3 0.5         0.7
#> 942          1  96665.077    final        1200       1, 1   0.3 0.5         0.7
#> 943          1  96692.971    final        1200       1, 1   0.3 0.5         0.7
#> 944          1  97118.684    final        1200       1, 1   0.3 0.5         0.7
#> 945          1  97227.626    final        1200       1, 1   0.3 0.5         0.7
#> 946          1  97272.131    final        1200       1, 1   0.3 0.5         0.7
#> 947          1  97292.952    final        1200       1, 1   0.3 0.5         0.7
#> 948          1  97428.044    final        1200       1, 1   0.3 0.5         0.7
#> 949          1  97519.610    final        1200       1, 1   0.3 0.5         0.7
#> 950          1  97523.292    final        1200       1, 1   0.3 0.5         0.7
#> 951          1  97568.225    final        1200       1, 1   0.3 0.5         0.7
#> 952          1  97861.557    final        1200       1, 1   0.3 0.5         0.7
#> 953          1  97893.792    final        1200       1, 1   0.3 0.5         0.7
#> 954          1  97986.194    final        1200       1, 1   0.3 0.5         0.7
#> 955          1  98107.121    final        1200       1, 1   0.3 0.5         0.7
#> 956          1  98175.904    final        1200       1, 1   0.3 0.5         0.7
#> 957          1  98231.630    final        1200       1, 1   0.3 0.5         0.7
#> 958          1  98282.323    final        1200       1, 1   0.3 0.5         0.7
#> 959          1  98316.892    final        1200       1, 1   0.3 0.5         0.7
#> 960          1  98401.219    final        1200       1, 1   0.3 0.5         0.7
#> 961          1  98403.299    final        1200       1, 1   0.3 0.5         0.7
#> 962          1  98414.928    final        1200       1, 1   0.3 0.5         0.7
#> 963          1  98510.463    final        1200       1, 1   0.3 0.5         0.7
#> 964          1  98597.832    final        1200       1, 1   0.3 0.5         0.7
#> 965          1  98656.515    final        1200       1, 1   0.3 0.5         0.7
#> 966          1  98740.289    final        1200       1, 1   0.3 0.5         0.7
#> 967          1  98768.493    final        1200       1, 1   0.3 0.5         0.7
#> 968          1  98803.327    final        1200       1, 1   0.3 0.5         0.7
#> 969          1  98902.923    final        1200       1, 1   0.3 0.5         0.7
#> 970          1  98942.162    final        1200       1, 1   0.3 0.5         0.7
#> 971          1  98954.852    final        1200       1, 1   0.3 0.5         0.7
#> 972          1  98961.708    final        1200       1, 1   0.3 0.5         0.7
#> 973          1  98972.720    final        1200       1, 1   0.3 0.5         0.7
#> 974          1  99049.107    final        1200       1, 1   0.3 0.5         0.7
#> 975          1  99140.596    final        1200       1, 1   0.3 0.5         0.7
#> 976          1  99222.340    final        1200       1, 1   0.3 0.5         0.7
#> 977          1  99254.600    final        1200       1, 1   0.3 0.5         0.7
#> 978          1  99283.431    final        1200       1, 1   0.3 0.5         0.7
#> 979          1  99323.913    final        1200       1, 1   0.3 0.5         0.7
#> 980          1  99330.644    final        1200       1, 1   0.3 0.5         0.7
#> 981          1  99342.016    final        1200       1, 1   0.3 0.5         0.7
#> 982          1  99467.212    final        1200       1, 1   0.3 0.5         0.7
#> 983          1  99645.770    final        1200       1, 1   0.3 0.5         0.7
#> 984          1  99740.947    final        1200       1, 1   0.3 0.5         0.7
#> 985          1  99780.664    final        1200       1, 1   0.3 0.5         0.7
#> 986          1  99803.380    final        1200       1, 1   0.3 0.5         0.7
#> 987          1  99828.369    final        1200       1, 1   0.3 0.5         0.7
#> 988          1 100088.262    final        1200       1, 1   0.3 0.5         0.7
#> 989          1 100109.054    final        1200       1, 1   0.3 0.5         0.7
#> 990          1 100222.318    final        1200       1, 1   0.3 0.5         0.7
#> 991          1 100257.060    final        1200       1, 1   0.3 0.5         0.7
#> 992          1 100405.207    final        1200       1, 1   0.3 0.5         0.7
#> 993          1 100447.849    final        1200       1, 1   0.3 0.5         0.7
#> 994          1 100465.243    final        1200       1, 1   0.3 0.5         0.7
#> 995          1 100475.682    final        1200       1, 1   0.3 0.5         0.7
#> 996          1 100523.050    final        1200       1, 1   0.3 0.5         0.7
#> 997          1 100639.930    final        1200       1, 1   0.3 0.5         0.7
#> 998          1 100702.729    final        1200       1, 1   0.3 0.5         0.7
#> 999          1 101119.416    final        1200       1, 1   0.3 0.5         0.7
#> 1000         1 101244.118    final        1200       1, 1   0.3 0.5         0.7
#> 1001         1 101257.593    final        1200       1, 1   0.3 0.5         0.7
#> 1002         1 101275.713    final        1200       1, 1   0.3 0.5         0.7
#> 1003         1 101414.355    final        1200       1, 1   0.3 0.5         0.7
#> 1004         1 101428.385    final        1200       1, 1   0.3 0.5         0.7
#> 1005         1 101484.624    final        1200       1, 1   0.3 0.5         0.7
#> 1006         1 101555.378    final        1200       1, 1   0.3 0.5         0.7
#> 1007         1 101634.961    final        1200       1, 1   0.3 0.5         0.7
#> 1008         1 101981.745    final        1200       1, 1   0.3 0.5         0.7
#> 1009         1 102104.615    final        1200       1, 1   0.3 0.5         0.7
#> 1010         1 102163.920    final        1200       1, 1   0.3 0.5         0.7
#> 1011         1 102327.530    final        1200       1, 1   0.3 0.5         0.7
#> 1012         1 102396.339    final        1200       1, 1   0.3 0.5         0.7
#> 1013         1 102549.657    final        1200       1, 1   0.3 0.5         0.7
#> 1014         1 102562.425    final        1200       1, 1   0.3 0.5         0.7
#> 1015         1 103249.790    final        1200       1, 1   0.3 0.5         0.7
#> 1016         1 103249.965    final        1200       1, 1   0.3 0.5         0.7
#> 1017         1 103347.267    final        1200       1, 1   0.3 0.5         0.7
#> 1018         1 103349.504    final        1200       1, 1   0.3 0.5         0.7
#> 1019         1 103486.401    final        1200       1, 1   0.3 0.5         0.7
#> 1020         1 103567.840    final        1200       1, 1   0.3 0.5         0.7
#> 1021         1 103626.929    final        1200       1, 1   0.3 0.5         0.7
#> 1022         1 103630.736    final        1200       1, 1   0.3 0.5         0.7
#> 1023         1 103733.523    final        1200       1, 1   0.3 0.5         0.7
#> 1024         1 104113.412    final        1200       1, 1   0.3 0.5         0.7
#> 1025         1 104240.377    final        1200       1, 1   0.3 0.5         0.7
#> 1026         1 104319.872    final        1200       1, 1   0.3 0.5         0.7
#> 1027         1 104495.322    final        1200       1, 1   0.3 0.5         0.7
#> 1028         1 104624.711    final        1200       1, 1   0.3 0.5         0.7
#> 1029         1 104662.803    final        1200       1, 1   0.3 0.5         0.7
#> 1030         1 104681.647    final        1200       1, 1   0.3 0.5         0.7
#> 1031         1 104688.168    final        1200       1, 1   0.3 0.5         0.7
#> 1032         1 104782.121    final        1200       1, 1   0.3 0.5         0.7
#> 1033         1 104854.251    final        1200       1, 1   0.3 0.5         0.7
#> 1034         1 104958.743    final        1200       1, 1   0.3 0.5         0.7
#> 1035         1 105115.828    final        1200       1, 1   0.3 0.5         0.7
#> 1036         1 105136.070    final        1200       1, 1   0.3 0.5         0.7
#> 1037         1 105219.371    final        1200       1, 1   0.3 0.5         0.7
#> 1038         1 105258.505    final        1200       1, 1   0.3 0.5         0.7
#> 1039         1 105281.523    final        1200       1, 1   0.3 0.5         0.7
#> 1040         1 105297.000    final        1200       1, 1   0.3 0.5         0.7
#> 1041         1 105418.816    final        1200       1, 1   0.3 0.5         0.7
#> 1042         1 105426.148    final        1200       1, 1   0.3 0.5         0.7
#> 1043         1 105428.825    final        1200       1, 1   0.3 0.5         0.7
#> 1044         1 105448.577    final        1200       1, 1   0.3 0.5         0.7
#> 1045         1 105543.231    final        1200       1, 1   0.3 0.5         0.7
#> 1046         1 105863.364    final        1200       1, 1   0.3 0.5         0.7
#> 1047         1 106023.661    final        1200       1, 1   0.3 0.5         0.7
#> 1048         1 106145.893    final        1200       1, 1   0.3 0.5         0.7
#> 1049         1 106168.866    final        1200       1, 1   0.3 0.5         0.7
#> 1050         1 106349.303    final        1200       1, 1   0.3 0.5         0.7
#> 1051         1 106444.382    final        1200       1, 1   0.3 0.5         0.7
#> 1052         1 106455.206    final        1200       1, 1   0.3 0.5         0.7
#> 1053         1 106469.702    final        1200       1, 1   0.3 0.5         0.7
#> 1054         1 106633.244    final        1200       1, 1   0.3 0.5         0.7
#> 1055         1 106733.778    final        1200       1, 1   0.3 0.5         0.7
#> 1056         1 106752.282    final        1200       1, 1   0.3 0.5         0.7
#> 1057         1 106759.758    final        1200       1, 1   0.3 0.5         0.7
#> 1058         1 106848.978    final        1200       1, 1   0.3 0.5         0.7
#> 1059         1 107175.665    final        1200       1, 1   0.3 0.5         0.7
#> 1060         1 107321.810    final        1200       1, 1   0.3 0.5         0.7
#> 1061         1 107335.037    final        1200       1, 1   0.3 0.5         0.7
#> 1062         1 107525.639    final        1200       1, 1   0.3 0.5         0.7
#> 1063         1 107534.440    final        1200       1, 1   0.3 0.5         0.7
#> 1064         1 107600.539    final        1200       1, 1   0.3 0.5         0.7
#> 1065         1 107634.871    final        1200       1, 1   0.3 0.5         0.7
#> 1066         1 107714.418    final        1200       1, 1   0.3 0.5         0.7
#> 1067         1 107845.425    final        1200       1, 1   0.3 0.5         0.7
#> 1068         1 107858.359    final        1200       1, 1   0.3 0.5         0.7
#> 1069         1 107936.616    final        1200       1, 1   0.3 0.5         0.7
#> 1070         1 108022.680    final        1200       1, 1   0.3 0.5         0.7
#> 1071         1 108198.371    final        1200       1, 1   0.3 0.5         0.7
#> 1072         1 108207.406    final        1200       1, 1   0.3 0.5         0.7
#> 1073         1 108248.969    final        1200       1, 1   0.3 0.5         0.7
#> 1074         1 108250.240    final        1200       1, 1   0.3 0.5         0.7
#> 1075         1 108266.313    final        1200       1, 1   0.3 0.5         0.7
#> 1076         1 108504.079    final        1200       1, 1   0.3 0.5         0.7
#> 1077         1 108614.615    final        1200       1, 1   0.3 0.5         0.7
#> 1078         1 108741.302    final        1200       1, 1   0.3 0.5         0.7
#> 1079         1 108761.152    final        1200       1, 1   0.3 0.5         0.7
#> 1080         1 108772.285    final        1200       1, 1   0.3 0.5         0.7
#> 1081         1 108782.409    final        1200       1, 1   0.3 0.5         0.7
#> 1082         1 108830.063    final        1200       1, 1   0.3 0.5         0.7
#> 1083         1 108922.839    final        1200       1, 1   0.3 0.5         0.7
#> 1084         1 108942.610    final        1200       1, 1   0.3 0.5         0.7
#> 1085         1 109046.777    final        1200       1, 1   0.3 0.5         0.7
#> 1086         1 109302.850    final        1200       1, 1   0.3 0.5         0.7
#> 1087         1 109432.811    final        1200       1, 1   0.3 0.5         0.7
#> 1088         1 109567.446    final        1200       1, 1   0.3 0.5         0.7
#> 1089         1 109578.014    final        1200       1, 1   0.3 0.5         0.7
#> 1090         1 109711.362    final        1200       1, 1   0.3 0.5         0.7
#> 1091         1 109851.967    final        1200       1, 1   0.3 0.5         0.7
#> 1092         1 109868.008    final        1200       1, 1   0.3 0.5         0.7
#> 1093         1 110050.822    final        1200       1, 1   0.3 0.5         0.7
#> 1094         1 110088.093    final        1200       1, 1   0.3 0.5         0.7
#> 1095         1 110356.185    final        1200       1, 1   0.3 0.5         0.7
#> 1096         1 110369.063    final        1200       1, 1   0.3 0.5         0.7
#> 1097         1 110492.280    final        1200       1, 1   0.3 0.5         0.7
#> 1098         1 110542.976    final        1200       1, 1   0.3 0.5         0.7
#> 1099         1 110639.465    final        1200       1, 1   0.3 0.5         0.7
#> 1100         1 110992.856    final        1200       1, 1   0.3 0.5         0.7
#> 1101         1 111056.341    final        1200       1, 1   0.3 0.5         0.7
#> 1102         1 111060.636    final        1200       1, 1   0.3 0.5         0.7
#> 1103         1 111140.702    final        1200       1, 1   0.3 0.5         0.7
#> 1104         1 111165.092    final        1200       1, 1   0.3 0.5         0.7
#> 1105         1 111273.791    final        1200       1, 1   0.3 0.5         0.7
#> 1106         1 111670.069    final        1200       1, 1   0.3 0.5         0.7
#> 1107         1 111692.574    final        1200       1, 1   0.3 0.5         0.7
#> 1108         1 111718.583    final        1200       1, 1   0.3 0.5         0.7
#> 1109         1 111783.212    final        1200       1, 1   0.3 0.5         0.7
#> 1110         1 111850.892    final        1200       1, 1   0.3 0.5         0.7
#> 1111         1 111929.475    final        1200       1, 1   0.3 0.5         0.7
#> 1112         1 111940.849    final        1200       1, 1   0.3 0.5         0.7
#> 1113         1 111982.372    final        1200       1, 1   0.3 0.5         0.7
#> 1114         1 112144.886    final        1200       1, 1   0.3 0.5         0.7
#> 1115         1 112360.838    final        1200       1, 1   0.3 0.5         0.7
#> 1116         1 112789.219    final        1200       1, 1   0.3 0.5         0.7
#> 1117         1 112802.193    final        1200       1, 1   0.3 0.5         0.7
#> 1118         1 112823.277    final        1200       1, 1   0.3 0.5         0.7
#> 1119         1 112893.708    final        1200       1, 1   0.3 0.5         0.7
#> 1120         1 112991.903    final        1200       1, 1   0.3 0.5         0.7
#> 1121         1 113055.231    final        1200       1, 1   0.3 0.5         0.7
#> 1122         1 113142.806    final        1200       1, 1   0.3 0.5         0.7
#> 1123         1 113418.943    final        1200       1, 1   0.3 0.5         0.7
#> 1124         1 113521.296    final        1200       1, 1   0.3 0.5         0.7
#> 1125         1 113671.805    final        1200       1, 1   0.3 0.5         0.7
#> 1126         1 113703.136    final        1200       1, 1   0.3 0.5         0.7
#> 1127         1 113764.584    final        1200       1, 1   0.3 0.5         0.7
#> 1128         1 113804.885    final        1200       1, 1   0.3 0.5         0.7
#> 1129         1 113950.178    final        1200       1, 1   0.3 0.5         0.7
#> 1130         1 114082.012    final        1200       1, 1   0.3 0.5         0.7
#> 1131         1 114087.477    final        1200       1, 1   0.3 0.5         0.7
#> 1132         1 114143.113    final        1200       1, 1   0.3 0.5         0.7
#> 1133         1 114559.314    final        1200       1, 1   0.3 0.5         0.7
#> 1134         1 114651.706    final        1200       1, 1   0.3 0.5         0.7
#> 1135         1 114694.807    final        1200       1, 1   0.3 0.5         0.7
#> 1136         1 114716.798    final        1200       1, 1   0.3 0.5         0.7
#> 1137         1 114740.769    final        1200       1, 1   0.3 0.5         0.7
#> 1138         1 114783.765    final        1200       1, 1   0.3 0.5         0.7
#> 1139         1 114824.705    final        1200       1, 1   0.3 0.5         0.7
#> 1140         1 114913.581    final        1200       1, 1   0.3 0.5         0.7
#> 1141         1 115284.625    final        1200       1, 1   0.3 0.5         0.7
#> 1142         1 115438.424    final        1200       1, 1   0.3 0.5         0.7
#> 1143         1 115803.184    final        1200       1, 1   0.3 0.5         0.7
#> 1144         1 115818.691    final        1200       1, 1   0.3 0.5         0.7
#> 1145         1 115941.549    final        1200       1, 1   0.3 0.5         0.7
#> 1146         1 116024.034    final        1200       1, 1   0.3 0.5         0.7
#> 1147         1 116072.440    final        1200       1, 1   0.3 0.5         0.7
#> 1148         1 116095.483    final        1200       1, 1   0.3 0.5         0.7
#> 1149         1 116238.238    final        1200       1, 1   0.3 0.5         0.7
#> 1150         1 116252.925    final        1200       1, 1   0.3 0.5         0.7
#> 1151         1 116256.080    final        1200       1, 1   0.3 0.5         0.7
#> 1152         1 116365.082    final        1200       1, 1   0.3 0.5         0.7
#> 1153         1 116500.789    final        1200       1, 1   0.3 0.5         0.7
#> 1154         1 116588.887    final        1200       1, 1   0.3 0.5         0.7
#> 1155         1 116598.745    final        1200       1, 1   0.3 0.5         0.7
#> 1156         1 116703.742    final        1200       1, 1   0.3 0.5         0.7
#> 1157         1 116791.997    final        1200       1, 1   0.3 0.5         0.7
#> 1158         1 116804.294    final        1200       1, 1   0.3 0.5         0.7
#> 1159         1 116927.411    final        1200       1, 1   0.3 0.5         0.7
#> 1160         1 117013.022    final        1200       1, 1   0.3 0.5         0.7
#> 1161         1 117123.695    final        1200       1, 1   0.3 0.5         0.7
#> 1162         1 117253.296    final        1200       1, 1   0.3 0.5         0.7
#> 1163         1 117406.701    final        1200       1, 1   0.3 0.5         0.7
#> 1164         1 117499.543    final        1200       1, 1   0.3 0.5         0.7
#> 1165         1 117503.317    final        1200       1, 1   0.3 0.5         0.7
#> 1166         1 117529.537    final        1200       1, 1   0.3 0.5         0.7
#> 1167         1 117569.013    final        1200       1, 1   0.3 0.5         0.7
#> 1168         1 117732.509    final        1200       1, 1   0.3 0.5         0.7
#> 1169         1 117969.790    final        1200       1, 1   0.3 0.5         0.7
#> 1170         1 117996.273    final        1200       1, 1   0.3 0.5         0.7
#> 1171         1 118014.707    final        1200       1, 1   0.3 0.5         0.7
#> 1172         1 118056.075    final        1200       1, 1   0.3 0.5         0.7
#> 1173         1 118242.236    final        1200       1, 1   0.3 0.5         0.7
#> 1174         1 118395.299    final        1200       1, 1   0.3 0.5         0.7
#> 1175         1 118409.497    final        1200       1, 1   0.3 0.5         0.7
#> 1176         1 118440.965    final        1200       1, 1   0.3 0.5         0.7
#> 1177         1 118649.276    final        1200       1, 1   0.3 0.5         0.7
#> 1178         1 118798.017    final        1200       1, 1   0.3 0.5         0.7
#> 1179         1 118998.542    final        1200       1, 1   0.3 0.5         0.7
#> 1180         1 119000.236    final        1200       1, 1   0.3 0.5         0.7
#> 1181         1 119108.205    final        1200       1, 1   0.3 0.5         0.7
#> 1182         1 119591.057    final        1200       1, 1   0.3 0.5         0.7
#> 1183         1 119608.843    final        1200       1, 1   0.3 0.5         0.7
#> 1184         1 119631.644    final        1200       1, 1   0.3 0.5         0.7
#> 1185         1 119668.178    final        1200       1, 1   0.3 0.5         0.7
#> 1186         1 119782.491    final        1200       1, 1   0.3 0.5         0.7
#> 1187         1 119850.749    final        1200       1, 1   0.3 0.5         0.7
#> 1188         1 119851.459    final        1200       1, 1   0.3 0.5         0.7
#> 1189         1 120050.378    final        1200       1, 1   0.3 0.5         0.7
#> 1190         1 120364.972    final        1200       1, 1   0.3 0.5         0.7
#> 1191         1 120737.366    final        1200       1, 1   0.3 0.5         0.7
#> 1192         1 120898.854    final        1200       1, 1   0.3 0.5         0.7
#> 1193         1 120985.285    final        1200       1, 1   0.3 0.5         0.7
#> 1194         2   1199.569    final        1200       1, 1   0.3 0.5         0.7
#> 1195         2   1301.395    final        1200       1, 1   0.3 0.5         0.7
#> 1196         2   1306.518    final        1200       1, 1   0.3 0.5         0.7
#> 1197         2   1376.316    final        1200       1, 1   0.3 0.5         0.7
#> 1198         2   1540.211    final        1200       1, 1   0.3 0.5         0.7
#> 1199         2   1582.420    final        1200       1, 1   0.3 0.5         0.7
#> 1200         2   1618.329    final        1200       1, 1   0.3 0.5         0.7
#> 1201         2   1642.896    final        1200       1, 1   0.3 0.5         0.7
#> 1202         2   1771.869    final        1200       1, 1   0.3 0.5         0.7
#> 1203         2   1920.191    final        1200       1, 1   0.3 0.5         0.7
#> 1204         2   1922.441    final        1200       1, 1   0.3 0.5         0.7
#> 1205         2   1926.363    final        1200       1, 1   0.3 0.5         0.7
#> 1206         2   1955.521    final        1200       1, 1   0.3 0.5         0.7
#> 1207         2   1994.133    final        1200       1, 1   0.3 0.5         0.7
#> 1208         2   2095.205    final        1200       1, 1   0.3 0.5         0.7
#> 1209         2   2231.999    final        1200       1, 1   0.3 0.5         0.7
#> 1210         2   2256.446    final        1200       1, 1   0.3 0.5         0.7
#> 1211         2   2358.267    final        1200       1, 1   0.3 0.5         0.7
#> 1212         2   2487.977    final        1200       1, 1   0.3 0.5         0.7
#> 1213         2   2498.284    final        1200       1, 1   0.3 0.5         0.7
#> 1214         2   2535.472    final        1200       1, 1   0.3 0.5         0.7
#> 1215         2   2568.346    final        1200       1, 1   0.3 0.5         0.7
#> 1216         2   2714.095    final        1200       1, 1   0.3 0.5         0.7
#> 1217         2   2843.770    final        1200       1, 1   0.3 0.5         0.7
#> 1218         2   2938.614    final        1200       1, 1   0.3 0.5         0.7
#> 1219         2   3022.327    final        1200       1, 1   0.3 0.5         0.7
#> 1220         2   3052.053    final        1200       1, 1   0.3 0.5         0.7
#> 1221         2   3198.202    final        1200       1, 1   0.3 0.5         0.7
#> 1222         2   3203.908    final        1200       1, 1   0.3 0.5         0.7
#> 1223         2   3208.932    final        1200       1, 1   0.3 0.5         0.7
#> 1224         2   3381.639    final        1200       1, 1   0.3 0.5         0.7
#> 1225         2   3437.564    final        1200       1, 1   0.3 0.5         0.7
#> 1226         2   3646.787    final        1200       1, 1   0.3 0.5         0.7
#> 1227         2   3653.872    final        1200       1, 1   0.3 0.5         0.7
#> 1228         2   3823.012    final        1200       1, 1   0.3 0.5         0.7
#> 1229         2   3886.451    final        1200       1, 1   0.3 0.5         0.7
#> 1230         2   4035.940    final        1200       1, 1   0.3 0.5         0.7
#> 1231         2   4115.949    final        1200       1, 1   0.3 0.5         0.7
#> 1232         2   4181.821    final        1200       1, 1   0.3 0.5         0.7
#> 1233         2   4453.577    final        1200       1, 1   0.3 0.5         0.7
#> 1234         2   4471.643    final        1200       1, 1   0.3 0.5         0.7
#> 1235         2   4479.306    final        1200       1, 1   0.3 0.5         0.7
#> 1236         2   4674.893    final        1200       1, 1   0.3 0.5         0.7
#> 1237         2   4770.809    final        1200       1, 1   0.3 0.5         0.7
#> 1238         2   5147.128    final        1200       1, 1   0.3 0.5         0.7
#> 1239         2   5278.726    final        1200       1, 1   0.3 0.5         0.7
#> 1240         2   5640.505    final        1200       1, 1   0.3 0.5         0.7
#> 1241         2   5771.841    final        1200       1, 1   0.3 0.5         0.7
#> 1242         2   5822.476    final        1200       1, 1   0.3 0.5         0.7
#> 1243         2   5854.723    final        1200       1, 1   0.3 0.5         0.7
#> 1244         2   5885.984    final        1200       1, 1   0.3 0.5         0.7
#> 1245         2   5948.282    final        1200       1, 1   0.3 0.5         0.7
#> 1246         2   6225.918    final        1200       1, 1   0.3 0.5         0.7
#> 1247         2   6247.964    final        1200       1, 1   0.3 0.5         0.7
#> 1248         2   6425.612    final        1200       1, 1   0.3 0.5         0.7
#> 1249         2   6490.655    final        1200       1, 1   0.3 0.5         0.7
#> 1250         2   6565.277    final        1200       1, 1   0.3 0.5         0.7
#> 1251         2   6590.207    final        1200       1, 1   0.3 0.5         0.7
#> 1252         2   6944.960    final        1200       1, 1   0.3 0.5         0.7
#> 1253         2   7031.020    final        1200       1, 1   0.3 0.5         0.7
#> 1254         2   7173.709    final        1200       1, 1   0.3 0.5         0.7
#> 1255         2   7218.016    final        1200       1, 1   0.3 0.5         0.7
#> 1256         2   7372.911    final        1200       1, 1   0.3 0.5         0.7
#> 1257         2   7374.209    final        1200       1, 1   0.3 0.5         0.7
#> 1258         2   7653.211    final        1200       1, 1   0.3 0.5         0.7
#> 1259         2   7717.314    final        1200       1, 1   0.3 0.5         0.7
#> 1260         2   7800.214    final        1200       1, 1   0.3 0.5         0.7
#> 1261         2   8007.810    final        1200       1, 1   0.3 0.5         0.7
#> 1262         2   8232.046    final        1200       1, 1   0.3 0.5         0.7
#> 1263         2   8266.772    final        1200       1, 1   0.3 0.5         0.7
#> 1264         2   8492.093    final        1200       1, 1   0.3 0.5         0.7
#> 1265         2   8516.065    final        1200       1, 1   0.3 0.5         0.7
#> 1266         2   8586.045    final        1200       1, 1   0.3 0.5         0.7
#> 1267         2   8617.734    final        1200       1, 1   0.3 0.5         0.7
#> 1268         2   8823.172    final        1200       1, 1   0.3 0.5         0.7
#> 1269         2   8954.236    final        1200       1, 1   0.3 0.5         0.7
#> 1270         2   9061.472    final        1200       1, 1   0.3 0.5         0.7
#> 1271         2   9102.075    final        1200       1, 1   0.3 0.5         0.7
#> 1272         2   9149.096    final        1200       1, 1   0.3 0.5         0.7
#> 1273         2   9274.947    final        1200       1, 1   0.3 0.5         0.7
#> 1274         2   9364.492    final        1200       1, 1   0.3 0.5         0.7
#> 1275         2   9875.398    final        1200       1, 1   0.3 0.5         0.7
#> 1276         2  10039.525    final        1200       1, 1   0.3 0.5         0.7
#> 1277         2  10352.278    final        1200       1, 1   0.3 0.5         0.7
#> 1278         2  10366.698    final        1200       1, 1   0.3 0.5         0.7
#> 1279         2  10509.301    final        1200       1, 1   0.3 0.5         0.7
#> 1280         2  10685.437    final        1200       1, 1   0.3 0.5         0.7
#> 1281         2  10716.271    final        1200       1, 1   0.3 0.5         0.7
#> 1282         2  10919.586    final        1200       1, 1   0.3 0.5         0.7
#> 1283         2  10997.527    final        1200       1, 1   0.3 0.5         0.7
#> 1284         2  11020.445    final        1200       1, 1   0.3 0.5         0.7
#> 1285         2  11314.862    final        1200       1, 1   0.3 0.5         0.7
#> 1286         2  11389.762    final        1200       1, 1   0.3 0.5         0.7
#> 1287         2  11405.231    final        1200       1, 1   0.3 0.5         0.7
#> 1288         2  11412.312    final        1200       1, 1   0.3 0.5         0.7
#> 1289         2  11446.720    final        1200       1, 1   0.3 0.5         0.7
#> 1290         2  11918.476    final        1200       1, 1   0.3 0.5         0.7
#> 1291         2  11951.639    final        1200       1, 1   0.3 0.5         0.7
#> 1292         2  11965.314    final        1200       1, 1   0.3 0.5         0.7
#> 1293         2  12033.816    final        1200       1, 1   0.3 0.5         0.7
#> 1294         2  12034.829    final        1200       1, 1   0.3 0.5         0.7
#> 1295         2  12137.425    final        1200       1, 1   0.3 0.5         0.7
#> 1296         2  12203.282    final        1200       1, 1   0.3 0.5         0.7
#> 1297         2  12242.081    final        1200       1, 1   0.3 0.5         0.7
#> 1298         2  12246.876    final        1200       1, 1   0.3 0.5         0.7
#> 1299         2  12314.813    final        1200       1, 1   0.3 0.5         0.7
#> 1300         2  12595.768    final        1200       1, 1   0.3 0.5         0.7
#> 1301         2  12683.811    final        1200       1, 1   0.3 0.5         0.7
#> 1302         2  12685.741    final        1200       1, 1   0.3 0.5         0.7
#> 1303         2  12722.997    final        1200       1, 1   0.3 0.5         0.7
#> 1304         2  12787.413    final        1200       1, 1   0.3 0.5         0.7
#> 1305         2  12986.442    final        1200       1, 1   0.3 0.5         0.7
#> 1306         2  12997.462    final        1200       1, 1   0.3 0.5         0.7
#> 1307         2  13185.706    final        1200       1, 1   0.3 0.5         0.7
#> 1308         2  13351.279    final        1200       1, 1   0.3 0.5         0.7
#> 1309         2  13389.770    final        1200       1, 1   0.3 0.5         0.7
#> 1310         2  13445.905    final        1200       1, 1   0.3 0.5         0.7
#> 1311         2  13446.773    final        1200       1, 1   0.3 0.5         0.7
#> 1312         2  13492.621    final        1200       1, 1   0.3 0.5         0.7
#> 1313         2  13697.265    final        1200       1, 1   0.3 0.5         0.7
#> 1314         2  13708.078    final        1200       1, 1   0.3 0.5         0.7
#> 1315         2  13762.130    final        1200       1, 1   0.3 0.5         0.7
#> 1316         2  14226.699    final        1200       1, 1   0.3 0.5         0.7
#> 1317         2  14278.804    final        1200       1, 1   0.3 0.5         0.7
#> 1318         2  14343.145    final        1200       1, 1   0.3 0.5         0.7
#> 1319         2  14371.831    final        1200       1, 1   0.3 0.5         0.7
#> 1320         2  14609.640    final        1200       1, 1   0.3 0.5         0.7
#> 1321         2  14620.041    final        1200       1, 1   0.3 0.5         0.7
#> 1322         2  14752.971    final        1200       1, 1   0.3 0.5         0.7
#> 1323         2  15109.338    final        1200       1, 1   0.3 0.5         0.7
#> 1324         2  15287.137    final        1200       1, 1   0.3 0.5         0.7
#> 1325         2  15370.965    final        1200       1, 1   0.3 0.5         0.7
#> 1326         2  15425.237    final        1200       1, 1   0.3 0.5         0.7
#> 1327         2  15514.928    final        1200       1, 1   0.3 0.5         0.7
#> 1328         2  15603.119    final        1200       1, 1   0.3 0.5         0.7
#> 1329         2  15887.368    final        1200       1, 1   0.3 0.5         0.7
#> 1330         2  15990.015    final        1200       1, 1   0.3 0.5         0.7
#> 1331         2  16199.285    final        1200       1, 1   0.3 0.5         0.7
#> 1332         2  16458.896    final        1200       1, 1   0.3 0.5         0.7
#> 1333         2  16749.555    final        1200       1, 1   0.3 0.5         0.7
#> 1334         2  16759.765    final        1200       1, 1   0.3 0.5         0.7
#> 1335         2  17162.912    final        1200       1, 1   0.3 0.5         0.7
#> 1336         2  17168.097    final        1200       1, 1   0.3 0.5         0.7
#> 1337         2  17197.871    final        1200       1, 1   0.3 0.5         0.7
#> 1338         2  17221.643    final        1200       1, 1   0.3 0.5         0.7
#> 1339         2  17332.743    final        1200       1, 1   0.3 0.5         0.7
#> 1340         2  17452.003    final        1200       1, 1   0.3 0.5         0.7
#> 1341         2  17528.889    final        1200       1, 1   0.3 0.5         0.7
#> 1342         2  17658.878    final        1200       1, 1   0.3 0.5         0.7
#> 1343         2  17837.149    final        1200       1, 1   0.3 0.5         0.7
#> 1344         2  17877.102    final        1200       1, 1   0.3 0.5         0.7
#> 1345         2  18075.852    final        1200       1, 1   0.3 0.5         0.7
#> 1346         2  18246.528    final        1200       1, 1   0.3 0.5         0.7
#> 1347         2  18308.608    final        1200       1, 1   0.3 0.5         0.7
#> 1348         2  18316.662    final        1200       1, 1   0.3 0.5         0.7
#> 1349         2  18519.330    final        1200       1, 1   0.3 0.5         0.7
#> 1350         2  18576.849    final        1200       1, 1   0.3 0.5         0.7
#> 1351         2  18613.536    final        1200       1, 1   0.3 0.5         0.7
#> 1352         2  18666.470    final        1200       1, 1   0.3 0.5         0.7
#> 1353         2  18908.706    final        1200       1, 1   0.3 0.5         0.7
#> 1354         2  19034.935    final        1200       1, 1   0.3 0.5         0.7
#> 1355         2  19065.305    final        1200       1, 1   0.3 0.5         0.7
#> 1356         2  19165.782    final        1200       1, 1   0.3 0.5         0.7
#> 1357         2  19447.793    final        1200       1, 1   0.3 0.5         0.7
#> 1358         2  19634.354    final        1200       1, 1   0.3 0.5         0.7
#> 1359         2  19751.968    final        1200       1, 1   0.3 0.5         0.7
#> 1360         2  19912.351    final        1200       1, 1   0.3 0.5         0.7
#> 1361         2  19948.947    final        1200       1, 1   0.3 0.5         0.7
#> 1362         2  20019.810    final        1200       1, 1   0.3 0.5         0.7
#> 1363         2  20264.468    final        1200       1, 1   0.3 0.5         0.7
#> 1364         2  20287.554    final        1200       1, 1   0.3 0.5         0.7
#> 1365         2  20353.018    final        1200       1, 1   0.3 0.5         0.7
#> 1366         2  20519.748    final        1200       1, 1   0.3 0.5         0.7
#> 1367         2  20685.690    final        1200       1, 1   0.3 0.5         0.7
#> 1368         2  20774.406    final        1200       1, 1   0.3 0.5         0.7
#> 1369         2  20923.994    final        1200       1, 1   0.3 0.5         0.7
#> 1370         2  20952.402    final        1200       1, 1   0.3 0.5         0.7
#> 1371         2  21078.739    final        1200       1, 1   0.3 0.5         0.7
#> 1372         2  21166.631    final        1200       1, 1   0.3 0.5         0.7
#> 1373         2  21226.616    final        1200       1, 1   0.3 0.5         0.7
#> 1374         2  21438.213    final        1200       1, 1   0.3 0.5         0.7
#> 1375         2  21533.234    final        1200       1, 1   0.3 0.5         0.7
#> 1376         2  21719.170    final        1200       1, 1   0.3 0.5         0.7
#> 1377         2  21776.908    final        1200       1, 1   0.3 0.5         0.7
#> 1378         2  21838.779    final        1200       1, 1   0.3 0.5         0.7
#> 1379         2  21884.860    final        1200       1, 1   0.3 0.5         0.7
#> 1380         2  21917.473    final        1200       1, 1   0.3 0.5         0.7
#> 1381         2  22242.174    final        1200       1, 1   0.3 0.5         0.7
#> 1382         2  22278.526    final        1200       1, 1   0.3 0.5         0.7
#> 1383         2  22309.735    final        1200       1, 1   0.3 0.5         0.7
#> 1384         2  22320.910    final        1200       1, 1   0.3 0.5         0.7
#> 1385         2  22366.943    final        1200       1, 1   0.3 0.5         0.7
#> 1386         2  22446.933    final        1200       1, 1   0.3 0.5         0.7
#> 1387         2  22468.045    final        1200       1, 1   0.3 0.5         0.7
#> 1388         2  22514.096    final        1200       1, 1   0.3 0.5         0.7
#> 1389         2  22529.533    final        1200       1, 1   0.3 0.5         0.7
#> 1390         2  22574.266    final        1200       1, 1   0.3 0.5         0.7
#> 1391         2  22583.500    final        1200       1, 1   0.3 0.5         0.7
#> 1392         2  22866.750    final        1200       1, 1   0.3 0.5         0.7
#> 1393         2  22981.652    final        1200       1, 1   0.3 0.5         0.7
#> 1394         2  23026.535    final        1200       1, 1   0.3 0.5         0.7
#> 1395         2  23164.863    final        1200       1, 1   0.3 0.5         0.7
#> 1396         2  23188.357    final        1200       1, 1   0.3 0.5         0.7
#> 1397         2  23406.841    final        1200       1, 1   0.3 0.5         0.7
#> 1398         2  23489.426    final        1200       1, 1   0.3 0.5         0.7
#> 1399         2  23567.309    final        1200       1, 1   0.3 0.5         0.7
#> 1400         2  23867.946    final        1200       1, 1   0.3 0.5         0.7
#> 1401         2  23997.809    final        1200       1, 1   0.3 0.5         0.7
#> 1402         2  24003.541    final        1200       1, 1   0.3 0.5         0.7
#> 1403         2  24031.886    final        1200       1, 1   0.3 0.5         0.7
#> 1404         2  24088.532    final        1200       1, 1   0.3 0.5         0.7
#> 1405         2  24283.169    final        1200       1, 1   0.3 0.5         0.7
#> 1406         2  24311.750    final        1200       1, 1   0.3 0.5         0.7
#> 1407         2  24346.756    final        1200       1, 1   0.3 0.5         0.7
#> 1408         2  24349.378    final        1200       1, 1   0.3 0.5         0.7
#> 1409         2  24378.660    final        1200       1, 1   0.3 0.5         0.7
#> 1410         2  24388.383    final        1200       1, 1   0.3 0.5         0.7
#> 1411         2  24510.122    final        1200       1, 1   0.3 0.5         0.7
#> 1412         2  24579.629    final        1200       1, 1   0.3 0.5         0.7
#> 1413         2  24707.526    final        1200       1, 1   0.3 0.5         0.7
#> 1414         2  24772.555    final        1200       1, 1   0.3 0.5         0.7
#> 1415         2  24900.579    final        1200       1, 1   0.3 0.5         0.7
#> 1416         2  25011.681    final        1200       1, 1   0.3 0.5         0.7
#> 1417         2  25281.269    final        1200       1, 1   0.3 0.5         0.7
#> 1418         2  25297.193    final        1200       1, 1   0.3 0.5         0.7
#> 1419         2  25317.425    final        1200       1, 1   0.3 0.5         0.7
#> 1420         2  25435.166    final        1200       1, 1   0.3 0.5         0.7
#> 1421         2  25604.746    final        1200       1, 1   0.3 0.5         0.7
#> 1422         2  25689.895    final        1200       1, 1   0.3 0.5         0.7
#> 1423         2  25826.354    final        1200       1, 1   0.3 0.5         0.7
#> 1424         2  26072.755    final        1200       1, 1   0.3 0.5         0.7
#> 1425         2  26266.540    final        1200       1, 1   0.3 0.5         0.7
#> 1426         2  26370.782    final        1200       1, 1   0.3 0.5         0.7
#> 1427         2  26541.025    final        1200       1, 1   0.3 0.5         0.7
#> 1428         2  26576.631    final        1200       1, 1   0.3 0.5         0.7
#> 1429         2  27064.669    final        1200       1, 1   0.3 0.5         0.7
#> 1430         2  27225.308    final        1200       1, 1   0.3 0.5         0.7
#> 1431         2  27263.986    final        1200       1, 1   0.3 0.5         0.7
#> 1432         2  27327.569    final        1200       1, 1   0.3 0.5         0.7
#> 1433         2  27337.162    final        1200       1, 1   0.3 0.5         0.7
#> 1434         2  27380.983    final        1200       1, 1   0.3 0.5         0.7
#> 1435         2  27408.371    final        1200       1, 1   0.3 0.5         0.7
#> 1436         2  27473.770    final        1200       1, 1   0.3 0.5         0.7
#> 1437         2  27512.965    final        1200       1, 1   0.3 0.5         0.7
#> 1438         2  27642.903    final        1200       1, 1   0.3 0.5         0.7
#> 1439         2  27650.082    final        1200       1, 1   0.3 0.5         0.7
#> 1440         2  27753.615    final        1200       1, 1   0.3 0.5         0.7
#> 1441         2  27968.984    final        1200       1, 1   0.3 0.5         0.7
#> 1442         2  27991.364    final        1200       1, 1   0.3 0.5         0.7
#> 1443         2  28032.754    final        1200       1, 1   0.3 0.5         0.7
#> 1444         2  28118.856    final        1200       1, 1   0.3 0.5         0.7
#> 1445         2  28194.724    final        1200       1, 1   0.3 0.5         0.7
#> 1446         2  28206.383    final        1200       1, 1   0.3 0.5         0.7
#> 1447         2  28234.938    final        1200       1, 1   0.3 0.5         0.7
#> 1448         2  28322.330    final        1200       1, 1   0.3 0.5         0.7
#> 1449         2  28384.297    final        1200       1, 1   0.3 0.5         0.7
#> 1450         2  28398.676    final        1200       1, 1   0.3 0.5         0.7
#> 1451         2  28610.075    final        1200       1, 1   0.3 0.5         0.7
#> 1452         2  28926.387    final        1200       1, 1   0.3 0.5         0.7
#> 1453         2  28929.395    final        1200       1, 1   0.3 0.5         0.7
#> 1454         2  28989.341    final        1200       1, 1   0.3 0.5         0.7
#> 1455         2  29117.003    final        1200       1, 1   0.3 0.5         0.7
#> 1456         2  29370.939    final        1200       1, 1   0.3 0.5         0.7
#> 1457         2  29473.145    final        1200       1, 1   0.3 0.5         0.7
#> 1458         2  29566.822    final        1200       1, 1   0.3 0.5         0.7
#> 1459         2  29653.407    final        1200       1, 1   0.3 0.5         0.7
#> 1460         2  29654.406    final        1200       1, 1   0.3 0.5         0.7
#> 1461         2  29686.899    final        1200       1, 1   0.3 0.5         0.7
#> 1462         2  29720.671    final        1200       1, 1   0.3 0.5         0.7
#> 1463         2  30170.946    final        1200       1, 1   0.3 0.5         0.7
#> 1464         2  30247.724    final        1200       1, 1   0.3 0.5         0.7
#> 1465         2  30445.881    final        1200       1, 1   0.3 0.5         0.7
#> 1466         2  30450.698    final        1200       1, 1   0.3 0.5         0.7
#> 1467         2  30687.506    final        1200       1, 1   0.3 0.5         0.7
#> 1468         2  30875.333    final        1200       1, 1   0.3 0.5         0.7
#> 1469         2  30968.373    final        1200       1, 1   0.3 0.5         0.7
#> 1470         2  31120.471    final        1200       1, 1   0.3 0.5         0.7
#> 1471         2  31204.514    final        1200       1, 1   0.3 0.5         0.7
#> 1472         2  31373.051    final        1200       1, 1   0.3 0.5         0.7
#> 1473         2  31484.054    final        1200       1, 1   0.3 0.5         0.7
#> 1474         2  31613.065    final        1200       1, 1   0.3 0.5         0.7
#> 1475         2  31758.915    final        1200       1, 1   0.3 0.5         0.7
#> 1476         2  31815.925    final        1200       1, 1   0.3 0.5         0.7
#> 1477         2  31974.858    final        1200       1, 1   0.3 0.5         0.7
#> 1478         2  32273.901    final        1200       1, 1   0.3 0.5         0.7
#> 1479         2  32357.797    final        1200       1, 1   0.3 0.5         0.7
#> 1480         2  32503.476    final        1200       1, 1   0.3 0.5         0.7
#> 1481         2  32544.025    final        1200       1, 1   0.3 0.5         0.7
#> 1482         2  32601.269    final        1200       1, 1   0.3 0.5         0.7
#> 1483         2  32712.650    final        1200       1, 1   0.3 0.5         0.7
#> 1484         2  32882.992    final        1200       1, 1   0.3 0.5         0.7
#> 1485         2  32955.969    final        1200       1, 1   0.3 0.5         0.7
#> 1486         2  33180.961    final        1200       1, 1   0.3 0.5         0.7
#> 1487         2  33187.850    final        1200       1, 1   0.3 0.5         0.7
#> 1488         2  33333.955    final        1200       1, 1   0.3 0.5         0.7
#> 1489         2  33450.975    final        1200       1, 1   0.3 0.5         0.7
#> 1490         2  33808.059    final        1200       1, 1   0.3 0.5         0.7
#> 1491         2  33832.004    final        1200       1, 1   0.3 0.5         0.7
#> 1492         2  33888.531    final        1200       1, 1   0.3 0.5         0.7
#> 1493         2  33930.618    final        1200       1, 1   0.3 0.5         0.7
#> 1494         2  33992.853    final        1200       1, 1   0.3 0.5         0.7
#> 1495         2  33998.788    final        1200       1, 1   0.3 0.5         0.7
#> 1496         2  34184.294    final        1200       1, 1   0.3 0.5         0.7
#> 1497         2  34294.762    final        1200       1, 1   0.3 0.5         0.7
#> 1498         2  34426.465    final        1200       1, 1   0.3 0.5         0.7
#> 1499         2  34505.999    final        1200       1, 1   0.3 0.5         0.7
#> 1500         2  34540.039    final        1200       1, 1   0.3 0.5         0.7
#> 1501         2  34629.458    final        1200       1, 1   0.3 0.5         0.7
#> 1502         2  34927.301    final        1200       1, 1   0.3 0.5         0.7
#> 1503         2  35006.444    final        1200       1, 1   0.3 0.5         0.7
#> 1504         2  35041.947    final        1200       1, 1   0.3 0.5         0.7
#> 1505         2  35247.916    final        1200       1, 1   0.3 0.5         0.7
#> 1506         2  35383.835    final        1200       1, 1   0.3 0.5         0.7
#> 1507         2  35505.259    final        1200       1, 1   0.3 0.5         0.7
#> 1508         2  35586.853    final        1200       1, 1   0.3 0.5         0.7
#> 1509         2  35663.735    final        1200       1, 1   0.3 0.5         0.7
#> 1510         2  35834.088    final        1200       1, 1   0.3 0.5         0.7
#> 1511         2  35936.129    final        1200       1, 1   0.3 0.5         0.7
#> 1512         2  36162.898    final        1200       1, 1   0.3 0.5         0.7
#> 1513         2  36248.399    final        1200       1, 1   0.3 0.5         0.7
#> 1514         2  36277.335    final        1200       1, 1   0.3 0.5         0.7
#> 1515         2  36474.620    final        1200       1, 1   0.3 0.5         0.7
#> 1516         2  36749.043    final        1200       1, 1   0.3 0.5         0.7
#> 1517         2  36789.067    final        1200       1, 1   0.3 0.5         0.7
#> 1518         2  36964.109    final        1200       1, 1   0.3 0.5         0.7
#> 1519         2  36980.604    final        1200       1, 1   0.3 0.5         0.7
#> 1520         2  37022.045    final        1200       1, 1   0.3 0.5         0.7
#> 1521         2  37068.447    final        1200       1, 1   0.3 0.5         0.7
#> 1522         2  37235.002    final        1200       1, 1   0.3 0.5         0.7
#> 1523         2  37349.855    final        1200       1, 1   0.3 0.5         0.7
#> 1524         2  37592.579    final        1200       1, 1   0.3 0.5         0.7
#> 1525         2  37598.955    final        1200       1, 1   0.3 0.5         0.7
#> 1526         2  37602.869    final        1200       1, 1   0.3 0.5         0.7
#> 1527         2  37876.974    final        1200       1, 1   0.3 0.5         0.7
#> 1528         2  37895.764    final        1200       1, 1   0.3 0.5         0.7
#> 1529         2  38005.742    final        1200       1, 1   0.3 0.5         0.7
#> 1530         2  38082.118    final        1200       1, 1   0.3 0.5         0.7
#> 1531         2  38119.021    final        1200       1, 1   0.3 0.5         0.7
#> 1532         2  38433.891    final        1200       1, 1   0.3 0.5         0.7
#> 1533         2  38486.073    final        1200       1, 1   0.3 0.5         0.7
#> 1534         2  38553.146    final        1200       1, 1   0.3 0.5         0.7
#> 1535         2  38689.328    final        1200       1, 1   0.3 0.5         0.7
#> 1536         2  38701.607    final        1200       1, 1   0.3 0.5         0.7
#> 1537         2  38718.562    final        1200       1, 1   0.3 0.5         0.7
#> 1538         2  38751.551    final        1200       1, 1   0.3 0.5         0.7
#> 1539         2  38810.343    final        1200       1, 1   0.3 0.5         0.7
#> 1540         2  38928.121    final        1200       1, 1   0.3 0.5         0.7
#> 1541         2  39069.190    final        1200       1, 1   0.3 0.5         0.7
#> 1542         2  39108.759    final        1200       1, 1   0.3 0.5         0.7
#> 1543         2  39108.775    final        1200       1, 1   0.3 0.5         0.7
#> 1544         2  39181.851    final        1200       1, 1   0.3 0.5         0.7
#> 1545         2  39194.156    final        1200       1, 1   0.3 0.5         0.7
#> 1546         2  39444.353    final        1200       1, 1   0.3 0.5         0.7
#> 1547         2  39650.425    final        1200       1, 1   0.3 0.5         0.7
#> 1548         2  39723.905    final        1200       1, 1   0.3 0.5         0.7
#> 1549         2  39754.489    final        1200       1, 1   0.3 0.5         0.7
#> 1550         2  40021.615    final        1200       1, 1   0.3 0.5         0.7
#> 1551         2  40068.554    final        1200       1, 1   0.3 0.5         0.7
#> 1552         2  40330.332    final        1200       1, 1   0.3 0.5         0.7
#> 1553         2  40472.578    final        1200       1, 1   0.3 0.5         0.7
#> 1554         2  40480.569    final        1200       1, 1   0.3 0.5         0.7
#> 1555         2  40550.901    final        1200       1, 1   0.3 0.5         0.7
#> 1556         2  40599.165    final        1200       1, 1   0.3 0.5         0.7
#> 1557         2  40749.568    final        1200       1, 1   0.3 0.5         0.7
#> 1558         2  40785.177    final        1200       1, 1   0.3 0.5         0.7
#> 1559         2  40820.330    final        1200       1, 1   0.3 0.5         0.7
#> 1560         2  40848.741    final        1200       1, 1   0.3 0.5         0.7
#> 1561         2  40862.981    final        1200       1, 1   0.3 0.5         0.7
#> 1562         2  40872.750    final        1200       1, 1   0.3 0.5         0.7
#> 1563         2  40976.612    final        1200       1, 1   0.3 0.5         0.7
#> 1564         2  41010.679    final        1200       1, 1   0.3 0.5         0.7
#> 1565         2  41040.068    final        1200       1, 1   0.3 0.5         0.7
#> 1566         2  41062.004    final        1200       1, 1   0.3 0.5         0.7
#> 1567         2  41378.984    final        1200       1, 1   0.3 0.5         0.7
#> 1568         2  41499.856    final        1200       1, 1   0.3 0.5         0.7
#> 1569         2  41677.653    final        1200       1, 1   0.3 0.5         0.7
#> 1570         2  41893.687    final        1200       1, 1   0.3 0.5         0.7
#> 1571         2  41953.281    final        1200       1, 1   0.3 0.5         0.7
#> 1572         2  41997.313    final        1200       1, 1   0.3 0.5         0.7
#> 1573         2  42090.326    final        1200       1, 1   0.3 0.5         0.7
#> 1574         2  42101.117    final        1200       1, 1   0.3 0.5         0.7
#> 1575         2  42135.651    final        1200       1, 1   0.3 0.5         0.7
#> 1576         2  42273.737    final        1200       1, 1   0.3 0.5         0.7
#> 1577         2  42426.834    final        1200       1, 1   0.3 0.5         0.7
#> 1578         2  42429.052    final        1200       1, 1   0.3 0.5         0.7
#> 1579         2  42552.469    final        1200       1, 1   0.3 0.5         0.7
#> 1580         2  42590.039    final        1200       1, 1   0.3 0.5         0.7
#> 1581         2  42591.902    final        1200       1, 1   0.3 0.5         0.7
#> 1582         2  42683.684    final        1200       1, 1   0.3 0.5         0.7
#> 1583         2  43057.228    final        1200       1, 1   0.3 0.5         0.7
#> 1584         2  43150.653    final        1200       1, 1   0.3 0.5         0.7
#> 1585         2  43277.764    final        1200       1, 1   0.3 0.5         0.7
#> 1586         2  43345.288    final        1200       1, 1   0.3 0.5         0.7
#> 1587         2  43466.221    final        1200       1, 1   0.3 0.5         0.7
#> 1588         2  43497.403    final        1200       1, 1   0.3 0.5         0.7
#> 1589         2  43531.781    final        1200       1, 1   0.3 0.5         0.7
#> 1590         2  43594.639    final        1200       1, 1   0.3 0.5         0.7
#> 1591         2  43684.648    final        1200       1, 1   0.3 0.5         0.7
#> 1592         2  43741.384    final        1200       1, 1   0.3 0.5         0.7
#> 1593         2  43761.595    final        1200       1, 1   0.3 0.5         0.7
#> 1594         2  43770.067    final        1200       1, 1   0.3 0.5         0.7
#> 1595         2  43841.939    final        1200       1, 1   0.3 0.5         0.7
#> 1596         2  44277.149    final        1200       1, 1   0.3 0.5         0.7
#> 1597         2  44349.737    final        1200       1, 1   0.3 0.5         0.7
#> 1598         2  44352.687    final        1200       1, 1   0.3 0.5         0.7
#> 1599         2  44467.441    final        1200       1, 1   0.3 0.5         0.7
#> 1600         2  44475.532    final        1200       1, 1   0.3 0.5         0.7
#> 1601         2  44857.707    final        1200       1, 1   0.3 0.5         0.7
#> 1602         2  44940.282    final        1200       1, 1   0.3 0.5         0.7
#> 1603         2  45159.736    final        1200       1, 1   0.3 0.5         0.7
#> 1604         2  45165.684    final        1200       1, 1   0.3 0.5         0.7
#> 1605         2  45165.840    final        1200       1, 1   0.3 0.5         0.7
#> 1606         2  45186.631    final        1200       1, 1   0.3 0.5         0.7
#> 1607         2  45187.442    final        1200       1, 1   0.3 0.5         0.7
#> 1608         2  45202.210    final        1200       1, 1   0.3 0.5         0.7
#> 1609         2  45268.321    final        1200       1, 1   0.3 0.5         0.7
#> 1610         2  45333.784    final        1200       1, 1   0.3 0.5         0.7
#> 1611         2  45342.159    final        1200       1, 1   0.3 0.5         0.7
#> 1612         2  45374.552    final        1200       1, 1   0.3 0.5         0.7
#> 1613         2  45488.506    final        1200       1, 1   0.3 0.5         0.7
#> 1614         2  45661.652    final        1200       1, 1   0.3 0.5         0.7
#> 1615         2  45771.381    final        1200       1, 1   0.3 0.5         0.7
#> 1616         2  45887.398    final        1200       1, 1   0.3 0.5         0.7
#> 1617         2  46044.776    final        1200       1, 1   0.3 0.5         0.7
#> 1618         2  46134.210    final        1200       1, 1   0.3 0.5         0.7
#> 1619         2  46491.123    final        1200       1, 1   0.3 0.5         0.7
#> 1620         2  46556.968    final        1200       1, 1   0.3 0.5         0.7
#> 1621         2  46739.226    final        1200       1, 1   0.3 0.5         0.7
#> 1622         2  46933.620    final        1200       1, 1   0.3 0.5         0.7
#> 1623         2  47057.592    final        1200       1, 1   0.3 0.5         0.7
#> 1624         2  47187.215    final        1200       1, 1   0.3 0.5         0.7
#> 1625         2  47188.282    final        1200       1, 1   0.3 0.5         0.7
#> 1626         2  47325.386    final        1200       1, 1   0.3 0.5         0.7
#> 1627         2  47485.047    final        1200       1, 1   0.3 0.5         0.7
#> 1628         2  47605.219    final        1200       1, 1   0.3 0.5         0.7
#> 1629         2  47607.161    final        1200       1, 1   0.3 0.5         0.7
#> 1630         2  47626.684    final        1200       1, 1   0.3 0.5         0.7
#> 1631         2  47632.444    final        1200       1, 1   0.3 0.5         0.7
#> 1632         2  47734.347    final        1200       1, 1   0.3 0.5         0.7
#> 1633         2  47774.091    final        1200       1, 1   0.3 0.5         0.7
#> 1634         2  48052.588    final        1200       1, 1   0.3 0.5         0.7
#> 1635         2  48065.644    final        1200       1, 1   0.3 0.5         0.7
#> 1636         2  48101.166    final        1200       1, 1   0.3 0.5         0.7
#> 1637         2  48134.626    final        1200       1, 1   0.3 0.5         0.7
#> 1638         2  48483.683    final        1200       1, 1   0.3 0.5         0.7
#> 1639         2  48650.199    final        1200       1, 1   0.3 0.5         0.7
#> 1640         2  48783.351    final        1200       1, 1   0.3 0.5         0.7
#> 1641         2  48793.580    final        1200       1, 1   0.3 0.5         0.7
#> 1642         2  48816.196    final        1200       1, 1   0.3 0.5         0.7
#> 1643         2  48821.666    final        1200       1, 1   0.3 0.5         0.7
#> 1644         2  48852.930    final        1200       1, 1   0.3 0.5         0.7
#> 1645         2  48948.638    final        1200       1, 1   0.3 0.5         0.7
#> 1646         2  48967.915    final        1200       1, 1   0.3 0.5         0.7
#> 1647         2  48998.343    final        1200       1, 1   0.3 0.5         0.7
#> 1648         2  49392.186    final        1200       1, 1   0.3 0.5         0.7
#> 1649         2  49491.594    final        1200       1, 1   0.3 0.5         0.7
#> 1650         2  49524.550    final        1200       1, 1   0.3 0.5         0.7
#> 1651         2  49652.800    final        1200       1, 1   0.3 0.5         0.7
#> 1652         2  49761.210    final        1200       1, 1   0.3 0.5         0.7
#> 1653         2  49805.779    final        1200       1, 1   0.3 0.5         0.7
#> 1654         2  49827.496    final        1200       1, 1   0.3 0.5         0.7
#> 1655         2  49853.082    final        1200       1, 1   0.3 0.5         0.7
#> 1656         2  49904.092    final        1200       1, 1   0.3 0.5         0.7
#> 1657         2  49931.572    final        1200       1, 1   0.3 0.5         0.7
#> 1658         2  50271.300    final        1200       1, 1   0.3 0.5         0.7
#> 1659         2  50892.064    final        1200       1, 1   0.3 0.5         0.7
#> 1660         2  51065.711    final        1200       1, 1   0.3 0.5         0.7
#> 1661         2  51277.931    final        1200       1, 1   0.3 0.5         0.7
#> 1662         2  51287.657    final        1200       1, 1   0.3 0.5         0.7
#> 1663         2  51602.235    final        1200       1, 1   0.3 0.5         0.7
#> 1664         2  51744.279    final        1200       1, 1   0.3 0.5         0.7
#> 1665         2  51879.589    final        1200       1, 1   0.3 0.5         0.7
#> 1666         2  51956.794    final        1200       1, 1   0.3 0.5         0.7
#> 1667         2  52095.253    final        1200       1, 1   0.3 0.5         0.7
#> 1668         2  52240.440    final        1200       1, 1   0.3 0.5         0.7
#> 1669         2  52270.225    final        1200       1, 1   0.3 0.5         0.7
#> 1670         2  52304.480    final        1200       1, 1   0.3 0.5         0.7
#> 1671         2  52356.607    final        1200       1, 1   0.3 0.5         0.7
#> 1672         2  52367.433    final        1200       1, 1   0.3 0.5         0.7
#> 1673         2  52368.079    final        1200       1, 1   0.3 0.5         0.7
#> 1674         2  52473.999    final        1200       1, 1   0.3 0.5         0.7
#> 1675         2  52521.809    final        1200       1, 1   0.3 0.5         0.7
#> 1676         2  52702.148    final        1200       1, 1   0.3 0.5         0.7
#> 1677         2  52840.202    final        1200       1, 1   0.3 0.5         0.7
#> 1678         2  52931.382    final        1200       1, 1   0.3 0.5         0.7
#> 1679         2  53189.474    final        1200       1, 1   0.3 0.5         0.7
#> 1680         2  53316.246    final        1200       1, 1   0.3 0.5         0.7
#> 1681         2  53380.122    final        1200       1, 1   0.3 0.5         0.7
#> 1682         2  53393.330    final        1200       1, 1   0.3 0.5         0.7
#> 1683         2  53427.907    final        1200       1, 1   0.3 0.5         0.7
#> 1684         2  53439.029    final        1200       1, 1   0.3 0.5         0.7
#> 1685         2  53445.772    final        1200       1, 1   0.3 0.5         0.7
#> 1686         2  53572.232    final        1200       1, 1   0.3 0.5         0.7
#> 1687         2  53597.539    final        1200       1, 1   0.3 0.5         0.7
#> 1688         2  53979.541    final        1200       1, 1   0.3 0.5         0.7
#> 1689         2  54004.116    final        1200       1, 1   0.3 0.5         0.7
#> 1690         2  54101.110    final        1200       1, 1   0.3 0.5         0.7
#> 1691         2  54242.781    final        1200       1, 1   0.3 0.5         0.7
#> 1692         2  54256.882    final        1200       1, 1   0.3 0.5         0.7
#> 1693         2  54416.639    final        1200       1, 1   0.3 0.5         0.7
#> 1694         2  54461.352    final        1200       1, 1   0.3 0.5         0.7
#> 1695         2  54780.818    final        1200       1, 1   0.3 0.5         0.7
#> 1696         2  54827.456    final        1200       1, 1   0.3 0.5         0.7
#> 1697         2  54955.169    final        1200       1, 1   0.3 0.5         0.7
#> 1698         2  54958.806    final        1200       1, 1   0.3 0.5         0.7
#> 1699         2  55133.687    final        1200       1, 1   0.3 0.5         0.7
#> 1700         2  55354.385    final        1200       1, 1   0.3 0.5         0.7
#> 1701         2  55676.742    final        1200       1, 1   0.3 0.5         0.7
#> 1702         2  55765.953    final        1200       1, 1   0.3 0.5         0.7
#> 1703         2  55926.190    final        1200       1, 1   0.3 0.5         0.7
#> 1704         2  55956.322    final        1200       1, 1   0.3 0.5         0.7
#> 1705         2  56091.607    final        1200       1, 1   0.3 0.5         0.7
#> 1706         2  56474.216    final        1200       1, 1   0.3 0.5         0.7
#> 1707         2  56769.945    final        1200       1, 1   0.3 0.5         0.7
#> 1708         2  56821.055    final        1200       1, 1   0.3 0.5         0.7
#> 1709         2  56859.213    final        1200       1, 1   0.3 0.5         0.7
#> 1710         2  56937.640    final        1200       1, 1   0.3 0.5         0.7
#> 1711         2  56938.542    final        1200       1, 1   0.3 0.5         0.7
#> 1712         2  56988.125    final        1200       1, 1   0.3 0.5         0.7
#> 1713         2  57002.889    final        1200       1, 1   0.3 0.5         0.7
#> 1714         2  57015.634    final        1200       1, 1   0.3 0.5         0.7
#> 1715         2  57056.380    final        1200       1, 1   0.3 0.5         0.7
#> 1716         2  57080.470    final        1200       1, 1   0.3 0.5         0.7
#> 1717         2  57147.254    final        1200       1, 1   0.3 0.5         0.7
#> 1718         2  57385.264    final        1200       1, 1   0.3 0.5         0.7
#> 1719         2  57656.537    final        1200       1, 1   0.3 0.5         0.7
#> 1720         2  57670.917    final        1200       1, 1   0.3 0.5         0.7
#> 1721         2  57763.529    final        1200       1, 1   0.3 0.5         0.7
#> 1722         2  57763.727    final        1200       1, 1   0.3 0.5         0.7
#> 1723         2  57843.376    final        1200       1, 1   0.3 0.5         0.7
#> 1724         2  57854.728    final        1200       1, 1   0.3 0.5         0.7
#> 1725         2  57870.672    final        1200       1, 1   0.3 0.5         0.7
#> 1726         2  58016.163    final        1200       1, 1   0.3 0.5         0.7
#> 1727         2  58136.706    final        1200       1, 1   0.3 0.5         0.7
#> 1728         2  58380.743    final        1200       1, 1   0.3 0.5         0.7
#> 1729         2  58442.074    final        1200       1, 1   0.3 0.5         0.7
#> 1730         2  58589.012    final        1200       1, 1   0.3 0.5         0.7
#> 1731         2  58676.453    final        1200       1, 1   0.3 0.5         0.7
#> 1732         2  58696.787    final        1200       1, 1   0.3 0.5         0.7
#> 1733         2  58740.270    final        1200       1, 1   0.3 0.5         0.7
#> 1734         2  58747.556    final        1200       1, 1   0.3 0.5         0.7
#> 1735         2  58815.365    final        1200       1, 1   0.3 0.5         0.7
#> 1736         2  58846.126    final        1200       1, 1   0.3 0.5         0.7
#> 1737         2  58904.836    final        1200       1, 1   0.3 0.5         0.7
#> 1738         2  59100.021    final        1200       1, 1   0.3 0.5         0.7
#> 1739         2  59129.057    final        1200       1, 1   0.3 0.5         0.7
#> 1740         2  59271.156    final        1200       1, 1   0.3 0.5         0.7
#> 1741         2  59296.972    final        1200       1, 1   0.3 0.5         0.7
#> 1742         2  59298.834    final        1200       1, 1   0.3 0.5         0.7
#> 1743         2  59399.730    final        1200       1, 1   0.3 0.5         0.7
#> 1744         2  59448.364    final        1200       1, 1   0.3 0.5         0.7
#> 1745         2  59474.695    final        1200       1, 1   0.3 0.5         0.7
#> 1746         2  59544.270    final        1200       1, 1   0.3 0.5         0.7
#> 1747         2  59570.994    final        1200       1, 1   0.3 0.5         0.7
#> 1748         2  59637.526    final        1200       1, 1   0.3 0.5         0.7
#> 1749         2  59844.922    final        1200       1, 1   0.3 0.5         0.7
#> 1750         2  59848.470    final        1200       1, 1   0.3 0.5         0.7
#> 1751         2  59926.512    final        1200       1, 1   0.3 0.5         0.7
#> 1752         2  60013.846    final        1200       1, 1   0.3 0.5         0.7
#> 1753         2  60151.370    final        1200       1, 1   0.3 0.5         0.7
#> 1754         2  60206.672    final        1200       1, 1   0.3 0.5         0.7
#> 1755         2  60215.999    final        1200       1, 1   0.3 0.5         0.7
#> 1756         2  60356.780    final        1200       1, 1   0.3 0.5         0.7
#> 1757         2  60357.989    final        1200       1, 1   0.3 0.5         0.7
#> 1758         2  60402.436    final        1200       1, 1   0.3 0.5         0.7
#> 1759         2  60428.563    final        1200       1, 1   0.3 0.5         0.7
#> 1760         2  60695.566    final        1200       1, 1   0.3 0.5         0.7
#> 1761         2  60828.300    final        1200       1, 1   0.3 0.5         0.7
#> 1762         2  61000.256    final        1200       1, 1   0.3 0.5         0.7
#> 1763         2  61126.402    final        1200       1, 1   0.3 0.5         0.7
#> 1764         2  61140.647    final        1200       1, 1   0.3 0.5         0.7
#> 1765         2  61154.122    final        1200       1, 1   0.3 0.5         0.7
#> 1766         2  61158.846    final        1200       1, 1   0.3 0.5         0.7
#> 1767         2  61296.816    final        1200       1, 1   0.3 0.5         0.7
#> 1768         2  61377.491    final        1200       1, 1   0.3 0.5         0.7
#> 1769         2  61433.906    final        1200       1, 1   0.3 0.5         0.7
#> 1770         2  61490.008    final        1200       1, 1   0.3 0.5         0.7
#> 1771         2  61558.541    final        1200       1, 1   0.3 0.5         0.7
#> 1772         2  61667.187    final        1200       1, 1   0.3 0.5         0.7
#> 1773         2  61816.890    final        1200       1, 1   0.3 0.5         0.7
#> 1774         2  61986.094    final        1200       1, 1   0.3 0.5         0.7
#> 1775         2  61990.462    final        1200       1, 1   0.3 0.5         0.7
#> 1776         2  62371.884    final        1200       1, 1   0.3 0.5         0.7
#> 1777         2  62476.283    final        1200       1, 1   0.3 0.5         0.7
#> 1778         2  62486.900    final        1200       1, 1   0.3 0.5         0.7
#> 1779         2  62521.883    final        1200       1, 1   0.3 0.5         0.7
#> 1780         2  62633.263    final        1200       1, 1   0.3 0.5         0.7
#> 1781         2  62756.147    final        1200       1, 1   0.3 0.5         0.7
#> 1782         2  62868.338    final        1200       1, 1   0.3 0.5         0.7
#> 1783         2  62934.869    final        1200       1, 1   0.3 0.5         0.7
#> 1784         2  63029.093    final        1200       1, 1   0.3 0.5         0.7
#> 1785         2  63073.835    final        1200       1, 1   0.3 0.5         0.7
#> 1786         2  63095.034    final        1200       1, 1   0.3 0.5         0.7
#> 1787         2  63113.769    final        1200       1, 1   0.3 0.5         0.7
#> 1788         2  63284.851    final        1200       1, 1   0.3 0.5         0.7
#> 1789         2  63291.033    final        1200       1, 1   0.3 0.5         0.7
#> 1790         2  63322.984    final        1200       1, 1   0.3 0.5         0.7
#> 1791         2  63527.998    final        1200       1, 1   0.3 0.5         0.7
#> 1792         2  63611.198    final        1200       1, 1   0.3 0.5         0.7
#> 1793         2  63734.384    final        1200       1, 1   0.3 0.5         0.7
#> 1794         2  63821.816    final        1200       1, 1   0.3 0.5         0.7
#> 1795         2  63837.668    final        1200       1, 1   0.3 0.5         0.7
#> 1796         2  63926.975    final        1200       1, 1   0.3 0.5         0.7
#> 1797         2  64057.518    final        1200       1, 1   0.3 0.5         0.7
#> 1798         2  64077.282    final        1200       1, 1   0.3 0.5         0.7
#> 1799         2  64084.603    final        1200       1, 1   0.3 0.5         0.7
#> 1800         2  64221.582    final        1200       1, 1   0.3 0.5         0.7
#> 1801         2  64240.412    final        1200       1, 1   0.3 0.5         0.7
#> 1802         2  64364.869    final        1200       1, 1   0.3 0.5         0.7
#> 1803         2  64441.476    final        1200       1, 1   0.3 0.5         0.7
#> 1804         2  64507.911    final        1200       1, 1   0.3 0.5         0.7
#> 1805         2  64520.372    final        1200       1, 1   0.3 0.5         0.7
#> 1806         2  64540.801    final        1200       1, 1   0.3 0.5         0.7
#> 1807         2  64802.128    final        1200       1, 1   0.3 0.5         0.7
#> 1808         2  64919.087    final        1200       1, 1   0.3 0.5         0.7
#> 1809         2  64996.037    final        1200       1, 1   0.3 0.5         0.7
#> 1810         2  65001.878    final        1200       1, 1   0.3 0.5         0.7
#> 1811         2  65066.813    final        1200       1, 1   0.3 0.5         0.7
#> 1812         2  65127.979    final        1200       1, 1   0.3 0.5         0.7
#> 1813         2  65191.262    final        1200       1, 1   0.3 0.5         0.7
#> 1814         2  65280.965    final        1200       1, 1   0.3 0.5         0.7
#> 1815         2  65295.307    final        1200       1, 1   0.3 0.5         0.7
#> 1816         2  65356.650    final        1200       1, 1   0.3 0.5         0.7
#> 1817         2  65438.353    final        1200       1, 1   0.3 0.5         0.7
#> 1818         2  65529.601    final        1200       1, 1   0.3 0.5         0.7
#> 1819         2  65752.895    final        1200       1, 1   0.3 0.5         0.7
#> 1820         2  65966.974    final        1200       1, 1   0.3 0.5         0.7
#> 1821         2  65990.705    final        1200       1, 1   0.3 0.5         0.7
#> 1822         2  66103.136    final        1200       1, 1   0.3 0.5         0.7
#> 1823         2  66129.012    final        1200       1, 1   0.3 0.5         0.7
#> 1824         2  66315.461    final        1200       1, 1   0.3 0.5         0.7
#> 1825         2  66329.342    final        1200       1, 1   0.3 0.5         0.7
#> 1826         2  66357.183    final        1200       1, 1   0.3 0.5         0.7
#> 1827         2  66376.746    final        1200       1, 1   0.3 0.5         0.7
#> 1828         2  66483.913    final        1200       1, 1   0.3 0.5         0.7
#> 1829         2  66508.889    final        1200       1, 1   0.3 0.5         0.7
#> 1830         2  66531.904    final        1200       1, 1   0.3 0.5         0.7
#> 1831         2  66564.910    final        1200       1, 1   0.3 0.5         0.7
#> 1832         2  66636.757    final        1200       1, 1   0.3 0.5         0.7
#> 1833         2  66659.283    final        1200       1, 1   0.3 0.5         0.7
#> 1834         2  66921.757    final        1200       1, 1   0.3 0.5         0.7
#> 1835         2  66925.290    final        1200       1, 1   0.3 0.5         0.7
#> 1836         2  67082.367    final        1200       1, 1   0.3 0.5         0.7
#> 1837         2  67135.977    final        1200       1, 1   0.3 0.5         0.7
#> 1838         2  67222.159    final        1200       1, 1   0.3 0.5         0.7
#> 1839         2  67325.922    final        1200       1, 1   0.3 0.5         0.7
#> 1840         2  67344.063    final        1200       1, 1   0.3 0.5         0.7
#> 1841         2  67396.029    final        1200       1, 1   0.3 0.5         0.7
#> 1842         2  67549.498    final        1200       1, 1   0.3 0.5         0.7
#> 1843         2  67569.663    final        1200       1, 1   0.3 0.5         0.7
#> 1844         2  67666.246    final        1200       1, 1   0.3 0.5         0.7
#> 1845         2  67695.859    final        1200       1, 1   0.3 0.5         0.7
#> 1846         2  67781.945    final        1200       1, 1   0.3 0.5         0.7
#> 1847         2  67823.632    final        1200       1, 1   0.3 0.5         0.7
#> 1848         2  67871.214    final        1200       1, 1   0.3 0.5         0.7
#> 1849         2  68010.915    final        1200       1, 1   0.3 0.5         0.7
#> 1850         2  68045.761    final        1200       1, 1   0.3 0.5         0.7
#> 1851         2  68182.533    final        1200       1, 1   0.3 0.5         0.7
#> 1852         2  68321.415    final        1200       1, 1   0.3 0.5         0.7
#> 1853         2  68344.474    final        1200       1, 1   0.3 0.5         0.7
#> 1854         2  68392.747    final        1200       1, 1   0.3 0.5         0.7
#> 1855         2  68424.341    final        1200       1, 1   0.3 0.5         0.7
#> 1856         2  68520.763    final        1200       1, 1   0.3 0.5         0.7
#> 1857         2  68561.299    final        1200       1, 1   0.3 0.5         0.7
#> 1858         2  68621.644    final        1200       1, 1   0.3 0.5         0.7
#> 1859         2  68774.240    final        1200       1, 1   0.3 0.5         0.7
#> 1860         2  68905.015    final        1200       1, 1   0.3 0.5         0.7
#> 1861         2  68951.375    final        1200       1, 1   0.3 0.5         0.7
#> 1862         2  69011.917    final        1200       1, 1   0.3 0.5         0.7
#> 1863         2  69063.912    final        1200       1, 1   0.3 0.5         0.7
#> 1864         2  69089.020    final        1200       1, 1   0.3 0.5         0.7
#> 1865         2  69185.886    final        1200       1, 1   0.3 0.5         0.7
#> 1866         2  69274.266    final        1200       1, 1   0.3 0.5         0.7
#> 1867         2  69360.586    final        1200       1, 1   0.3 0.5         0.7
#> 1868         2  69549.081    final        1200       1, 1   0.3 0.5         0.7
#> 1869         2  69557.082    final        1200       1, 1   0.3 0.5         0.7
#> 1870         2  69639.648    final        1200       1, 1   0.3 0.5         0.7
#> 1871         2  69682.004    final        1200       1, 1   0.3 0.5         0.7
#> 1872         2  69703.612    final        1200       1, 1   0.3 0.5         0.7
#> 1873         2  69762.028    final        1200       1, 1   0.3 0.5         0.7
#> 1874         2  69944.505    final        1200       1, 1   0.3 0.5         0.7
#> 1875         2  70078.666    final        1200       1, 1   0.3 0.5         0.7
#> 1876         2  70141.291    final        1200       1, 1   0.3 0.5         0.7
#> 1877         2  70190.588    final        1200       1, 1   0.3 0.5         0.7
#> 1878         2  70338.512    final        1200       1, 1   0.3 0.5         0.7
#> 1879         2  70540.399    final        1200       1, 1   0.3 0.5         0.7
#> 1880         2  70695.552    final        1200       1, 1   0.3 0.5         0.7
#> 1881         2  70830.056    final        1200       1, 1   0.3 0.5         0.7
#> 1882         2  70881.880    final        1200       1, 1   0.3 0.5         0.7
#> 1883         2  70959.109    final        1200       1, 1   0.3 0.5         0.7
#> 1884         2  71150.323    final        1200       1, 1   0.3 0.5         0.7
#> 1885         2  71202.234    final        1200       1, 1   0.3 0.5         0.7
#> 1886         2  71342.930    final        1200       1, 1   0.3 0.5         0.7
#> 1887         2  71412.052    final        1200       1, 1   0.3 0.5         0.7
#> 1888         2  71602.645    final        1200       1, 1   0.3 0.5         0.7
#> 1889         2  71610.200    final        1200       1, 1   0.3 0.5         0.7
#> 1890         2  71725.214    final        1200       1, 1   0.3 0.5         0.7
#> 1891         2  71855.053    final        1200       1, 1   0.3 0.5         0.7
#> 1892         2  71916.384    final        1200       1, 1   0.3 0.5         0.7
#> 1893         2  71988.074    final        1200       1, 1   0.3 0.5         0.7
#> 1894         2  72162.745    final        1200       1, 1   0.3 0.5         0.7
#> 1895         2  72194.739    final        1200       1, 1   0.3 0.5         0.7
#> 1896         2  72201.462    final        1200       1, 1   0.3 0.5         0.7
#> 1897         2  72251.660    final        1200       1, 1   0.3 0.5         0.7
#> 1898         2  72303.420    final        1200       1, 1   0.3 0.5         0.7
#> 1899         2  72371.518    final        1200       1, 1   0.3 0.5         0.7
#> 1900         2  72411.089    final        1200       1, 1   0.3 0.5         0.7
#> 1901         2  72474.622    final        1200       1, 1   0.3 0.5         0.7
#> 1902         2  72959.070    final        1200       1, 1   0.3 0.5         0.7
#> 1903         2  73035.266    final        1200       1, 1   0.3 0.5         0.7
#> 1904         2  73154.963    final        1200       1, 1   0.3 0.5         0.7
#> 1905         2  73184.961    final        1200       1, 1   0.3 0.5         0.7
#> 1906         2  73202.539    final        1200       1, 1   0.3 0.5         0.7
#> 1907         2  73233.495    final        1200       1, 1   0.3 0.5         0.7
#> 1908         2  73233.781    final        1200       1, 1   0.3 0.5         0.7
#> 1909         2  73352.910    final        1200       1, 1   0.3 0.5         0.7
#> 1910         2  73490.140    final        1200       1, 1   0.3 0.5         0.7
#> 1911         2  73530.540    final        1200       1, 1   0.3 0.5         0.7
#> 1912         2  73545.941    final        1200       1, 1   0.3 0.5         0.7
#> 1913         2  73726.318    final        1200       1, 1   0.3 0.5         0.7
#> 1914         2  73905.790    final        1200       1, 1   0.3 0.5         0.7
#> 1915         2  74104.408    final        1200       1, 1   0.3 0.5         0.7
#> 1916         2  74227.339    final        1200       1, 1   0.3 0.5         0.7
#> 1917         2  74254.806    final        1200       1, 1   0.3 0.5         0.7
#> 1918         2  74461.765    final        1200       1, 1   0.3 0.5         0.7
#> 1919         2  74498.934    final        1200       1, 1   0.3 0.5         0.7
#> 1920         2  74521.411    final        1200       1, 1   0.3 0.5         0.7
#> 1921         2  74604.092    final        1200       1, 1   0.3 0.5         0.7
#> 1922         2  74611.331    final        1200       1, 1   0.3 0.5         0.7
#> 1923         2  74862.941    final        1200       1, 1   0.3 0.5         0.7
#> 1924         2  75160.758    final        1200       1, 1   0.3 0.5         0.7
#> 1925         2  75277.243    final        1200       1, 1   0.3 0.5         0.7
#> 1926         2  75348.750    final        1200       1, 1   0.3 0.5         0.7
#> 1927         2  75461.812    final        1200       1, 1   0.3 0.5         0.7
#> 1928         2  75656.882    final        1200       1, 1   0.3 0.5         0.7
#> 1929         2  75844.518    final        1200       1, 1   0.3 0.5         0.7
#> 1930         2  75895.842    final        1200       1, 1   0.3 0.5         0.7
#> 1931         2  75993.987    final        1200       1, 1   0.3 0.5         0.7
#> 1932         2  75999.856    final        1200       1, 1   0.3 0.5         0.7
#> 1933         2  76034.884    final        1200       1, 1   0.3 0.5         0.7
#> 1934         2  76206.184    final        1200       1, 1   0.3 0.5         0.7
#> 1935         2  76244.987    final        1200       1, 1   0.3 0.5         0.7
#> 1936         2  76370.694    final        1200       1, 1   0.3 0.5         0.7
#> 1937         2  76402.586    final        1200       1, 1   0.3 0.5         0.7
#> 1938         2  76483.419    final        1200       1, 1   0.3 0.5         0.7
#> 1939         2  76541.992    final        1200       1, 1   0.3 0.5         0.7
#> 1940         2  76603.643    final        1200       1, 1   0.3 0.5         0.7
#> 1941         2  76618.791    final        1200       1, 1   0.3 0.5         0.7
#> 1942         2  76725.853    final        1200       1, 1   0.3 0.5         0.7
#> 1943         2  76889.531    final        1200       1, 1   0.3 0.5         0.7
#> 1944         2  76898.500    final        1200       1, 1   0.3 0.5         0.7
#> 1945         2  77020.516    final        1200       1, 1   0.3 0.5         0.7
#> 1946         2  77191.666    final        1200       1, 1   0.3 0.5         0.7
#> 1947         2  77376.432    final        1200       1, 1   0.3 0.5         0.7
#> 1948         2  77435.763    final        1200       1, 1   0.3 0.5         0.7
#> 1949         2  77541.515    final        1200       1, 1   0.3 0.5         0.7
#> 1950         2  77551.241    final        1200       1, 1   0.3 0.5         0.7
#> 1951         2  77604.081    final        1200       1, 1   0.3 0.5         0.7
#> 1952         2  77695.702    final        1200       1, 1   0.3 0.5         0.7
#> 1953         2  77815.881    final        1200       1, 1   0.3 0.5         0.7
#> 1954         2  77872.676    final        1200       1, 1   0.3 0.5         0.7
#> 1955         2  78066.351    final        1200       1, 1   0.3 0.5         0.7
#> 1956         2  78135.841    final        1200       1, 1   0.3 0.5         0.7
#> 1957         2  78142.733    final        1200       1, 1   0.3 0.5         0.7
#> 1958         2  78166.443    final        1200       1, 1   0.3 0.5         0.7
#> 1959         2  78308.564    final        1200       1, 1   0.3 0.5         0.7
#> 1960         2  78432.313    final        1200       1, 1   0.3 0.5         0.7
#> 1961         2  78483.849    final        1200       1, 1   0.3 0.5         0.7
#> 1962         2  78516.729    final        1200       1, 1   0.3 0.5         0.7
#> 1963         2  78748.838    final        1200       1, 1   0.3 0.5         0.7
#> 1964         2  78810.756    final        1200       1, 1   0.3 0.5         0.7
#> 1965         2  78839.160    final        1200       1, 1   0.3 0.5         0.7
#> 1966         2  78899.744    final        1200       1, 1   0.3 0.5         0.7
#> 1967         2  79037.183    final        1200       1, 1   0.3 0.5         0.7
#> 1968         2  79126.112    final        1200       1, 1   0.3 0.5         0.7
#> 1969         2  79167.152    final        1200       1, 1   0.3 0.5         0.7
#> 1970         2  79202.058    final        1200       1, 1   0.3 0.5         0.7
#> 1971         2  79273.998    final        1200       1, 1   0.3 0.5         0.7
#> 1972         2  79360.015    final        1200       1, 1   0.3 0.5         0.7
#> 1973         2  79426.894    final        1200       1, 1   0.3 0.5         0.7
#> 1974         2  79697.034    final        1200       1, 1   0.3 0.5         0.7
#> 1975         2  79725.296    final        1200       1, 1   0.3 0.5         0.7
#> 1976         2  79930.218    final        1200       1, 1   0.3 0.5         0.7
#> 1977         2  79959.282    final        1200       1, 1   0.3 0.5         0.7
#> 1978         2  80008.413    final        1200       1, 1   0.3 0.5         0.7
#> 1979         2  80010.451    final        1200       1, 1   0.3 0.5         0.7
#> 1980         2  80185.816    final        1200       1, 1   0.3 0.5         0.7
#> 1981         2  80404.979    final        1200       1, 1   0.3 0.5         0.7
#> 1982         2  80412.444    final        1200       1, 1   0.3 0.5         0.7
#> 1983         2  80504.609    final        1200       1, 1   0.3 0.5         0.7
#> 1984         2  80560.524    final        1200       1, 1   0.3 0.5         0.7
#> 1985         2  80654.516    final        1200       1, 1   0.3 0.5         0.7
#> 1986         2  80723.457    final        1200       1, 1   0.3 0.5         0.7
#> 1987         2  80828.558    final        1200       1, 1   0.3 0.5         0.7
#> 1988         2  80980.344    final        1200       1, 1   0.3 0.5         0.7
#> 1989         2  81100.837    final        1200       1, 1   0.3 0.5         0.7
#> 1990         2  81165.701    final        1200       1, 1   0.3 0.5         0.7
#> 1991         2  81425.604    final        1200       1, 1   0.3 0.5         0.7
#> 1992         2  81485.697    final        1200       1, 1   0.3 0.5         0.7
#> 1993         2  81632.352    final        1200       1, 1   0.3 0.5         0.7
#> 1994         2  81879.003    final        1200       1, 1   0.3 0.5         0.7
#> 1995         2  81897.285    final        1200       1, 1   0.3 0.5         0.7
#> 1996         2  81949.052    final        1200       1, 1   0.3 0.5         0.7
#> 1997         2  82024.180    final        1200       1, 1   0.3 0.5         0.7
#> 1998         2  82135.843    final        1200       1, 1   0.3 0.5         0.7
#> 1999         2  82153.478    final        1200       1, 1   0.3 0.5         0.7
#> 2000         2  82236.721    final        1200       1, 1   0.3 0.5         0.7
#> 2001         2  82279.548    final        1200       1, 1   0.3 0.5         0.7
#> 2002         2  82311.828    final        1200       1, 1   0.3 0.5         0.7
#> 2003         2  82346.218    final        1200       1, 1   0.3 0.5         0.7
#> 2004         2  82412.871    final        1200       1, 1   0.3 0.5         0.7
#> 2005         2  82458.862    final        1200       1, 1   0.3 0.5         0.7
#> 2006         2  82499.516    final        1200       1, 1   0.3 0.5         0.7
#> 2007         2  82603.163    final        1200       1, 1   0.3 0.5         0.7
#> 2008         2  82623.543    final        1200       1, 1   0.3 0.5         0.7
#> 2009         2  82639.075    final        1200       1, 1   0.3 0.5         0.7
#> 2010         2  82687.933    final        1200       1, 1   0.3 0.5         0.7
#> 2011         2  82839.339    final        1200       1, 1   0.3 0.5         0.7
#> 2012         2  83048.446    final        1200       1, 1   0.3 0.5         0.7
#> 2013         2  83117.937    final        1200       1, 1   0.3 0.5         0.7
#> 2014         2  83143.719    final        1200       1, 1   0.3 0.5         0.7
#> 2015         2  83287.315    final        1200       1, 1   0.3 0.5         0.7
#> 2016         2  83466.581    final        1200       1, 1   0.3 0.5         0.7
#> 2017         2  83485.468    final        1200       1, 1   0.3 0.5         0.7
#> 2018         2  83585.284    final        1200       1, 1   0.3 0.5         0.7
#> 2019         2  83612.325    final        1200       1, 1   0.3 0.5         0.7
#> 2020         2  83890.172    final        1200       1, 1   0.3 0.5         0.7
#> 2021         2  83906.994    final        1200       1, 1   0.3 0.5         0.7
#> 2022         2  83986.676    final        1200       1, 1   0.3 0.5         0.7
#> 2023         2  84014.222    final        1200       1, 1   0.3 0.5         0.7
#> 2024         2  84288.122    final        1200       1, 1   0.3 0.5         0.7
#> 2025         2  84547.838    final        1200       1, 1   0.3 0.5         0.7
#> 2026         2  84582.849    final        1200       1, 1   0.3 0.5         0.7
#> 2027         2  84584.674    final        1200       1, 1   0.3 0.5         0.7
#> 2028         2  84670.212    final        1200       1, 1   0.3 0.5         0.7
#> 2029         2  84712.269    final        1200       1, 1   0.3 0.5         0.7
#> 2030         2  84716.808    final        1200       1, 1   0.3 0.5         0.7
#> 2031         2  84752.864    final        1200       1, 1   0.3 0.5         0.7
#> 2032         2  84786.430    final        1200       1, 1   0.3 0.5         0.7
#> 2033         2  84801.116    final        1200       1, 1   0.3 0.5         0.7
#> 2034         2  84817.380    final        1200       1, 1   0.3 0.5         0.7
#> 2035         2  85017.610    final        1200       1, 1   0.3 0.5         0.7
#> 2036         2  85212.343    final        1200       1, 1   0.3 0.5         0.7
#> 2037         2  85295.156    final        1200       1, 1   0.3 0.5         0.7
#> 2038         2  85312.658    final        1200       1, 1   0.3 0.5         0.7
#> 2039         2  85332.278    final        1200       1, 1   0.3 0.5         0.7
#> 2040         2  85465.634    final        1200       1, 1   0.3 0.5         0.7
#> 2041         2  85521.123    final        1200       1, 1   0.3 0.5         0.7
#> 2042         2  85578.044    final        1200       1, 1   0.3 0.5         0.7
#> 2043         2  85615.532    final        1200       1, 1   0.3 0.5         0.7
#> 2044         2  85686.558    final        1200       1, 1   0.3 0.5         0.7
#> 2045         2  85813.972    final        1200       1, 1   0.3 0.5         0.7
#> 2046         2  85844.340    final        1200       1, 1   0.3 0.5         0.7
#> 2047         2  85865.236    final        1200       1, 1   0.3 0.5         0.7
#> 2048         2  85879.905    final        1200       1, 1   0.3 0.5         0.7
#> 2049         2  85947.578    final        1200       1, 1   0.3 0.5         0.7
#> 2050         2  86019.810    final        1200       1, 1   0.3 0.5         0.7
#> 2051         2  86050.390    final        1200       1, 1   0.3 0.5         0.7
#> 2052         2  86094.103    final        1200       1, 1   0.3 0.5         0.7
#> 2053         2  86169.429    final        1200       1, 1   0.3 0.5         0.7
#> 2054         2  86192.936    final        1200       1, 1   0.3 0.5         0.7
#> 2055         2  86304.131    final        1200       1, 1   0.3 0.5         0.7
#> 2056         2  86367.472    final        1200       1, 1   0.3 0.5         0.7
#> 2057         2  86415.444    final        1200       1, 1   0.3 0.5         0.7
#> 2058         2  86463.365    final        1200       1, 1   0.3 0.5         0.7
#> 2059         2  86507.063    final        1200       1, 1   0.3 0.5         0.7
#> 2060         2  86542.261    final        1200       1, 1   0.3 0.5         0.7
#> 2061         2  86776.901    final        1200       1, 1   0.3 0.5         0.7
#> 2062         2  86966.173    final        1200       1, 1   0.3 0.5         0.7
#> 2063         2  87050.785    final        1200       1, 1   0.3 0.5         0.7
#> 2064         2  87066.356    final        1200       1, 1   0.3 0.5         0.7
#> 2065         2  87320.422    final        1200       1, 1   0.3 0.5         0.7
#> 2066         2  87386.949    final        1200       1, 1   0.3 0.5         0.7
#> 2067         2  87404.219    final        1200       1, 1   0.3 0.5         0.7
#> 2068         2  87499.565    final        1200       1, 1   0.3 0.5         0.7
#> 2069         2  87716.977    final        1200       1, 1   0.3 0.5         0.7
#> 2070         2  87913.326    final        1200       1, 1   0.3 0.5         0.7
#> 2071         2  88192.571    final        1200       1, 1   0.3 0.5         0.7
#> 2072         2  88309.850    final        1200       1, 1   0.3 0.5         0.7
#> 2073         2  88377.126    final        1200       1, 1   0.3 0.5         0.7
#> 2074         2  88438.511    final        1200       1, 1   0.3 0.5         0.7
#> 2075         2  88443.004    final        1200       1, 1   0.3 0.5         0.7
#> 2076         2  88481.102    final        1200       1, 1   0.3 0.5         0.7
#> 2077         2  88770.166    final        1200       1, 1   0.3 0.5         0.7
#> 2078         2  88939.400    final        1200       1, 1   0.3 0.5         0.7
#> 2079         2  88957.442    final        1200       1, 1   0.3 0.5         0.7
#> 2080         2  88970.804    final        1200       1, 1   0.3 0.5         0.7
#> 2081         2  89024.875    final        1200       1, 1   0.3 0.5         0.7
#> 2082         2  89104.347    final        1200       1, 1   0.3 0.5         0.7
#> 2083         2  89143.749    final        1200       1, 1   0.3 0.5         0.7
#> 2084         2  89333.907    final        1200       1, 1   0.3 0.5         0.7
#> 2085         2  89625.097    final        1200       1, 1   0.3 0.5         0.7
#> 2086         2  89743.977    final        1200       1, 1   0.3 0.5         0.7
#> 2087         2  89936.239    final        1200       1, 1   0.3 0.5         0.7
#> 2088         2  90017.042    final        1200       1, 1   0.3 0.5         0.7
#> 2089         2  90216.316    final        1200       1, 1   0.3 0.5         0.7
#> 2090         2  90217.269    final        1200       1, 1   0.3 0.5         0.7
#> 2091         2  90238.715    final        1200       1, 1   0.3 0.5         0.7
#> 2092         2  90373.188    final        1200       1, 1   0.3 0.5         0.7
#> 2093         2  90376.250    final        1200       1, 1   0.3 0.5         0.7
#> 2094         2  90440.010    final        1200       1, 1   0.3 0.5         0.7
#> 2095         2  90513.918    final        1200       1, 1   0.3 0.5         0.7
#> 2096         2  90550.500    final        1200       1, 1   0.3 0.5         0.7
#> 2097         2  90759.121    final        1200       1, 1   0.3 0.5         0.7
#> 2098         2  90797.024    final        1200       1, 1   0.3 0.5         0.7
#> 2099         2  91107.144    final        1200       1, 1   0.3 0.5         0.7
#> 2100         2  91191.665    final        1200       1, 1   0.3 0.5         0.7
#> 2101         2  91330.104    final        1200       1, 1   0.3 0.5         0.7
#> 2102         2  91413.372    final        1200       1, 1   0.3 0.5         0.7
#> 2103         2  91491.174    final        1200       1, 1   0.3 0.5         0.7
#> 2104         2  91711.578    final        1200       1, 1   0.3 0.5         0.7
#> 2105         2  91735.636    final        1200       1, 1   0.3 0.5         0.7
#> 2106         2  91759.310    final        1200       1, 1   0.3 0.5         0.7
#> 2107         2  91923.973    final        1200       1, 1   0.3 0.5         0.7
#> 2108         2  92063.227    final        1200       1, 1   0.3 0.5         0.7
#> 2109         2  92080.392    final        1200       1, 1   0.3 0.5         0.7
#> 2110         2  92391.431    final        1200       1, 1   0.3 0.5         0.7
#> 2111         2  92516.915    final        1200       1, 1   0.3 0.5         0.7
#> 2112         2  92757.524    final        1200       1, 1   0.3 0.5         0.7
#> 2113         2  93049.793    final        1200       1, 1   0.3 0.5         0.7
#> 2114         2  93227.213    final        1200       1, 1   0.3 0.5         0.7
#> 2115         2  93353.164    final        1200       1, 1   0.3 0.5         0.7
#> 2116         2  93409.215    final        1200       1, 1   0.3 0.5         0.7
#> 2117         2  93409.979    final        1200       1, 1   0.3 0.5         0.7
#> 2118         2  93607.234    final        1200       1, 1   0.3 0.5         0.7
#> 2119         2  93610.559    final        1200       1, 1   0.3 0.5         0.7
#> 2120         2  93617.555    final        1200       1, 1   0.3 0.5         0.7
#> 2121         2  93712.647    final        1200       1, 1   0.3 0.5         0.7
#> 2122         2  93772.870    final        1200       1, 1   0.3 0.5         0.7
#> 2123         2  93807.432    final        1200       1, 1   0.3 0.5         0.7
#> 2124         2  93939.970    final        1200       1, 1   0.3 0.5         0.7
#> 2125         2  93969.725    final        1200       1, 1   0.3 0.5         0.7
#> 2126         2  93982.593    final        1200       1, 1   0.3 0.5         0.7
#> 2127         2  94002.151    final        1200       1, 1   0.3 0.5         0.7
#> 2128         2  94124.495    final        1200       1, 1   0.3 0.5         0.7
#> 2129         2  94300.120    final        1200       1, 1   0.3 0.5         0.7
#> 2130         2  94481.379    final        1200       1, 1   0.3 0.5         0.7
#> 2131         2  94520.845    final        1200       1, 1   0.3 0.5         0.7
#> 2132         2  94789.490    final        1200       1, 1   0.3 0.5         0.7
#> 2133         2  94935.597    final        1200       1, 1   0.3 0.5         0.7
#> 2134         2  95045.690    final        1200       1, 1   0.3 0.5         0.7
#> 2135         2  95091.285    final        1200       1, 1   0.3 0.5         0.7
#> 2136         2  95126.186    final        1200       1, 1   0.3 0.5         0.7
#> 2137         2  95178.600    final        1200       1, 1   0.3 0.5         0.7
#> 2138         2  95288.707    final        1200       1, 1   0.3 0.5         0.7
#> 2139         2  95888.650    final        1200       1, 1   0.3 0.5         0.7
#> 2140         2  95969.726    final        1200       1, 1   0.3 0.5         0.7
#> 2141         2  96007.987    final        1200       1, 1   0.3 0.5         0.7
#> 2142         2  96105.554    final        1200       1, 1   0.3 0.5         0.7
#> 2143         2  96286.609    final        1200       1, 1   0.3 0.5         0.7
#> 2144         2  96433.312    final        1200       1, 1   0.3 0.5         0.7
#> 2145         2  96573.312    final        1200       1, 1   0.3 0.5         0.7
#> 2146         2  96726.414    final        1200       1, 1   0.3 0.5         0.7
#> 2147         2  96837.279    final        1200       1, 1   0.3 0.5         0.7
#> 2148         2  96933.002    final        1200       1, 1   0.3 0.5         0.7
#> 2149         2  97059.167    final        1200       1, 1   0.3 0.5         0.7
#> 2150         2  97087.277    final        1200       1, 1   0.3 0.5         0.7
#> 2151         2  97094.685    final        1200       1, 1   0.3 0.5         0.7
#> 2152         2  97134.335    final        1200       1, 1   0.3 0.5         0.7
#> 2153         2  97253.328    final        1200       1, 1   0.3 0.5         0.7
#> 2154         2  97322.061    final        1200       1, 1   0.3 0.5         0.7
#> 2155         2  97522.062    final        1200       1, 1   0.3 0.5         0.7
#> 2156         2  97540.494    final        1200       1, 1   0.3 0.5         0.7
#> 2157         2  97548.192    final        1200       1, 1   0.3 0.5         0.7
#> 2158         2  97953.400    final        1200       1, 1   0.3 0.5         0.7
#> 2159         2  98023.112    final        1200       1, 1   0.3 0.5         0.7
#> 2160         2  98035.114    final        1200       1, 1   0.3 0.5         0.7
#> 2161         2  98069.249    final        1200       1, 1   0.3 0.5         0.7
#> 2162         2  98105.873    final        1200       1, 1   0.3 0.5         0.7
#> 2163         2  98186.210    final        1200       1, 1   0.3 0.5         0.7
#> 2164         2  98232.344    final        1200       1, 1   0.3 0.5         0.7
#> 2165         2  98311.546    final        1200       1, 1   0.3 0.5         0.7
#> 2166         2  98386.911    final        1200       1, 1   0.3 0.5         0.7
#> 2167         2  98737.643    final        1200       1, 1   0.3 0.5         0.7
#> 2168         2  98747.819    final        1200       1, 1   0.3 0.5         0.7
#> 2169         2  98813.391    final        1200       1, 1   0.3 0.5         0.7
#> 2170         2  98826.441    final        1200       1, 1   0.3 0.5         0.7
#> 2171         2  99022.330    final        1200       1, 1   0.3 0.5         0.7
#> 2172         2  99075.024    final        1200       1, 1   0.3 0.5         0.7
#> 2173         2  99589.574    final        1200       1, 1   0.3 0.5         0.7
#> 2174         2  99691.185    final        1200       1, 1   0.3 0.5         0.7
#> 2175         2  99716.407    final        1200       1, 1   0.3 0.5         0.7
#> 2176         2  99994.961    final        1200       1, 1   0.3 0.5         0.7
#> 2177         2 100016.038    final        1200       1, 1   0.3 0.5         0.7
#> 2178         2 100166.468    final        1200       1, 1   0.3 0.5         0.7
#> 2179         2 100180.138    final        1200       1, 1   0.3 0.5         0.7
#> 2180         2 100224.981    final        1200       1, 1   0.3 0.5         0.7
#> 2181         2 100293.882    final        1200       1, 1   0.3 0.5         0.7
#> 2182         2 100393.208    final        1200       1, 1   0.3 0.5         0.7
#> 2183         2 100465.576    final        1200       1, 1   0.3 0.5         0.7
#> 2184         2 100552.143    final        1200       1, 1   0.3 0.5         0.7
#> 2185         2 100685.326    final        1200       1, 1   0.3 0.5         0.7
#> 2186         2 100848.159    final        1200       1, 1   0.3 0.5         0.7
#> 2187         2 100938.180    final        1200       1, 1   0.3 0.5         0.7
#> 2188         2 100965.086    final        1200       1, 1   0.3 0.5         0.7
#> 2189         2 101035.470    final        1200       1, 1   0.3 0.5         0.7
#> 2190         2 101082.465    final        1200       1, 1   0.3 0.5         0.7
#> 2191         2 101217.249    final        1200       1, 1   0.3 0.5         0.7
#> 2192         2 101251.331    final        1200       1, 1   0.3 0.5         0.7
#> 2193         2 101314.941    final        1200       1, 1   0.3 0.5         0.7
#> 2194         2 101332.008    final        1200       1, 1   0.3 0.5         0.7
#> 2195         2 101378.201    final        1200       1, 1   0.3 0.5         0.7
#> 2196         2 101459.404    final        1200       1, 1   0.3 0.5         0.7
#> 2197         2 101510.139    final        1200       1, 1   0.3 0.5         0.7
#> 2198         2 101677.832    final        1200       1, 1   0.3 0.5         0.7
#> 2199         2 101693.991    final        1200       1, 1   0.3 0.5         0.7
#> 2200         2 101839.548    final        1200       1, 1   0.3 0.5         0.7
#> 2201         2 101938.547    final        1200       1, 1   0.3 0.5         0.7
#> 2202         2 101949.069    final        1200       1, 1   0.3 0.5         0.7
#> 2203         2 102213.011    final        1200       1, 1   0.3 0.5         0.7
#> 2204         2 102485.973    final        1200       1, 1   0.3 0.5         0.7
#> 2205         2 102520.633    final        1200       1, 1   0.3 0.5         0.7
#> 2206         2 102524.839    final        1200       1, 1   0.3 0.5         0.7
#> 2207         2 102611.726    final        1200       1, 1   0.3 0.5         0.7
#> 2208         2 102617.744    final        1200       1, 1   0.3 0.5         0.7
#> 2209         2 102971.781    final        1200       1, 1   0.3 0.5         0.7
#> 2210         2 102991.339    final        1200       1, 1   0.3 0.5         0.7
#> 2211         2 103102.927    final        1200       1, 1   0.3 0.5         0.7
#> 2212         2 103323.785    final        1200       1, 1   0.3 0.5         0.7
#> 2213         2 103427.084    final        1200       1, 1   0.3 0.5         0.7
#> 2214         2 103651.521    final        1200       1, 1   0.3 0.5         0.7
#> 2215         2 103817.448    final        1200       1, 1   0.3 0.5         0.7
#> 2216         2 103888.569    final        1200       1, 1   0.3 0.5         0.7
#> 2217         2 103915.852    final        1200       1, 1   0.3 0.5         0.7
#> 2218         2 103961.885    final        1200       1, 1   0.3 0.5         0.7
#> 2219         2 104192.573    final        1200       1, 1   0.3 0.5         0.7
#> 2220         2 104316.588    final        1200       1, 1   0.3 0.5         0.7
#> 2221         2 104382.409    final        1200       1, 1   0.3 0.5         0.7
#> 2222         2 104432.513    final        1200       1, 1   0.3 0.5         0.7
#> 2223         2 104475.754    final        1200       1, 1   0.3 0.5         0.7
#> 2224         2 104672.574    final        1200       1, 1   0.3 0.5         0.7
#> 2225         2 104900.123    final        1200       1, 1   0.3 0.5         0.7
#> 2226         2 104903.386    final        1200       1, 1   0.3 0.5         0.7
#> 2227         2 104908.479    final        1200       1, 1   0.3 0.5         0.7
#> 2228         2 105005.002    final        1200       1, 1   0.3 0.5         0.7
#> 2229         2 105011.043    final        1200       1, 1   0.3 0.5         0.7
#> 2230         2 105014.338    final        1200       1, 1   0.3 0.5         0.7
#> 2231         2 105056.447    final        1200       1, 1   0.3 0.5         0.7
#> 2232         2 105127.606    final        1200       1, 1   0.3 0.5         0.7
#> 2233         2 105212.979    final        1200       1, 1   0.3 0.5         0.7
#> 2234         2 105287.998    final        1200       1, 1   0.3 0.5         0.7
#> 2235         2 105303.366    final        1200       1, 1   0.3 0.5         0.7
#> 2236         2 105371.631    final        1200       1, 1   0.3 0.5         0.7
#> 2237         2 105533.836    final        1200       1, 1   0.3 0.5         0.7
#> 2238         2 105550.592    final        1200       1, 1   0.3 0.5         0.7
#> 2239         2 105714.404    final        1200       1, 1   0.3 0.5         0.7
#> 2240         2 105866.876    final        1200       1, 1   0.3 0.5         0.7
#> 2241         2 105932.900    final        1200       1, 1   0.3 0.5         0.7
#> 2242         2 106250.649    final        1200       1, 1   0.3 0.5         0.7
#> 2243         2 106286.183    final        1200       1, 1   0.3 0.5         0.7
#> 2244         2 106329.867    final        1200       1, 1   0.3 0.5         0.7
#> 2245         2 106337.186    final        1200       1, 1   0.3 0.5         0.7
#> 2246         2 106365.240    final        1200       1, 1   0.3 0.5         0.7
#> 2247         2 106436.645    final        1200       1, 1   0.3 0.5         0.7
#> 2248         2 106510.942    final        1200       1, 1   0.3 0.5         0.7
#> 2249         2 106557.062    final        1200       1, 1   0.3 0.5         0.7
#> 2250         2 106581.166    final        1200       1, 1   0.3 0.5         0.7
#> 2251         2 106982.972    final        1200       1, 1   0.3 0.5         0.7
#> 2252         2 107037.957    final        1200       1, 1   0.3 0.5         0.7
#> 2253         2 107122.139    final        1200       1, 1   0.3 0.5         0.7
#> 2254         2 107150.906    final        1200       1, 1   0.3 0.5         0.7
#> 2255         2 107323.880    final        1200       1, 1   0.3 0.5         0.7
#> 2256         2 107352.115    final        1200       1, 1   0.3 0.5         0.7
#> 2257         2 107381.740    final        1200       1, 1   0.3 0.5         0.7
#> 2258         2 107433.342    final        1200       1, 1   0.3 0.5         0.7
#> 2259         2 107451.755    final        1200       1, 1   0.3 0.5         0.7
#> 2260         2 107495.667    final        1200       1, 1   0.3 0.5         0.7
#> 2261         2 107646.998    final        1200       1, 1   0.3 0.5         0.7
#> 2262         2 107662.191    final        1200       1, 1   0.3 0.5         0.7
#> 2263         2 107896.959    final        1200       1, 1   0.3 0.5         0.7
#> 2264         2 107898.633    final        1200       1, 1   0.3 0.5         0.7
#> 2265         2 107937.921    final        1200       1, 1   0.3 0.5         0.7
#> 2266         2 108328.319    final        1200       1, 1   0.3 0.5         0.7
#> 2267         2 108492.921    final        1200       1, 1   0.3 0.5         0.7
#> 2268         2 108740.227    final        1200       1, 1   0.3 0.5         0.7
#> 2269         2 108758.224    final        1200       1, 1   0.3 0.5         0.7
#> 2270         2 108789.438    final        1200       1, 1   0.3 0.5         0.7
#> 2271         2 108954.608    final        1200       1, 1   0.3 0.5         0.7
#> 2272         2 109175.990    final        1200       1, 1   0.3 0.5         0.7
#> 2273         2 109230.554    final        1200       1, 1   0.3 0.5         0.7
#> 2274         2 109348.164    final        1200       1, 1   0.3 0.5         0.7
#> 2275         2 109497.043    final        1200       1, 1   0.3 0.5         0.7
#> 2276         2 109649.839    final        1200       1, 1   0.3 0.5         0.7
#> 2277         2 109689.406    final        1200       1, 1   0.3 0.5         0.7
#> 2278         2 110236.908    final        1200       1, 1   0.3 0.5         0.7
#> 2279         2 110297.144    final        1200       1, 1   0.3 0.5         0.7
#> 2280         2 110539.633    final        1200       1, 1   0.3 0.5         0.7
#> 2281         2 110558.288    final        1200       1, 1   0.3 0.5         0.7
#> 2282         2 110615.801    final        1200       1, 1   0.3 0.5         0.7
#> 2283         2 110815.259    final        1200       1, 1   0.3 0.5         0.7
#> 2284         2 110943.188    final        1200       1, 1   0.3 0.5         0.7
#> 2285         2 111003.280    final        1200       1, 1   0.3 0.5         0.7
#> 2286         2 111014.285    final        1200       1, 1   0.3 0.5         0.7
#> 2287         2 111100.260    final        1200       1, 1   0.3 0.5         0.7
#> 2288         2 111381.505    final        1200       1, 1   0.3 0.5         0.7
#> 2289         2 111442.528    final        1200       1, 1   0.3 0.5         0.7
#> 2290         2 111459.416    final        1200       1, 1   0.3 0.5         0.7
#> 2291         2 111474.099    final        1200       1, 1   0.3 0.5         0.7
#> 2292         2 111610.985    final        1200       1, 1   0.3 0.5         0.7
#> 2293         2 111616.355    final        1200       1, 1   0.3 0.5         0.7
#> 2294         2 111759.337    final        1200       1, 1   0.3 0.5         0.7
#> 2295         2 111772.473    final        1200       1, 1   0.3 0.5         0.7
#> 2296         2 111825.255    final        1200       1, 1   0.3 0.5         0.7
#> 2297         2 112399.245    final        1200       1, 1   0.3 0.5         0.7
#> 2298         2 112494.933    final        1200       1, 1   0.3 0.5         0.7
#> 2299         2 112522.079    final        1200       1, 1   0.3 0.5         0.7
#> 2300         2 112634.848    final        1200       1, 1   0.3 0.5         0.7
#> 2301         2 112675.023    final        1200       1, 1   0.3 0.5         0.7
#> 2302         2 112809.991    final        1200       1, 1   0.3 0.5         0.7
#> 2303         2 112817.691    final        1200       1, 1   0.3 0.5         0.7
#> 2304         2 112839.319    final        1200       1, 1   0.3 0.5         0.7
#> 2305         2 113069.882    final        1200       1, 1   0.3 0.5         0.7
#> 2306         2 113132.574    final        1200       1, 1   0.3 0.5         0.7
#> 2307         2 113292.395    final        1200       1, 1   0.3 0.5         0.7
#> 2308         2 113395.418    final        1200       1, 1   0.3 0.5         0.7
#> 2309         2 113498.478    final        1200       1, 1   0.3 0.5         0.7
#> 2310         2 113575.751    final        1200       1, 1   0.3 0.5         0.7
#> 2311         2 113641.098    final        1200       1, 1   0.3 0.5         0.7
#> 2312         2 114061.466    final        1200       1, 1   0.3 0.5         0.7
#> 2313         2 114109.353    final        1200       1, 1   0.3 0.5         0.7
#> 2314         2 114158.440    final        1200       1, 1   0.3 0.5         0.7
#> 2315         2 114233.399    final        1200       1, 1   0.3 0.5         0.7
#> 2316         2 114247.553    final        1200       1, 1   0.3 0.5         0.7
#> 2317         2 114548.257    final        1200       1, 1   0.3 0.5         0.7
#> 2318         2 114577.193    final        1200       1, 1   0.3 0.5         0.7
#> 2319         2 114644.575    final        1200       1, 1   0.3 0.5         0.7
#> 2320         2 114667.976    final        1200       1, 1   0.3 0.5         0.7
#> 2321         2 114801.181    final        1200       1, 1   0.3 0.5         0.7
#> 2322         2 115068.611    final        1200       1, 1   0.3 0.5         0.7
#> 2323         2 115113.855    final        1200       1, 1   0.3 0.5         0.7
#> 2324         2 115353.988    final        1200       1, 1   0.3 0.5         0.7
#> 2325         2 115365.526    final        1200       1, 1   0.3 0.5         0.7
#> 2326         2 115451.942    final        1200       1, 1   0.3 0.5         0.7
#> 2327         2 115477.539    final        1200       1, 1   0.3 0.5         0.7
#> 2328         2 115497.321    final        1200       1, 1   0.3 0.5         0.7
#> 2329         2 115501.530    final        1200       1, 1   0.3 0.5         0.7
#> 2330         2 115510.273    final        1200       1, 1   0.3 0.5         0.7
#> 2331         2 115565.471    final        1200       1, 1   0.3 0.5         0.7
#> 2332         2 115681.307    final        1200       1, 1   0.3 0.5         0.7
#> 2333         2 115764.761    final        1200       1, 1   0.3 0.5         0.7
#> 2334         2 115770.874    final        1200       1, 1   0.3 0.5         0.7
#> 2335         2 115913.919    final        1200       1, 1   0.3 0.5         0.7
#> 2336         2 116265.815    final        1200       1, 1   0.3 0.5         0.7
#> 2337         2 116372.968    final        1200       1, 1   0.3 0.5         0.7
#> 2338         2 116398.413    final        1200       1, 1   0.3 0.5         0.7
#> 2339         2 116446.811    final        1200       1, 1   0.3 0.5         0.7
#> 2340         2 116482.696    final        1200       1, 1   0.3 0.5         0.7
#> 2341         2 116543.009    final        1200       1, 1   0.3 0.5         0.7
#> 2342         2 116603.674    final        1200       1, 1   0.3 0.5         0.7
#> 2343         2 116818.335    final        1200       1, 1   0.3 0.5         0.7
#> 2344         2 117012.330    final        1200       1, 1   0.3 0.5         0.7
#> 2345         2 117037.238    final        1200       1, 1   0.3 0.5         0.7
#> 2346         2 117163.366    final        1200       1, 1   0.3 0.5         0.7
#> 2347         2 117191.833    final        1200       1, 1   0.3 0.5         0.7
#> 2348         2 117215.838    final        1200       1, 1   0.3 0.5         0.7
#> 2349         2 117219.495    final        1200       1, 1   0.3 0.5         0.7
#> 2350         2 117438.530    final        1200       1, 1   0.3 0.5         0.7
#> 2351         2 117504.979    final        1200       1, 1   0.3 0.5         0.7
#> 2352         2 117594.970    final        1200       1, 1   0.3 0.5         0.7
#> 2353         2 117595.553    final        1200       1, 1   0.3 0.5         0.7
#> 2354         2 117657.998    final        1200       1, 1   0.3 0.5         0.7
#> 2355         2 117704.093    final        1200       1, 1   0.3 0.5         0.7
#> 2356         2 117839.805    final        1200       1, 1   0.3 0.5         0.7
#> 2357         2 117895.988    final        1200       1, 1   0.3 0.5         0.7
#> 2358         2 117991.548    final        1200       1, 1   0.3 0.5         0.7
#> 2359         2 117999.923    final        1200       1, 1   0.3 0.5         0.7
#> 2360         2 118003.438    final        1200       1, 1   0.3 0.5         0.7
#> 2361         2 118108.168    final        1200       1, 1   0.3 0.5         0.7
#> 2362         2 118161.691    final        1200       1, 1   0.3 0.5         0.7
#> 2363         2 118175.365    final        1200       1, 1   0.3 0.5         0.7
#> 2364         2 118234.803    final        1200       1, 1   0.3 0.5         0.7
#> 2365         2 118298.134    final        1200       1, 1   0.3 0.5         0.7
#> 2366         2 118505.340    final        1200       1, 1   0.3 0.5         0.7
#> 2367         2 118521.502    final        1200       1, 1   0.3 0.5         0.7
#> 2368         2 118712.549    final        1200       1, 1   0.3 0.5         0.7
#> 2369         2 118744.323    final        1200       1, 1   0.3 0.5         0.7
#> 2370         2 118823.736    final        1200       1, 1   0.3 0.5         0.7
#> 2371         2 118992.959    final        1200       1, 1   0.3 0.5         0.7
#> 2372         2 119208.347    final        1200       1, 1   0.3 0.5         0.7
#> 2373         2 119402.289    final        1200       1, 1   0.3 0.5         0.7
#> 2374         2 119538.342    final        1200       1, 1   0.3 0.5         0.7
#> 2375         2 119622.502    final        1200       1, 1   0.3 0.5         0.7
#> 2376         2 119629.639    final        1200       1, 1   0.3 0.5         0.7
#> 2377         2 119656.255    final        1200       1, 1   0.3 0.5         0.7
#> 2378         2 119889.850    final        1200       1, 1   0.3 0.5         0.7
#> 2379         2 119935.519    final        1200       1, 1   0.3 0.5         0.7
#> 2380         2 120143.326    final        1200       1, 1   0.3 0.5         0.7
#> 2381         2 120349.838    final        1200       1, 1   0.3 0.5         0.7
#> 2382         2 120363.241    final        1200       1, 1   0.3 0.5         0.7
#> 2383         2 120390.863    final        1200       1, 1   0.3 0.5         0.7
#> 2384         3   1214.035    final        1200       1, 1   0.3 0.5         0.7
#> 2385         3   1360.230    final        1200       1, 1   0.3 0.5         0.7
#> 2386         3   1415.124    final        1200       1, 1   0.3 0.5         0.7
#> 2387         3   1518.131    final        1200       1, 1   0.3 0.5         0.7
#> 2388         3   1543.670    final        1200       1, 1   0.3 0.5         0.7
#> 2389         3   1553.015    final        1200       1, 1   0.3 0.5         0.7
#> 2390         3   1627.097    final        1200       1, 1   0.3 0.5         0.7
#> 2391         3   1729.297    final        1200       1, 1   0.3 0.5         0.7
#> 2392         3   1979.902    final        1200       1, 1   0.3 0.5         0.7
#> 2393         3   2017.107    final        1200       1, 1   0.3 0.5         0.7
#> 2394         3   2111.685    final        1200       1, 1   0.3 0.5         0.7
#> 2395         3   2317.668    final        1200       1, 1   0.3 0.5         0.7
#> 2396         3   2751.713    final        1200       1, 1   0.3 0.5         0.7
#> 2397         3   2988.736    final        1200       1, 1   0.3 0.5         0.7
#> 2398         3   3001.622    final        1200       1, 1   0.3 0.5         0.7
#> 2399         3   3115.149    final        1200       1, 1   0.3 0.5         0.7
#> 2400         3   3127.977    final        1200       1, 1   0.3 0.5         0.7
#> 2401         3   3146.146    final        1200       1, 1   0.3 0.5         0.7
#> 2402         3   3198.608    final        1200       1, 1   0.3 0.5         0.7
#> 2403         3   3272.582    final        1200       1, 1   0.3 0.5         0.7
#> 2404         3   3409.028    final        1200       1, 1   0.3 0.5         0.7
#> 2405         3   3621.014    final        1200       1, 1   0.3 0.5         0.7
#> 2406         3   3626.559    final        1200       1, 1   0.3 0.5         0.7
#> 2407         3   3655.147    final        1200       1, 1   0.3 0.5         0.7
#> 2408         3   3861.313    final        1200       1, 1   0.3 0.5         0.7
#> 2409         3   3869.905    final        1200       1, 1   0.3 0.5         0.7
#> 2410         3   4003.176    final        1200       1, 1   0.3 0.5         0.7
#> 2411         3   4198.162    final        1200       1, 1   0.3 0.5         0.7
#> 2412         3   4304.729    final        1200       1, 1   0.3 0.5         0.7
#> 2413         3   4311.130    final        1200       1, 1   0.3 0.5         0.7
#> 2414         3   4321.882    final        1200       1, 1   0.3 0.5         0.7
#> 2415         3   4386.106    final        1200       1, 1   0.3 0.5         0.7
#> 2416         3   4528.266    final        1200       1, 1   0.3 0.5         0.7
#> 2417         3   4631.038    final        1200       1, 1   0.3 0.5         0.7
#> 2418         3   4653.494    final        1200       1, 1   0.3 0.5         0.7
#> 2419         3   4907.747    final        1200       1, 1   0.3 0.5         0.7
#> 2420         3   4912.262    final        1200       1, 1   0.3 0.5         0.7
#> 2421         3   4921.603    final        1200       1, 1   0.3 0.5         0.7
#> 2422         3   4965.233    final        1200       1, 1   0.3 0.5         0.7
#> 2423         3   5021.577    final        1200       1, 1   0.3 0.5         0.7
#> 2424         3   5081.336    final        1200       1, 1   0.3 0.5         0.7
#> 2425         3   5214.466    final        1200       1, 1   0.3 0.5         0.7
#> 2426         3   5318.106    final        1200       1, 1   0.3 0.5         0.7
#> 2427         3   5354.894    final        1200       1, 1   0.3 0.5         0.7
#> 2428         3   5541.459    final        1200       1, 1   0.3 0.5         0.7
#> 2429         3   5547.134    final        1200       1, 1   0.3 0.5         0.7
#> 2430         3   5619.009    final        1200       1, 1   0.3 0.5         0.7
#> 2431         3   5635.094    final        1200       1, 1   0.3 0.5         0.7
#> 2432         3   5690.101    final        1200       1, 1   0.3 0.5         0.7
#> 2433         3   5881.186    final        1200       1, 1   0.3 0.5         0.7
#> 2434         3   5997.217    final        1200       1, 1   0.3 0.5         0.7
#> 2435         3   6029.201    final        1200       1, 1   0.3 0.5         0.7
#> 2436         3   6136.327    final        1200       1, 1   0.3 0.5         0.7
#> 2437         3   6336.239    final        1200       1, 1   0.3 0.5         0.7
#> 2438         3   6422.074    final        1200       1, 1   0.3 0.5         0.7
#> 2439         3   6448.762    final        1200       1, 1   0.3 0.5         0.7
#> 2440         3   6563.319    final        1200       1, 1   0.3 0.5         0.7
#> 2441         3   6616.239    final        1200       1, 1   0.3 0.5         0.7
#> 2442         3   6642.216    final        1200       1, 1   0.3 0.5         0.7
#> 2443         3   6945.448    final        1200       1, 1   0.3 0.5         0.7
#> 2444         3   7093.226    final        1200       1, 1   0.3 0.5         0.7
#> 2445         3   7532.758    final        1200       1, 1   0.3 0.5         0.7
#> 2446         3   7555.578    final        1200       1, 1   0.3 0.5         0.7
#> 2447         3   7569.601    final        1200       1, 1   0.3 0.5         0.7
#> 2448         3   7579.870    final        1200       1, 1   0.3 0.5         0.7
#> 2449         3   7754.020    final        1200       1, 1   0.3 0.5         0.7
#> 2450         3   7826.951    final        1200       1, 1   0.3 0.5         0.7
#> 2451         3   7869.123    final        1200       1, 1   0.3 0.5         0.7
#> 2452         3   8044.363    final        1200       1, 1   0.3 0.5         0.7
#> 2453         3   8134.710    final        1200       1, 1   0.3 0.5         0.7
#> 2454         3   8143.034    final        1200       1, 1   0.3 0.5         0.7
#> 2455         3   8150.712    final        1200       1, 1   0.3 0.5         0.7
#> 2456         3   8250.471    final        1200       1, 1   0.3 0.5         0.7
#> 2457         3   8322.642    final        1200       1, 1   0.3 0.5         0.7
#> 2458         3   8339.664    final        1200       1, 1   0.3 0.5         0.7
#> 2459         3   8422.954    final        1200       1, 1   0.3 0.5         0.7
#> 2460         3   8455.265    final        1200       1, 1   0.3 0.5         0.7
#> 2461         3   8825.524    final        1200       1, 1   0.3 0.5         0.7
#> 2462         3   8859.882    final        1200       1, 1   0.3 0.5         0.7
#> 2463         3   8882.441    final        1200       1, 1   0.3 0.5         0.7
#> 2464         3   9113.051    final        1200       1, 1   0.3 0.5         0.7
#> 2465         3   9128.056    final        1200       1, 1   0.3 0.5         0.7
#> 2466         3   9382.706    final        1200       1, 1   0.3 0.5         0.7
#> 2467         3   9410.618    final        1200       1, 1   0.3 0.5         0.7
#> 2468         3   9750.775    final        1200       1, 1   0.3 0.5         0.7
#> 2469         3   9759.819    final        1200       1, 1   0.3 0.5         0.7
#> 2470         3   9812.853    final        1200       1, 1   0.3 0.5         0.7
#> 2471         3   9872.369    final        1200       1, 1   0.3 0.5         0.7
#> 2472         3  10038.328    final        1200       1, 1   0.3 0.5         0.7
#> 2473         3  10075.298    final        1200       1, 1   0.3 0.5         0.7
#> 2474         3  10373.364    final        1200       1, 1   0.3 0.5         0.7
#> 2475         3  10448.219    final        1200       1, 1   0.3 0.5         0.7
#> 2476         3  10530.749    final        1200       1, 1   0.3 0.5         0.7
#> 2477         3  10737.300    final        1200       1, 1   0.3 0.5         0.7
#> 2478         3  10780.377    final        1200       1, 1   0.3 0.5         0.7
#> 2479         3  10956.015    final        1200       1, 1   0.3 0.5         0.7
#> 2480         3  11290.543    final        1200       1, 1   0.3 0.5         0.7
#> 2481         3  11372.068    final        1200       1, 1   0.3 0.5         0.7
#> 2482         3  11386.897    final        1200       1, 1   0.3 0.5         0.7
#> 2483         3  11460.017    final        1200       1, 1   0.3 0.5         0.7
#> 2484         3  11619.947    final        1200       1, 1   0.3 0.5         0.7
#> 2485         3  11750.076    final        1200       1, 1   0.3 0.5         0.7
#> 2486         3  11790.996    final        1200       1, 1   0.3 0.5         0.7
#> 2487         3  11854.500    final        1200       1, 1   0.3 0.5         0.7
#> 2488         3  11920.329    final        1200       1, 1   0.3 0.5         0.7
#> 2489         3  11923.684    final        1200       1, 1   0.3 0.5         0.7
#> 2490         3  11965.751    final        1200       1, 1   0.3 0.5         0.7
#> 2491         3  12100.060    final        1200       1, 1   0.3 0.5         0.7
#> 2492         3  12130.683    final        1200       1, 1   0.3 0.5         0.7
#> 2493         3  12158.682    final        1200       1, 1   0.3 0.5         0.7
#> 2494         3  12160.923    final        1200       1, 1   0.3 0.5         0.7
#> 2495         3  12162.111    final        1200       1, 1   0.3 0.5         0.7
#> 2496         3  12275.247    final        1200       1, 1   0.3 0.5         0.7
#> 2497         3  12418.149    final        1200       1, 1   0.3 0.5         0.7
#> 2498         3  12501.147    final        1200       1, 1   0.3 0.5         0.7
#> 2499         3  12502.431    final        1200       1, 1   0.3 0.5         0.7
#> 2500         3  12598.347    final        1200       1, 1   0.3 0.5         0.7
#> 2501         3  12656.555    final        1200       1, 1   0.3 0.5         0.7
#> 2502         3  13328.499    final        1200       1, 1   0.3 0.5         0.7
#> 2503         3  13367.528    final        1200       1, 1   0.3 0.5         0.7
#> 2504         3  13555.102    final        1200       1, 1   0.3 0.5         0.7
#> 2505         3  13599.552    final        1200       1, 1   0.3 0.5         0.7
#> 2506         3  13669.731    final        1200       1, 1   0.3 0.5         0.7
#> 2507         3  13671.096    final        1200       1, 1   0.3 0.5         0.7
#> 2508         3  13730.162    final        1200       1, 1   0.3 0.5         0.7
#> 2509         3  13767.226    final        1200       1, 1   0.3 0.5         0.7
#> 2510         3  13778.210    final        1200       1, 1   0.3 0.5         0.7
#> 2511         3  13788.305    final        1200       1, 1   0.3 0.5         0.7
#> 2512         3  13840.415    final        1200       1, 1   0.3 0.5         0.7
#> 2513         3  13902.567    final        1200       1, 1   0.3 0.5         0.7
#> 2514         3  13917.500    final        1200       1, 1   0.3 0.5         0.7
#> 2515         3  14091.987    final        1200       1, 1   0.3 0.5         0.7
#> 2516         3  14199.872    final        1200       1, 1   0.3 0.5         0.7
#> 2517         3  14328.105    final        1200       1, 1   0.3 0.5         0.7
#> 2518         3  14373.364    final        1200       1, 1   0.3 0.5         0.7
#> 2519         3  14518.942    final        1200       1, 1   0.3 0.5         0.7
#> 2520         3  14552.928    final        1200       1, 1   0.3 0.5         0.7
#> 2521         3  14650.393    final        1200       1, 1   0.3 0.5         0.7
#> 2522         3  14655.248    final        1200       1, 1   0.3 0.5         0.7
#> 2523         3  14703.789    final        1200       1, 1   0.3 0.5         0.7
#> 2524         3  14737.383    final        1200       1, 1   0.3 0.5         0.7
#> 2525         3  14744.899    final        1200       1, 1   0.3 0.5         0.7
#> 2526         3  14766.661    final        1200       1, 1   0.3 0.5         0.7
#> 2527         3  14866.049    final        1200       1, 1   0.3 0.5         0.7
#> 2528         3  14886.481    final        1200       1, 1   0.3 0.5         0.7
#> 2529         3  14952.683    final        1200       1, 1   0.3 0.5         0.7
#> 2530         3  15238.094    final        1200       1, 1   0.3 0.5         0.7
#> 2531         3  15321.254    final        1200       1, 1   0.3 0.5         0.7
#> 2532         3  15422.639    final        1200       1, 1   0.3 0.5         0.7
#> 2533         3  15658.882    final        1200       1, 1   0.3 0.5         0.7
#> 2534         3  15663.171    final        1200       1, 1   0.3 0.5         0.7
#> 2535         3  15911.679    final        1200       1, 1   0.3 0.5         0.7
#> 2536         3  15960.283    final        1200       1, 1   0.3 0.5         0.7
#> 2537         3  15971.265    final        1200       1, 1   0.3 0.5         0.7
#> 2538         3  16067.618    final        1200       1, 1   0.3 0.5         0.7
#> 2539         3  16167.913    final        1200       1, 1   0.3 0.5         0.7
#> 2540         3  16279.147    final        1200       1, 1   0.3 0.5         0.7
#> 2541         3  16436.734    final        1200       1, 1   0.3 0.5         0.7
#> 2542         3  16456.482    final        1200       1, 1   0.3 0.5         0.7
#> 2543         3  16559.569    final        1200       1, 1   0.3 0.5         0.7
#> 2544         3  16574.052    final        1200       1, 1   0.3 0.5         0.7
#> 2545         3  16641.430    final        1200       1, 1   0.3 0.5         0.7
#> 2546         3  16696.861    final        1200       1, 1   0.3 0.5         0.7
#> 2547         3  16793.492    final        1200       1, 1   0.3 0.5         0.7
#> 2548         3  16905.478    final        1200       1, 1   0.3 0.5         0.7
#> 2549         3  17166.337    final        1200       1, 1   0.3 0.5         0.7
#> 2550         3  17227.864    final        1200       1, 1   0.3 0.5         0.7
#> 2551         3  17418.329    final        1200       1, 1   0.3 0.5         0.7
#> 2552         3  17421.696    final        1200       1, 1   0.3 0.5         0.7
#> 2553         3  17587.668    final        1200       1, 1   0.3 0.5         0.7
#> 2554         3  17612.478    final        1200       1, 1   0.3 0.5         0.7
#> 2555         3  17739.734    final        1200       1, 1   0.3 0.5         0.7
#> 2556         3  17839.592    final        1200       1, 1   0.3 0.5         0.7
#> 2557         3  17854.959    final        1200       1, 1   0.3 0.5         0.7
#> 2558         3  18175.119    final        1200       1, 1   0.3 0.5         0.7
#> 2559         3  18228.092    final        1200       1, 1   0.3 0.5         0.7
#> 2560         3  18289.540    final        1200       1, 1   0.3 0.5         0.7
#> 2561         3  18303.559    final        1200       1, 1   0.3 0.5         0.7
#> 2562         3  18451.514    final        1200       1, 1   0.3 0.5         0.7
#> 2563         3  18586.281    final        1200       1, 1   0.3 0.5         0.7
#> 2564         3  18737.753    final        1200       1, 1   0.3 0.5         0.7
#> 2565         3  19224.055    final        1200       1, 1   0.3 0.5         0.7
#> 2566         3  19315.316    final        1200       1, 1   0.3 0.5         0.7
#> 2567         3  19364.623    final        1200       1, 1   0.3 0.5         0.7
#> 2568         3  19441.530    final        1200       1, 1   0.3 0.5         0.7
#> 2569         3  19645.558    final        1200       1, 1   0.3 0.5         0.7
#> 2570         3  19811.408    final        1200       1, 1   0.3 0.5         0.7
#> 2571         3  19867.120    final        1200       1, 1   0.3 0.5         0.7
#> 2572         3  19928.066    final        1200       1, 1   0.3 0.5         0.7
#> 2573         3  19963.583    final        1200       1, 1   0.3 0.5         0.7
#> 2574         3  19967.227    final        1200       1, 1   0.3 0.5         0.7
#> 2575         3  19979.011    final        1200       1, 1   0.3 0.5         0.7
#> 2576         3  20023.644    final        1200       1, 1   0.3 0.5         0.7
#> 2577         3  20051.676    final        1200       1, 1   0.3 0.5         0.7
#> 2578         3  20076.752    final        1200       1, 1   0.3 0.5         0.7
#> 2579         3  20152.247    final        1200       1, 1   0.3 0.5         0.7
#> 2580         3  20167.529    final        1200       1, 1   0.3 0.5         0.7
#> 2581         3  20199.529    final        1200       1, 1   0.3 0.5         0.7
#> 2582         3  20330.364    final        1200       1, 1   0.3 0.5         0.7
#> 2583         3  20360.967    final        1200       1, 1   0.3 0.5         0.7
#> 2584         3  20635.105    final        1200       1, 1   0.3 0.5         0.7
#> 2585         3  20678.656    final        1200       1, 1   0.3 0.5         0.7
#> 2586         3  20783.317    final        1200       1, 1   0.3 0.5         0.7
#> 2587         3  20902.123    final        1200       1, 1   0.3 0.5         0.7
#> 2588         3  21052.557    final        1200       1, 1   0.3 0.5         0.7
#> 2589         3  21060.097    final        1200       1, 1   0.3 0.5         0.7
#> 2590         3  21081.100    final        1200       1, 1   0.3 0.5         0.7
#> 2591         3  21511.339    final        1200       1, 1   0.3 0.5         0.7
#> 2592         3  21612.209    final        1200       1, 1   0.3 0.5         0.7
#> 2593         3  21786.664    final        1200       1, 1   0.3 0.5         0.7
#> 2594         3  21849.011    final        1200       1, 1   0.3 0.5         0.7
#> 2595         3  21954.174    final        1200       1, 1   0.3 0.5         0.7
#> 2596         3  21999.330    final        1200       1, 1   0.3 0.5         0.7
#> 2597         3  22085.296    final        1200       1, 1   0.3 0.5         0.7
#> 2598         3  22141.599    final        1200       1, 1   0.3 0.5         0.7
#> 2599         3  22249.541    final        1200       1, 1   0.3 0.5         0.7
#> 2600         3  22299.637    final        1200       1, 1   0.3 0.5         0.7
#> 2601         3  22343.384    final        1200       1, 1   0.3 0.5         0.7
#> 2602         3  22387.172    final        1200       1, 1   0.3 0.5         0.7
#> 2603         3  22401.036    final        1200       1, 1   0.3 0.5         0.7
#> 2604         3  22484.996    final        1200       1, 1   0.3 0.5         0.7
#> 2605         3  22499.223    final        1200       1, 1   0.3 0.5         0.7
#> 2606         3  22600.914    final        1200       1, 1   0.3 0.5         0.7
#> 2607         3  22696.580    final        1200       1, 1   0.3 0.5         0.7
#> 2608         3  22701.167    final        1200       1, 1   0.3 0.5         0.7
#> 2609         3  22953.534    final        1200       1, 1   0.3 0.5         0.7
#> 2610         3  22981.943    final        1200       1, 1   0.3 0.5         0.7
#> 2611         3  23029.210    final        1200       1, 1   0.3 0.5         0.7
#> 2612         3  23101.012    final        1200       1, 1   0.3 0.5         0.7
#> 2613         3  23111.663    final        1200       1, 1   0.3 0.5         0.7
#> 2614         3  23135.182    final        1200       1, 1   0.3 0.5         0.7
#> 2615         3  23440.156    final        1200       1, 1   0.3 0.5         0.7
#> 2616         3  23501.716    final        1200       1, 1   0.3 0.5         0.7
#> 2617         3  23827.387    final        1200       1, 1   0.3 0.5         0.7
#> 2618         3  24033.142    final        1200       1, 1   0.3 0.5         0.7
#> 2619         3  24177.175    final        1200       1, 1   0.3 0.5         0.7
#> 2620         3  24435.367    final        1200       1, 1   0.3 0.5         0.7
#> 2621         3  24547.121    final        1200       1, 1   0.3 0.5         0.7
#> 2622         3  24629.981    final        1200       1, 1   0.3 0.5         0.7
#> 2623         3  24708.938    final        1200       1, 1   0.3 0.5         0.7
#> 2624         3  24865.995    final        1200       1, 1   0.3 0.5         0.7
#> 2625         3  24887.638    final        1200       1, 1   0.3 0.5         0.7
#> 2626         3  24887.710    final        1200       1, 1   0.3 0.5         0.7
#> 2627         3  24993.022    final        1200       1, 1   0.3 0.5         0.7
#> 2628         3  25084.381    final        1200       1, 1   0.3 0.5         0.7
#> 2629         3  25158.614    final        1200       1, 1   0.3 0.5         0.7
#> 2630         3  25164.032    final        1200       1, 1   0.3 0.5         0.7
#> 2631         3  25228.185    final        1200       1, 1   0.3 0.5         0.7
#> 2632         3  25319.907    final        1200       1, 1   0.3 0.5         0.7
#> 2633         3  25356.624    final        1200       1, 1   0.3 0.5         0.7
#> 2634         3  25415.565    final        1200       1, 1   0.3 0.5         0.7
#> 2635         3  25431.243    final        1200       1, 1   0.3 0.5         0.7
#> 2636         3  25503.356    final        1200       1, 1   0.3 0.5         0.7
#> 2637         3  25718.827    final        1200       1, 1   0.3 0.5         0.7
#> 2638         3  25789.852    final        1200       1, 1   0.3 0.5         0.7
#> 2639         3  25790.197    final        1200       1, 1   0.3 0.5         0.7
#> 2640         3  25802.977    final        1200       1, 1   0.3 0.5         0.7
#> 2641         3  25888.519    final        1200       1, 1   0.3 0.5         0.7
#> 2642         3  25968.239    final        1200       1, 1   0.3 0.5         0.7
#> 2643         3  26294.328    final        1200       1, 1   0.3 0.5         0.7
#> 2644         3  26314.078    final        1200       1, 1   0.3 0.5         0.7
#> 2645         3  26558.496    final        1200       1, 1   0.3 0.5         0.7
#> 2646         3  26583.134    final        1200       1, 1   0.3 0.5         0.7
#> 2647         3  26602.006    final        1200       1, 1   0.3 0.5         0.7
#> 2648         3  26605.117    final        1200       1, 1   0.3 0.5         0.7
#> 2649         3  26616.891    final        1200       1, 1   0.3 0.5         0.7
#> 2650         3  26734.267    final        1200       1, 1   0.3 0.5         0.7
#> 2651         3  26814.475    final        1200       1, 1   0.3 0.5         0.7
#> 2652         3  26833.729    final        1200       1, 1   0.3 0.5         0.7
#> 2653         3  26914.927    final        1200       1, 1   0.3 0.5         0.7
#> 2654         3  26938.140    final        1200       1, 1   0.3 0.5         0.7
#> 2655         3  26970.135    final        1200       1, 1   0.3 0.5         0.7
#> 2656         3  27038.599    final        1200       1, 1   0.3 0.5         0.7
#> 2657         3  27068.641    final        1200       1, 1   0.3 0.5         0.7
#> 2658         3  27205.045    final        1200       1, 1   0.3 0.5         0.7
#> 2659         3  27334.583    final        1200       1, 1   0.3 0.5         0.7
#> 2660         3  27362.749    final        1200       1, 1   0.3 0.5         0.7
#> 2661         3  27582.612    final        1200       1, 1   0.3 0.5         0.7
#> 2662         3  27632.565    final        1200       1, 1   0.3 0.5         0.7
#> 2663         3  27767.990    final        1200       1, 1   0.3 0.5         0.7
#> 2664         3  27778.597    final        1200       1, 1   0.3 0.5         0.7
#> 2665         3  27855.593    final        1200       1, 1   0.3 0.5         0.7
#> 2666         3  28022.895    final        1200       1, 1   0.3 0.5         0.7
#> 2667         3  28113.427    final        1200       1, 1   0.3 0.5         0.7
#> 2668         3  28135.482    final        1200       1, 1   0.3 0.5         0.7
#> 2669         3  28203.359    final        1200       1, 1   0.3 0.5         0.7
#> 2670         3  28423.864    final        1200       1, 1   0.3 0.5         0.7
#> 2671         3  28957.800    final        1200       1, 1   0.3 0.5         0.7
#> 2672         3  29063.755    final        1200       1, 1   0.3 0.5         0.7
#> 2673         3  29667.923    final        1200       1, 1   0.3 0.5         0.7
#> 2674         3  29835.820    final        1200       1, 1   0.3 0.5         0.7
#> 2675         3  29954.452    final        1200       1, 1   0.3 0.5         0.7
#> 2676         3  30029.726    final        1200       1, 1   0.3 0.5         0.7
#> 2677         3  30034.395    final        1200       1, 1   0.3 0.5         0.7
#> 2678         3  30064.173    final        1200       1, 1   0.3 0.5         0.7
#> 2679         3  30071.567    final        1200       1, 1   0.3 0.5         0.7
#> 2680         3  30172.675    final        1200       1, 1   0.3 0.5         0.7
#> 2681         3  30285.342    final        1200       1, 1   0.3 0.5         0.7
#> 2682         3  30407.762    final        1200       1, 1   0.3 0.5         0.7
#> 2683         3  30426.936    final        1200       1, 1   0.3 0.5         0.7
#> 2684         3  30465.584    final        1200       1, 1   0.3 0.5         0.7
#> 2685         3  30782.902    final        1200       1, 1   0.3 0.5         0.7
#> 2686         3  30849.682    final        1200       1, 1   0.3 0.5         0.7
#> 2687         3  30867.089    final        1200       1, 1   0.3 0.5         0.7
#> 2688         3  31202.131    final        1200       1, 1   0.3 0.5         0.7
#> 2689         3  31210.906    final        1200       1, 1   0.3 0.5         0.7
#> 2690         3  31470.743    final        1200       1, 1   0.3 0.5         0.7
#> 2691         3  31479.150    final        1200       1, 1   0.3 0.5         0.7
#> 2692         3  31482.456    final        1200       1, 1   0.3 0.5         0.7
#> 2693         3  31593.764    final        1200       1, 1   0.3 0.5         0.7
#> 2694         3  31851.531    final        1200       1, 1   0.3 0.5         0.7
#> 2695         3  32097.275    final        1200       1, 1   0.3 0.5         0.7
#> 2696         3  32158.043    final        1200       1, 1   0.3 0.5         0.7
#> 2697         3  32239.281    final        1200       1, 1   0.3 0.5         0.7
#> 2698         3  32549.582    final        1200       1, 1   0.3 0.5         0.7
#> 2699         3  32678.252    final        1200       1, 1   0.3 0.5         0.7
#> 2700         3  32724.088    final        1200       1, 1   0.3 0.5         0.7
#> 2701         3  32813.206    final        1200       1, 1   0.3 0.5         0.7
#> 2702         3  32882.866    final        1200       1, 1   0.3 0.5         0.7
#> 2703         3  32951.605    final        1200       1, 1   0.3 0.5         0.7
#> 2704         3  33027.924    final        1200       1, 1   0.3 0.5         0.7
#> 2705         3  33085.140    final        1200       1, 1   0.3 0.5         0.7
#> 2706         3  33145.194    final        1200       1, 1   0.3 0.5         0.7
#> 2707         3  33358.437    final        1200       1, 1   0.3 0.5         0.7
#> 2708         3  33459.129    final        1200       1, 1   0.3 0.5         0.7
#> 2709         3  33518.482    final        1200       1, 1   0.3 0.5         0.7
#> 2710         3  33661.757    final        1200       1, 1   0.3 0.5         0.7
#> 2711         3  33668.095    final        1200       1, 1   0.3 0.5         0.7
#> 2712         3  33673.423    final        1200       1, 1   0.3 0.5         0.7
#> 2713         3  33868.544    final        1200       1, 1   0.3 0.5         0.7
#> 2714         3  33894.488    final        1200       1, 1   0.3 0.5         0.7
#> 2715         3  33977.656    final        1200       1, 1   0.3 0.5         0.7
#> 2716         3  34012.612    final        1200       1, 1   0.3 0.5         0.7
#> 2717         3  34055.485    final        1200       1, 1   0.3 0.5         0.7
#> 2718         3  34275.433    final        1200       1, 1   0.3 0.5         0.7
#> 2719         3  34298.485    final        1200       1, 1   0.3 0.5         0.7
#> 2720         3  34315.722    final        1200       1, 1   0.3 0.5         0.7
#> 2721         3  34318.527    final        1200       1, 1   0.3 0.5         0.7
#> 2722         3  34423.148    final        1200       1, 1   0.3 0.5         0.7
#> 2723         3  34644.258    final        1200       1, 1   0.3 0.5         0.7
#> 2724         3  34711.329    final        1200       1, 1   0.3 0.5         0.7
#> 2725         3  34742.601    final        1200       1, 1   0.3 0.5         0.7
#> 2726         3  34744.157    final        1200       1, 1   0.3 0.5         0.7
#> 2727         3  34940.077    final        1200       1, 1   0.3 0.5         0.7
#> 2728         3  35060.834    final        1200       1, 1   0.3 0.5         0.7
#> 2729         3  35130.666    final        1200       1, 1   0.3 0.5         0.7
#> 2730         3  35197.562    final        1200       1, 1   0.3 0.5         0.7
#> 2731         3  35422.299    final        1200       1, 1   0.3 0.5         0.7
#> 2732         3  35440.172    final        1200       1, 1   0.3 0.5         0.7
#> 2733         3  35460.523    final        1200       1, 1   0.3 0.5         0.7
#> 2734         3  35580.658    final        1200       1, 1   0.3 0.5         0.7
#> 2735         3  35693.436    final        1200       1, 1   0.3 0.5         0.7
#> 2736         3  36130.430    final        1200       1, 1   0.3 0.5         0.7
#> 2737         3  36257.046    final        1200       1, 1   0.3 0.5         0.7
#> 2738         3  36364.923    final        1200       1, 1   0.3 0.5         0.7
#> 2739         3  36404.386    final        1200       1, 1   0.3 0.5         0.7
#> 2740         3  36827.394    final        1200       1, 1   0.3 0.5         0.7
#> 2741         3  36869.447    final        1200       1, 1   0.3 0.5         0.7
#> 2742         3  36883.326    final        1200       1, 1   0.3 0.5         0.7
#> 2743         3  37020.782    final        1200       1, 1   0.3 0.5         0.7
#> 2744         3  37129.902    final        1200       1, 1   0.3 0.5         0.7
#> 2745         3  37211.590    final        1200       1, 1   0.3 0.5         0.7
#> 2746         3  37251.599    final        1200       1, 1   0.3 0.5         0.7
#> 2747         3  37397.625    final        1200       1, 1   0.3 0.5         0.7
#> 2748         3  37471.814    final        1200       1, 1   0.3 0.5         0.7
#> 2749         3  37813.264    final        1200       1, 1   0.3 0.5         0.7
#> 2750         3  38048.889    final        1200       1, 1   0.3 0.5         0.7
#> 2751         3  38193.095    final        1200       1, 1   0.3 0.5         0.7
#> 2752         3  38201.230    final        1200       1, 1   0.3 0.5         0.7
#> 2753         3  38212.238    final        1200       1, 1   0.3 0.5         0.7
#> 2754         3  38341.696    final        1200       1, 1   0.3 0.5         0.7
#> 2755         3  38383.476    final        1200       1, 1   0.3 0.5         0.7
#> 2756         3  38395.344    final        1200       1, 1   0.3 0.5         0.7
#> 2757         3  38434.889    final        1200       1, 1   0.3 0.5         0.7
#> 2758         3  38463.416    final        1200       1, 1   0.3 0.5         0.7
#> 2759         3  38476.895    final        1200       1, 1   0.3 0.5         0.7
#> 2760         3  38510.438    final        1200       1, 1   0.3 0.5         0.7
#> 2761         3  38514.744    final        1200       1, 1   0.3 0.5         0.7
#> 2762         3  38726.632    final        1200       1, 1   0.3 0.5         0.7
#> 2763         3  38905.074    final        1200       1, 1   0.3 0.5         0.7
#> 2764         3  38916.055    final        1200       1, 1   0.3 0.5         0.7
#> 2765         3  38920.293    final        1200       1, 1   0.3 0.5         0.7
#> 2766         3  38972.983    final        1200       1, 1   0.3 0.5         0.7
#> 2767         3  39117.057    final        1200       1, 1   0.3 0.5         0.7
#> 2768         3  39135.728    final        1200       1, 1   0.3 0.5         0.7
#> 2769         3  39178.040    final        1200       1, 1   0.3 0.5         0.7
#> 2770         3  39201.298    final        1200       1, 1   0.3 0.5         0.7
#> 2771         3  39287.200    final        1200       1, 1   0.3 0.5         0.7
#> 2772         3  39379.801    final        1200       1, 1   0.3 0.5         0.7
#> 2773         3  39389.932    final        1200       1, 1   0.3 0.5         0.7
#> 2774         3  39634.415    final        1200       1, 1   0.3 0.5         0.7
#> 2775         3  39790.766    final        1200       1, 1   0.3 0.5         0.7
#> 2776         3  39979.159    final        1200       1, 1   0.3 0.5         0.7
#> 2777         3  40300.715    final        1200       1, 1   0.3 0.5         0.7
#> 2778         3  40424.097    final        1200       1, 1   0.3 0.5         0.7
#> 2779         3  40433.103    final        1200       1, 1   0.3 0.5         0.7
#> 2780         3  40554.306    final        1200       1, 1   0.3 0.5         0.7
#> 2781         3  40569.594    final        1200       1, 1   0.3 0.5         0.7
#> 2782         3  40733.466    final        1200       1, 1   0.3 0.5         0.7
#> 2783         3  40742.557    final        1200       1, 1   0.3 0.5         0.7
#> 2784         3  41017.740    final        1200       1, 1   0.3 0.5         0.7
#> 2785         3  41018.297    final        1200       1, 1   0.3 0.5         0.7
#> 2786         3  41426.265    final        1200       1, 1   0.3 0.5         0.7
#> 2787         3  41496.441    final        1200       1, 1   0.3 0.5         0.7
#> 2788         3  41538.750    final        1200       1, 1   0.3 0.5         0.7
#> 2789         3  41555.443    final        1200       1, 1   0.3 0.5         0.7
#> 2790         3  41595.385    final        1200       1, 1   0.3 0.5         0.7
#> 2791         3  41736.501    final        1200       1, 1   0.3 0.5         0.7
#> 2792         3  42043.546    final        1200       1, 1   0.3 0.5         0.7
#> 2793         3  42141.590    final        1200       1, 1   0.3 0.5         0.7
#> 2794         3  42227.438    final        1200       1, 1   0.3 0.5         0.7
#> 2795         3  42318.838    final        1200       1, 1   0.3 0.5         0.7
#> 2796         3  42648.813    final        1200       1, 1   0.3 0.5         0.7
#> 2797         3  42675.755    final        1200       1, 1   0.3 0.5         0.7
#> 2798         3  42772.905    final        1200       1, 1   0.3 0.5         0.7
#> 2799         3  42820.249    final        1200       1, 1   0.3 0.5         0.7
#> 2800         3  42919.244    final        1200       1, 1   0.3 0.5         0.7
#> 2801         3  42963.110    final        1200       1, 1   0.3 0.5         0.7
#> 2802         3  43055.852    final        1200       1, 1   0.3 0.5         0.7
#> 2803         3  43164.317    final        1200       1, 1   0.3 0.5         0.7
#> 2804         3  43348.265    final        1200       1, 1   0.3 0.5         0.7
#> 2805         3  43401.564    final        1200       1, 1   0.3 0.5         0.7
#> 2806         3  43444.052    final        1200       1, 1   0.3 0.5         0.7
#> 2807         3  43493.130    final        1200       1, 1   0.3 0.5         0.7
#> 2808         3  43735.364    final        1200       1, 1   0.3 0.5         0.7
#> 2809         3  43771.366    final        1200       1, 1   0.3 0.5         0.7
#> 2810         3  43954.568    final        1200       1, 1   0.3 0.5         0.7
#> 2811         3  44114.459    final        1200       1, 1   0.3 0.5         0.7
#> 2812         3  44165.982    final        1200       1, 1   0.3 0.5         0.7
#> 2813         3  44319.721    final        1200       1, 1   0.3 0.5         0.7
#> 2814         3  44480.180    final        1200       1, 1   0.3 0.5         0.7
#> 2815         3  44558.896    final        1200       1, 1   0.3 0.5         0.7
#> 2816         3  44633.104    final        1200       1, 1   0.3 0.5         0.7
#> 2817         3  45088.667    final        1200       1, 1   0.3 0.5         0.7
#> 2818         3  45172.706    final        1200       1, 1   0.3 0.5         0.7
#> 2819         3  45203.348    final        1200       1, 1   0.3 0.5         0.7
#> 2820         3  45338.775    final        1200       1, 1   0.3 0.5         0.7
#> 2821         3  45415.097    final        1200       1, 1   0.3 0.5         0.7
#> 2822         3  45473.942    final        1200       1, 1   0.3 0.5         0.7
#> 2823         3  45532.450    final        1200       1, 1   0.3 0.5         0.7
#> 2824         3  45578.973    final        1200       1, 1   0.3 0.5         0.7
#> 2825         3  45665.904    final        1200       1, 1   0.3 0.5         0.7
#> 2826         3  45674.578    final        1200       1, 1   0.3 0.5         0.7
#> 2827         3  45735.115    final        1200       1, 1   0.3 0.5         0.7
#> 2828         3  45738.104    final        1200       1, 1   0.3 0.5         0.7
#> 2829         3  45922.846    final        1200       1, 1   0.3 0.5         0.7
#> 2830         3  45932.852    final        1200       1, 1   0.3 0.5         0.7
#> 2831         3  46061.120    final        1200       1, 1   0.3 0.5         0.7
#> 2832         3  46111.709    final        1200       1, 1   0.3 0.5         0.7
#> 2833         3  46170.082    final        1200       1, 1   0.3 0.5         0.7
#> 2834         3  46178.672    final        1200       1, 1   0.3 0.5         0.7
#> 2835         3  46313.230    final        1200       1, 1   0.3 0.5         0.7
#> 2836         3  46330.706    final        1200       1, 1   0.3 0.5         0.7
#> 2837         3  46429.299    final        1200       1, 1   0.3 0.5         0.7
#> 2838         3  46584.831    final        1200       1, 1   0.3 0.5         0.7
#> 2839         3  46598.363    final        1200       1, 1   0.3 0.5         0.7
#> 2840         3  46689.670    final        1200       1, 1   0.3 0.5         0.7
#> 2841         3  46842.303    final        1200       1, 1   0.3 0.5         0.7
#> 2842         3  46882.105    final        1200       1, 1   0.3 0.5         0.7
#> 2843         3  46922.078    final        1200       1, 1   0.3 0.5         0.7
#> 2844         3  46952.812    final        1200       1, 1   0.3 0.5         0.7
#> 2845         3  47374.951    final        1200       1, 1   0.3 0.5         0.7
#> 2846         3  47500.916    final        1200       1, 1   0.3 0.5         0.7
#> 2847         3  47618.630    final        1200       1, 1   0.3 0.5         0.7
#> 2848         3  47624.106    final        1200       1, 1   0.3 0.5         0.7
#> 2849         3  47702.224    final        1200       1, 1   0.3 0.5         0.7
#> 2850         3  47887.457    final        1200       1, 1   0.3 0.5         0.7
#> 2851         3  48029.673    final        1200       1, 1   0.3 0.5         0.7
#> 2852         3  48147.988    final        1200       1, 1   0.3 0.5         0.7
#> 2853         3  48206.509    final        1200       1, 1   0.3 0.5         0.7
#> 2854         3  48221.183    final        1200       1, 1   0.3 0.5         0.7
#> 2855         3  48326.865    final        1200       1, 1   0.3 0.5         0.7
#> 2856         3  48536.999    final        1200       1, 1   0.3 0.5         0.7
#> 2857         3  48552.524    final        1200       1, 1   0.3 0.5         0.7
#> 2858         3  48583.442    final        1200       1, 1   0.3 0.5         0.7
#> 2859         3  48799.917    final        1200       1, 1   0.3 0.5         0.7
#> 2860         3  49022.904    final        1200       1, 1   0.3 0.5         0.7
#> 2861         3  49053.838    final        1200       1, 1   0.3 0.5         0.7
#> 2862         3  49186.469    final        1200       1, 1   0.3 0.5         0.7
#> 2863         3  49308.263    final        1200       1, 1   0.3 0.5         0.7
#> 2864         3  49374.829    final        1200       1, 1   0.3 0.5         0.7
#> 2865         3  49595.320    final        1200       1, 1   0.3 0.5         0.7
#> 2866         3  49738.198    final        1200       1, 1   0.3 0.5         0.7
#> 2867         3  50000.884    final        1200       1, 1   0.3 0.5         0.7
#> 2868         3  50091.960    final        1200       1, 1   0.3 0.5         0.7
#> 2869         3  50210.183    final        1200       1, 1   0.3 0.5         0.7
#> 2870         3  50218.092    final        1200       1, 1   0.3 0.5         0.7
#> 2871         3  50298.457    final        1200       1, 1   0.3 0.5         0.7
#> 2872         3  50313.906    final        1200       1, 1   0.3 0.5         0.7
#> 2873         3  50398.708    final        1200       1, 1   0.3 0.5         0.7
#> 2874         3  50448.626    final        1200       1, 1   0.3 0.5         0.7
#> 2875         3  50468.637    final        1200       1, 1   0.3 0.5         0.7
#> 2876         3  50712.446    final        1200       1, 1   0.3 0.5         0.7
#> 2877         3  50716.211    final        1200       1, 1   0.3 0.5         0.7
#> 2878         3  50753.693    final        1200       1, 1   0.3 0.5         0.7
#> 2879         3  50776.757    final        1200       1, 1   0.3 0.5         0.7
#> 2880         3  50842.784    final        1200       1, 1   0.3 0.5         0.7
#> 2881         3  50949.296    final        1200       1, 1   0.3 0.5         0.7
#> 2882         3  50963.204    final        1200       1, 1   0.3 0.5         0.7
#> 2883         3  51078.815    final        1200       1, 1   0.3 0.5         0.7
#> 2884         3  51157.268    final        1200       1, 1   0.3 0.5         0.7
#> 2885         3  51256.996    final        1200       1, 1   0.3 0.5         0.7
#> 2886         3  51296.708    final        1200       1, 1   0.3 0.5         0.7
#> 2887         3  51384.710    final        1200       1, 1   0.3 0.5         0.7
#> 2888         3  51427.110    final        1200       1, 1   0.3 0.5         0.7
#> 2889         3  51564.143    final        1200       1, 1   0.3 0.5         0.7
#> 2890         3  51634.550    final        1200       1, 1   0.3 0.5         0.7
#> 2891         3  51731.469    final        1200       1, 1   0.3 0.5         0.7
#> 2892         3  51770.821    final        1200       1, 1   0.3 0.5         0.7
#> 2893         3  51870.789    final        1200       1, 1   0.3 0.5         0.7
#> 2894         3  51920.675    final        1200       1, 1   0.3 0.5         0.7
#> 2895         3  52038.248    final        1200       1, 1   0.3 0.5         0.7
#> 2896         3  52199.559    final        1200       1, 1   0.3 0.5         0.7
#> 2897         3  52208.617    final        1200       1, 1   0.3 0.5         0.7
#> 2898         3  52389.414    final        1200       1, 1   0.3 0.5         0.7
#> 2899         3  52416.365    final        1200       1, 1   0.3 0.5         0.7
#> 2900         3  52464.549    final        1200       1, 1   0.3 0.5         0.7
#> 2901         3  52513.976    final        1200       1, 1   0.3 0.5         0.7
#> 2902         3  52745.822    final        1200       1, 1   0.3 0.5         0.7
#> 2903         3  52852.771    final        1200       1, 1   0.3 0.5         0.7
#> 2904         3  52944.573    final        1200       1, 1   0.3 0.5         0.7
#> 2905         3  52988.985    final        1200       1, 1   0.3 0.5         0.7
#> 2906         3  53040.321    final        1200       1, 1   0.3 0.5         0.7
#> 2907         3  53053.483    final        1200       1, 1   0.3 0.5         0.7
#> 2908         3  53214.925    final        1200       1, 1   0.3 0.5         0.7
#> 2909         3  53259.587    final        1200       1, 1   0.3 0.5         0.7
#> 2910         3  53326.161    final        1200       1, 1   0.3 0.5         0.7
#> 2911         3  53360.172    final        1200       1, 1   0.3 0.5         0.7
#> 2912         3  53522.832    final        1200       1, 1   0.3 0.5         0.7
#> 2913         3  53588.384    final        1200       1, 1   0.3 0.5         0.7
#> 2914         3  53820.471    final        1200       1, 1   0.3 0.5         0.7
#> 2915         3  53859.862    final        1200       1, 1   0.3 0.5         0.7
#> 2916         3  53869.069    final        1200       1, 1   0.3 0.5         0.7
#> 2917         3  53934.472    final        1200       1, 1   0.3 0.5         0.7
#> 2918         3  54040.996    final        1200       1, 1   0.3 0.5         0.7
#> 2919         3  54413.963    final        1200       1, 1   0.3 0.5         0.7
#> 2920         3  54438.946    final        1200       1, 1   0.3 0.5         0.7
#> 2921         3  54594.716    final        1200       1, 1   0.3 0.5         0.7
#> 2922         3  54731.504    final        1200       1, 1   0.3 0.5         0.7
#> 2923         3  54798.720    final        1200       1, 1   0.3 0.5         0.7
#> 2924         3  54854.955    final        1200       1, 1   0.3 0.5         0.7
#> 2925         3  54943.378    final        1200       1, 1   0.3 0.5         0.7
#> 2926         3  55040.247    final        1200       1, 1   0.3 0.5         0.7
#> 2927         3  55057.592    final        1200       1, 1   0.3 0.5         0.7
#> 2928         3  55143.033    final        1200       1, 1   0.3 0.5         0.7
#> 2929         3  55370.246    final        1200       1, 1   0.3 0.5         0.7
#> 2930         3  55409.373    final        1200       1, 1   0.3 0.5         0.7
#> 2931         3  55441.421    final        1200       1, 1   0.3 0.5         0.7
#> 2932         3  55663.124    final        1200       1, 1   0.3 0.5         0.7
#> 2933         3  55832.458    final        1200       1, 1   0.3 0.5         0.7
#> 2934         3  55848.393    final        1200       1, 1   0.3 0.5         0.7
#> 2935         3  55921.734    final        1200       1, 1   0.3 0.5         0.7
#> 2936         3  55996.503    final        1200       1, 1   0.3 0.5         0.7
#> 2937         3  55998.021    final        1200       1, 1   0.3 0.5         0.7
#> 2938         3  56034.443    final        1200       1, 1   0.3 0.5         0.7
#> 2939         3  56051.899    final        1200       1, 1   0.3 0.5         0.7
#> 2940         3  56055.920    final        1200       1, 1   0.3 0.5         0.7
#> 2941         3  56130.505    final        1200       1, 1   0.3 0.5         0.7
#> 2942         3  56142.680    final        1200       1, 1   0.3 0.5         0.7
#> 2943         3  56289.061    final        1200       1, 1   0.3 0.5         0.7
#> 2944         3  56319.931    final        1200       1, 1   0.3 0.5         0.7
#> 2945         3  56325.714    final        1200       1, 1   0.3 0.5         0.7
#> 2946         3  56365.758    final        1200       1, 1   0.3 0.5         0.7
#> 2947         3  56457.507    final        1200       1, 1   0.3 0.5         0.7
#> 2948         3  56542.102    final        1200       1, 1   0.3 0.5         0.7
#> 2949         3  56622.557    final        1200       1, 1   0.3 0.5         0.7
#> 2950         3  56782.834    final        1200       1, 1   0.3 0.5         0.7
#> 2951         3  56790.592    final        1200       1, 1   0.3 0.5         0.7
#> 2952         3  56813.648    final        1200       1, 1   0.3 0.5         0.7
#> 2953         3  56876.515    final        1200       1, 1   0.3 0.5         0.7
#> 2954         3  56961.848    final        1200       1, 1   0.3 0.5         0.7
#> 2955         3  57102.794    final        1200       1, 1   0.3 0.5         0.7
#> 2956         3  57157.064    final        1200       1, 1   0.3 0.5         0.7
#> 2957         3  57184.348    final        1200       1, 1   0.3 0.5         0.7
#> 2958         3  57568.381    final        1200       1, 1   0.3 0.5         0.7
#> 2959         3  57683.956    final        1200       1, 1   0.3 0.5         0.7
#> 2960         3  57730.460    final        1200       1, 1   0.3 0.5         0.7
#> 2961         3  57898.426    final        1200       1, 1   0.3 0.5         0.7
#> 2962         3  57976.609    final        1200       1, 1   0.3 0.5         0.7
#> 2963         3  58032.367    final        1200       1, 1   0.3 0.5         0.7
#> 2964         3  58100.530    final        1200       1, 1   0.3 0.5         0.7
#> 2965         3  58431.444    final        1200       1, 1   0.3 0.5         0.7
#> 2966         3  58434.088    final        1200       1, 1   0.3 0.5         0.7
#> 2967         3  58450.795    final        1200       1, 1   0.3 0.5         0.7
#> 2968         3  58479.020    final        1200       1, 1   0.3 0.5         0.7
#> 2969         3  58648.084    final        1200       1, 1   0.3 0.5         0.7
#> 2970         3  58768.408    final        1200       1, 1   0.3 0.5         0.7
#> 2971         3  58838.225    final        1200       1, 1   0.3 0.5         0.7
#> 2972         3  58956.507    final        1200       1, 1   0.3 0.5         0.7
#> 2973         3  59109.684    final        1200       1, 1   0.3 0.5         0.7
#> 2974         3  59131.047    final        1200       1, 1   0.3 0.5         0.7
#> 2975         3  59397.449    final        1200       1, 1   0.3 0.5         0.7
#> 2976         3  59461.331    final        1200       1, 1   0.3 0.5         0.7
#> 2977         3  59469.270    final        1200       1, 1   0.3 0.5         0.7
#> 2978         3  59635.543    final        1200       1, 1   0.3 0.5         0.7
#> 2979         3  59696.358    final        1200       1, 1   0.3 0.5         0.7
#> 2980         3  60224.831    final        1200       1, 1   0.3 0.5         0.7
#> 2981         3  60334.408    final        1200       1, 1   0.3 0.5         0.7
#> 2982         3  60366.454    final        1200       1, 1   0.3 0.5         0.7
#> 2983         3  60513.639    final        1200       1, 1   0.3 0.5         0.7
#> 2984         3  60684.264    final        1200       1, 1   0.3 0.5         0.7
#> 2985         3  60717.440    final        1200       1, 1   0.3 0.5         0.7
#> 2986         3  60740.770    final        1200       1, 1   0.3 0.5         0.7
#> 2987         3  60891.175    final        1200       1, 1   0.3 0.5         0.7
#> 2988         3  61053.915    final        1200       1, 1   0.3 0.5         0.7
#> 2989         3  61130.447    final        1200       1, 1   0.3 0.5         0.7
#> 2990         3  61160.274    final        1200       1, 1   0.3 0.5         0.7
#> 2991         3  61196.199    final        1200       1, 1   0.3 0.5         0.7
#> 2992         3  61323.180    final        1200       1, 1   0.3 0.5         0.7
#> 2993         3  61395.036    final        1200       1, 1   0.3 0.5         0.7
#> 2994         3  61443.471    final        1200       1, 1   0.3 0.5         0.7
#> 2995         3  61732.205    final        1200       1, 1   0.3 0.5         0.7
#> 2996         3  61898.907    final        1200       1, 1   0.3 0.5         0.7
#> 2997         3  61913.738    final        1200       1, 1   0.3 0.5         0.7
#> 2998         3  62012.213    final        1200       1, 1   0.3 0.5         0.7
#> 2999         3  62050.187    final        1200       1, 1   0.3 0.5         0.7
#> 3000         3  62120.006    final        1200       1, 1   0.3 0.5         0.7
#> 3001         3  62216.368    final        1200       1, 1   0.3 0.5         0.7
#> 3002         3  62376.743    final        1200       1, 1   0.3 0.5         0.7
#> 3003         3  62380.488    final        1200       1, 1   0.3 0.5         0.7
#> 3004         3  62391.256    final        1200       1, 1   0.3 0.5         0.7
#> 3005         3  62421.765    final        1200       1, 1   0.3 0.5         0.7
#> 3006         3  62700.311    final        1200       1, 1   0.3 0.5         0.7
#> 3007         3  62916.514    final        1200       1, 1   0.3 0.5         0.7
#> 3008         3  63174.195    final        1200       1, 1   0.3 0.5         0.7
#> 3009         3  63221.040    final        1200       1, 1   0.3 0.5         0.7
#> 3010         3  63409.114    final        1200       1, 1   0.3 0.5         0.7
#> 3011         3  63430.830    final        1200       1, 1   0.3 0.5         0.7
#> 3012         3  63445.174    final        1200       1, 1   0.3 0.5         0.7
#> 3013         3  63521.581    final        1200       1, 1   0.3 0.5         0.7
#> 3014         3  63586.994    final        1200       1, 1   0.3 0.5         0.7
#> 3015         3  63725.209    final        1200       1, 1   0.3 0.5         0.7
#> 3016         3  63772.874    final        1200       1, 1   0.3 0.5         0.7
#> 3017         3  63774.205    final        1200       1, 1   0.3 0.5         0.7
#> 3018         3  63925.127    final        1200       1, 1   0.3 0.5         0.7
#> 3019         3  63928.922    final        1200       1, 1   0.3 0.5         0.7
#> 3020         3  63932.426    final        1200       1, 1   0.3 0.5         0.7
#> 3021         3  63983.238    final        1200       1, 1   0.3 0.5         0.7
#> 3022         3  64203.189    final        1200       1, 1   0.3 0.5         0.7
#> 3023         3  64276.060    final        1200       1, 1   0.3 0.5         0.7
#> 3024         3  64476.641    final        1200       1, 1   0.3 0.5         0.7
#> 3025         3  64616.777    final        1200       1, 1   0.3 0.5         0.7
#> 3026         3  64718.571    final        1200       1, 1   0.3 0.5         0.7
#> 3027         3  64901.755    final        1200       1, 1   0.3 0.5         0.7
#> 3028         3  64990.549    final        1200       1, 1   0.3 0.5         0.7
#> 3029         3  65020.402    final        1200       1, 1   0.3 0.5         0.7
#> 3030         3  65023.044    final        1200       1, 1   0.3 0.5         0.7
#> 3031         3  65144.001    final        1200       1, 1   0.3 0.5         0.7
#> 3032         3  65169.546    final        1200       1, 1   0.3 0.5         0.7
#> 3033         3  65235.550    final        1200       1, 1   0.3 0.5         0.7
#> 3034         3  65598.359    final        1200       1, 1   0.3 0.5         0.7
#> 3035         3  65808.759    final        1200       1, 1   0.3 0.5         0.7
#> 3036         3  65882.755    final        1200       1, 1   0.3 0.5         0.7
#> 3037         3  65908.296    final        1200       1, 1   0.3 0.5         0.7
#> 3038         3  65909.427    final        1200       1, 1   0.3 0.5         0.7
#> 3039         3  65971.682    final        1200       1, 1   0.3 0.5         0.7
#> 3040         3  65989.712    final        1200       1, 1   0.3 0.5         0.7
#> 3041         3  66046.897    final        1200       1, 1   0.3 0.5         0.7
#> 3042         3  66115.860    final        1200       1, 1   0.3 0.5         0.7
#> 3043         3  66185.456    final        1200       1, 1   0.3 0.5         0.7
#> 3044         3  66235.416    final        1200       1, 1   0.3 0.5         0.7
#> 3045         3  66366.714    final        1200       1, 1   0.3 0.5         0.7
#> 3046         3  66382.812    final        1200       1, 1   0.3 0.5         0.7
#> 3047         3  66387.624    final        1200       1, 1   0.3 0.5         0.7
#> 3048         3  66568.419    final        1200       1, 1   0.3 0.5         0.7
#> 3049         3  66573.515    final        1200       1, 1   0.3 0.5         0.7
#> 3050         3  66647.457    final        1200       1, 1   0.3 0.5         0.7
#> 3051         3  66653.714    final        1200       1, 1   0.3 0.5         0.7
#> 3052         3  66661.483    final        1200       1, 1   0.3 0.5         0.7
#> 3053         3  66825.277    final        1200       1, 1   0.3 0.5         0.7
#> 3054         3  66926.831    final        1200       1, 1   0.3 0.5         0.7
#> 3055         3  67008.971    final        1200       1, 1   0.3 0.5         0.7
#> 3056         3  67196.177    final        1200       1, 1   0.3 0.5         0.7
#> 3057         3  67290.573    final        1200       1, 1   0.3 0.5         0.7
#> 3058         3  67346.039    final        1200       1, 1   0.3 0.5         0.7
#> 3059         3  67404.527    final        1200       1, 1   0.3 0.5         0.7
#> 3060         3  67846.689    final        1200       1, 1   0.3 0.5         0.7
#> 3061         3  67895.610    final        1200       1, 1   0.3 0.5         0.7
#> 3062         3  68040.008    final        1200       1, 1   0.3 0.5         0.7
#> 3063         3  68115.483    final        1200       1, 1   0.3 0.5         0.7
#> 3064         3  68136.213    final        1200       1, 1   0.3 0.5         0.7
#> 3065         3  68187.213    final        1200       1, 1   0.3 0.5         0.7
#> 3066         3  68202.315    final        1200       1, 1   0.3 0.5         0.7
#> 3067         3  68259.286    final        1200       1, 1   0.3 0.5         0.7
#> 3068         3  68450.100    final        1200       1, 1   0.3 0.5         0.7
#> 3069         3  68464.687    final        1200       1, 1   0.3 0.5         0.7
#> 3070         3  68885.994    final        1200       1, 1   0.3 0.5         0.7
#> 3071         3  68967.579    final        1200       1, 1   0.3 0.5         0.7
#> 3072         3  69012.612    final        1200       1, 1   0.3 0.5         0.7
#> 3073         3  69024.584    final        1200       1, 1   0.3 0.5         0.7
#> 3074         3  69228.633    final        1200       1, 1   0.3 0.5         0.7
#> 3075         3  69261.196    final        1200       1, 1   0.3 0.5         0.7
#> 3076         3  69552.741    final        1200       1, 1   0.3 0.5         0.7
#> 3077         3  69614.865    final        1200       1, 1   0.3 0.5         0.7
#> 3078         3  69628.425    final        1200       1, 1   0.3 0.5         0.7
#> 3079         3  69728.982    final        1200       1, 1   0.3 0.5         0.7
#> 3080         3  70456.202    final        1200       1, 1   0.3 0.5         0.7
#> 3081         3  70519.786    final        1200       1, 1   0.3 0.5         0.7
#> 3082         3  70594.851    final        1200       1, 1   0.3 0.5         0.7
#> 3083         3  70677.005    final        1200       1, 1   0.3 0.5         0.7
#> 3084         3  70883.522    final        1200       1, 1   0.3 0.5         0.7
#> 3085         3  71033.201    final        1200       1, 1   0.3 0.5         0.7
#> 3086         3  71076.748    final        1200       1, 1   0.3 0.5         0.7
#> 3087         3  71188.248    final        1200       1, 1   0.3 0.5         0.7
#> 3088         3  71283.908    final        1200       1, 1   0.3 0.5         0.7
#> 3089         3  71361.337    final        1200       1, 1   0.3 0.5         0.7
#> 3090         3  71396.594    final        1200       1, 1   0.3 0.5         0.7
#> 3091         3  71424.596    final        1200       1, 1   0.3 0.5         0.7
#> 3092         3  71623.146    final        1200       1, 1   0.3 0.5         0.7
#> 3093         3  71709.874    final        1200       1, 1   0.3 0.5         0.7
#> 3094         3  71723.585    final        1200       1, 1   0.3 0.5         0.7
#> 3095         3  71732.728    final        1200       1, 1   0.3 0.5         0.7
#> 3096         3  71757.985    final        1200       1, 1   0.3 0.5         0.7
#> 3097         3  71764.670    final        1200       1, 1   0.3 0.5         0.7
#> 3098         3  71880.227    final        1200       1, 1   0.3 0.5         0.7
#> 3099         3  72172.934    final        1200       1, 1   0.3 0.5         0.7
#> 3100         3  72179.008    final        1200       1, 1   0.3 0.5         0.7
#> 3101         3  72192.210    final        1200       1, 1   0.3 0.5         0.7
#> 3102         3  72336.088    final        1200       1, 1   0.3 0.5         0.7
#> 3103         3  72352.421    final        1200       1, 1   0.3 0.5         0.7
#> 3104         3  72864.226    final        1200       1, 1   0.3 0.5         0.7
#> 3105         3  72996.274    final        1200       1, 1   0.3 0.5         0.7
#> 3106         3  73339.441    final        1200       1, 1   0.3 0.5         0.7
#> 3107         3  73721.059    final        1200       1, 1   0.3 0.5         0.7
#> 3108         3  73863.422    final        1200       1, 1   0.3 0.5         0.7
#> 3109         3  73871.168    final        1200       1, 1   0.3 0.5         0.7
#> 3110         3  74167.633    final        1200       1, 1   0.3 0.5         0.7
#> 3111         3  74483.560    final        1200       1, 1   0.3 0.5         0.7
#> 3112         3  74623.548    final        1200       1, 1   0.3 0.5         0.7
#> 3113         3  74716.084    final        1200       1, 1   0.3 0.5         0.7
#> 3114         3  74745.332    final        1200       1, 1   0.3 0.5         0.7
#> 3115         3  74747.768    final        1200       1, 1   0.3 0.5         0.7
#> 3116         3  74981.301    final        1200       1, 1   0.3 0.5         0.7
#> 3117         3  74993.280    final        1200       1, 1   0.3 0.5         0.7
#> 3118         3  75293.717    final        1200       1, 1   0.3 0.5         0.7
#> 3119         3  75294.222    final        1200       1, 1   0.3 0.5         0.7
#> 3120         3  75396.627    final        1200       1, 1   0.3 0.5         0.7
#> 3121         3  75424.798    final        1200       1, 1   0.3 0.5         0.7
#> 3122         3  75467.076    final        1200       1, 1   0.3 0.5         0.7
#> 3123         3  75719.801    final        1200       1, 1   0.3 0.5         0.7
#> 3124         3  75770.753    final        1200       1, 1   0.3 0.5         0.7
#> 3125         3  75811.931    final        1200       1, 1   0.3 0.5         0.7
#> 3126         3  75877.874    final        1200       1, 1   0.3 0.5         0.7
#> 3127         3  75956.812    final        1200       1, 1   0.3 0.5         0.7
#> 3128         3  75965.279    final        1200       1, 1   0.3 0.5         0.7
#> 3129         3  76028.871    final        1200       1, 1   0.3 0.5         0.7
#> 3130         3  76199.494    final        1200       1, 1   0.3 0.5         0.7
#> 3131         3  76403.454    final        1200       1, 1   0.3 0.5         0.7
#> 3132         3  76850.971    final        1200       1, 1   0.3 0.5         0.7
#> 3133         3  76852.008    final        1200       1, 1   0.3 0.5         0.7
#> 3134         3  76915.468    final        1200       1, 1   0.3 0.5         0.7
#> 3135         3  77347.145    final        1200       1, 1   0.3 0.5         0.7
#> 3136         3  77470.242    final        1200       1, 1   0.3 0.5         0.7
#> 3137         3  77637.767    final        1200       1, 1   0.3 0.5         0.7
#> 3138         3  77648.630    final        1200       1, 1   0.3 0.5         0.7
#> 3139         3  77777.042    final        1200       1, 1   0.3 0.5         0.7
#> 3140         3  78010.964    final        1200       1, 1   0.3 0.5         0.7
#> 3141         3  78181.993    final        1200       1, 1   0.3 0.5         0.7
#> 3142         3  78193.512    final        1200       1, 1   0.3 0.5         0.7
#> 3143         3  78255.188    final        1200       1, 1   0.3 0.5         0.7
#> 3144         3  78331.042    final        1200       1, 1   0.3 0.5         0.7
#> 3145         3  78623.596    final        1200       1, 1   0.3 0.5         0.7
#> 3146         3  78672.740    final        1200       1, 1   0.3 0.5         0.7
#> 3147         3  79065.831    final        1200       1, 1   0.3 0.5         0.7
#> 3148         3  79090.095    final        1200       1, 1   0.3 0.5         0.7
#> 3149         3  79172.478    final        1200       1, 1   0.3 0.5         0.7
#> 3150         3  79577.481    final        1200       1, 1   0.3 0.5         0.7
#> 3151         3  79692.069    final        1200       1, 1   0.3 0.5         0.7
#> 3152         3  79828.833    final        1200       1, 1   0.3 0.5         0.7
#> 3153         3  79846.354    final        1200       1, 1   0.3 0.5         0.7
#> 3154         3  80007.299    final        1200       1, 1   0.3 0.5         0.7
#> 3155         3  80010.482    final        1200       1, 1   0.3 0.5         0.7
#> 3156         3  80012.779    final        1200       1, 1   0.3 0.5         0.7
#> 3157         3  80053.638    final        1200       1, 1   0.3 0.5         0.7
#> 3158         3  80075.191    final        1200       1, 1   0.3 0.5         0.7
#> 3159         3  80159.628    final        1200       1, 1   0.3 0.5         0.7
#> 3160         3  80182.122    final        1200       1, 1   0.3 0.5         0.7
#> 3161         3  80297.766    final        1200       1, 1   0.3 0.5         0.7
#> 3162         3  80336.635    final        1200       1, 1   0.3 0.5         0.7
#> 3163         3  80364.892    final        1200       1, 1   0.3 0.5         0.7
#> 3164         3  80499.460    final        1200       1, 1   0.3 0.5         0.7
#> 3165         3  80536.698    final        1200       1, 1   0.3 0.5         0.7
#> 3166         3  80600.567    final        1200       1, 1   0.3 0.5         0.7
#> 3167         3  80660.256    final        1200       1, 1   0.3 0.5         0.7
#> 3168         3  80685.296    final        1200       1, 1   0.3 0.5         0.7
#> 3169         3  80702.768    final        1200       1, 1   0.3 0.5         0.7
#> 3170         3  80734.137    final        1200       1, 1   0.3 0.5         0.7
#> 3171         3  80840.768    final        1200       1, 1   0.3 0.5         0.7
#> 3172         3  80860.075    final        1200       1, 1   0.3 0.5         0.7
#> 3173         3  80895.450    final        1200       1, 1   0.3 0.5         0.7
#> 3174         3  80971.168    final        1200       1, 1   0.3 0.5         0.7
#> 3175         3  81096.929    final        1200       1, 1   0.3 0.5         0.7
#> 3176         3  81214.173    final        1200       1, 1   0.3 0.5         0.7
#> 3177         3  81234.504    final        1200       1, 1   0.3 0.5         0.7
#> 3178         3  81454.090    final        1200       1, 1   0.3 0.5         0.7
#> 3179         3  81458.437    final        1200       1, 1   0.3 0.5         0.7
#> 3180         3  81466.723    final        1200       1, 1   0.3 0.5         0.7
#> 3181         3  81469.426    final        1200       1, 1   0.3 0.5         0.7
#> 3182         3  81556.411    final        1200       1, 1   0.3 0.5         0.7
#> 3183         3  81952.595    final        1200       1, 1   0.3 0.5         0.7
#> 3184         3  81972.900    final        1200       1, 1   0.3 0.5         0.7
#> 3185         3  82105.292    final        1200       1, 1   0.3 0.5         0.7
#> 3186         3  82189.580    final        1200       1, 1   0.3 0.5         0.7
#> 3187         3  82222.585    final        1200       1, 1   0.3 0.5         0.7
#> 3188         3  82356.870    final        1200       1, 1   0.3 0.5         0.7
#> 3189         3  82508.818    final        1200       1, 1   0.3 0.5         0.7
#> 3190         3  82561.482    final        1200       1, 1   0.3 0.5         0.7
#> 3191         3  82620.616    final        1200       1, 1   0.3 0.5         0.7
#> 3192         3  82901.224    final        1200       1, 1   0.3 0.5         0.7
#> 3193         3  83661.992    final        1200       1, 1   0.3 0.5         0.7
#> 3194         3  83691.526    final        1200       1, 1   0.3 0.5         0.7
#> 3195         3  83723.197    final        1200       1, 1   0.3 0.5         0.7
#> 3196         3  83780.946    final        1200       1, 1   0.3 0.5         0.7
#> 3197         3  83913.416    final        1200       1, 1   0.3 0.5         0.7
#> 3198         3  83935.301    final        1200       1, 1   0.3 0.5         0.7
#> 3199         3  84104.746    final        1200       1, 1   0.3 0.5         0.7
#> 3200         3  84106.839    final        1200       1, 1   0.3 0.5         0.7
#> 3201         3  84393.476    final        1200       1, 1   0.3 0.5         0.7
#> 3202         3  84401.539    final        1200       1, 1   0.3 0.5         0.7
#> 3203         3  84470.707    final        1200       1, 1   0.3 0.5         0.7
#> 3204         3  84478.614    final        1200       1, 1   0.3 0.5         0.7
#> 3205         3  84507.712    final        1200       1, 1   0.3 0.5         0.7
#> 3206         3  84583.604    final        1200       1, 1   0.3 0.5         0.7
#> 3207         3  84597.827    final        1200       1, 1   0.3 0.5         0.7
#> 3208         3  84784.255    final        1200       1, 1   0.3 0.5         0.7
#> 3209         3  84818.876    final        1200       1, 1   0.3 0.5         0.7
#> 3210         3  84820.672    final        1200       1, 1   0.3 0.5         0.7
#> 3211         3  84920.552    final        1200       1, 1   0.3 0.5         0.7
#> 3212         3  84993.164    final        1200       1, 1   0.3 0.5         0.7
#> 3213         3  85153.464    final        1200       1, 1   0.3 0.5         0.7
#> 3214         3  85176.892    final        1200       1, 1   0.3 0.5         0.7
#> 3215         3  85244.197    final        1200       1, 1   0.3 0.5         0.7
#> 3216         3  85365.409    final        1200       1, 1   0.3 0.5         0.7
#> 3217         3  85391.442    final        1200       1, 1   0.3 0.5         0.7
#> 3218         3  85401.391    final        1200       1, 1   0.3 0.5         0.7
#> 3219         3  85571.235    final        1200       1, 1   0.3 0.5         0.7
#> 3220         3  85617.119    final        1200       1, 1   0.3 0.5         0.7
#> 3221         3  85761.081    final        1200       1, 1   0.3 0.5         0.7
#> 3222         3  85838.537    final        1200       1, 1   0.3 0.5         0.7
#> 3223         3  85885.906    final        1200       1, 1   0.3 0.5         0.7
#> 3224         3  85965.234    final        1200       1, 1   0.3 0.5         0.7
#> 3225         3  85984.875    final        1200       1, 1   0.3 0.5         0.7
#> 3226         3  85988.883    final        1200       1, 1   0.3 0.5         0.7
#> 3227         3  85991.366    final        1200       1, 1   0.3 0.5         0.7
#> 3228         3  86215.619    final        1200       1, 1   0.3 0.5         0.7
#> 3229         3  86239.827    final        1200       1, 1   0.3 0.5         0.7
#> 3230         3  86258.140    final        1200       1, 1   0.3 0.5         0.7
#> 3231         3  86336.579    final        1200       1, 1   0.3 0.5         0.7
#> 3232         3  86537.106    final        1200       1, 1   0.3 0.5         0.7
#> 3233         3  86560.023    final        1200       1, 1   0.3 0.5         0.7
#> 3234         3  86569.418    final        1200       1, 1   0.3 0.5         0.7
#> 3235         3  86575.622    final        1200       1, 1   0.3 0.5         0.7
#> 3236         3  86581.715    final        1200       1, 1   0.3 0.5         0.7
#> 3237         3  86603.584    final        1200       1, 1   0.3 0.5         0.7
#> 3238         3  86758.118    final        1200       1, 1   0.3 0.5         0.7
#> 3239         3  86877.916    final        1200       1, 1   0.3 0.5         0.7
#> 3240         3  87154.963    final        1200       1, 1   0.3 0.5         0.7
#> 3241         3  87361.811    final        1200       1, 1   0.3 0.5         0.7
#> 3242         3  87402.465    final        1200       1, 1   0.3 0.5         0.7
#> 3243         3  87512.249    final        1200       1, 1   0.3 0.5         0.7
#> 3244         3  87610.853    final        1200       1, 1   0.3 0.5         0.7
#> 3245         3  87626.128    final        1200       1, 1   0.3 0.5         0.7
#> 3246         3  87643.470    final        1200       1, 1   0.3 0.5         0.7
#> 3247         3  87648.024    final        1200       1, 1   0.3 0.5         0.7
#> 3248         3  87656.284    final        1200       1, 1   0.3 0.5         0.7
#> 3249         3  87824.092    final        1200       1, 1   0.3 0.5         0.7
#> 3250         3  87940.257    final        1200       1, 1   0.3 0.5         0.7
#> 3251         3  87978.366    final        1200       1, 1   0.3 0.5         0.7
#> 3252         3  88009.081    final        1200       1, 1   0.3 0.5         0.7
#> 3253         3  88082.969    final        1200       1, 1   0.3 0.5         0.7
#> 3254         3  88096.295    final        1200       1, 1   0.3 0.5         0.7
#> 3255         3  88291.414    final        1200       1, 1   0.3 0.5         0.7
#> 3256         3  88420.686    final        1200       1, 1   0.3 0.5         0.7
#> 3257         3  88437.463    final        1200       1, 1   0.3 0.5         0.7
#> 3258         3  88483.352    final        1200       1, 1   0.3 0.5         0.7
#> 3259         3  88487.824    final        1200       1, 1   0.3 0.5         0.7
#> 3260         3  88796.864    final        1200       1, 1   0.3 0.5         0.7
#> 3261         3  88846.091    final        1200       1, 1   0.3 0.5         0.7
#> 3262         3  88865.815    final        1200       1, 1   0.3 0.5         0.7
#> 3263         3  89143.339    final        1200       1, 1   0.3 0.5         0.7
#> 3264         3  89208.204    final        1200       1, 1   0.3 0.5         0.7
#> 3265         3  89260.742    final        1200       1, 1   0.3 0.5         0.7
#> 3266         3  89309.109    final        1200       1, 1   0.3 0.5         0.7
#> 3267         3  89386.131    final        1200       1, 1   0.3 0.5         0.7
#> 3268         3  89478.150    final        1200       1, 1   0.3 0.5         0.7
#> 3269         3  89680.054    final        1200       1, 1   0.3 0.5         0.7
#> 3270         3  89764.628    final        1200       1, 1   0.3 0.5         0.7
#> 3271         3  89787.290    final        1200       1, 1   0.3 0.5         0.7
#> 3272         3  89787.436    final        1200       1, 1   0.3 0.5         0.7
#> 3273         3  89867.016    final        1200       1, 1   0.3 0.5         0.7
#> 3274         3  89910.117    final        1200       1, 1   0.3 0.5         0.7
#> 3275         3  89911.158    final        1200       1, 1   0.3 0.5         0.7
#> 3276         3  89974.110    final        1200       1, 1   0.3 0.5         0.7
#> 3277         3  90045.421    final        1200       1, 1   0.3 0.5         0.7
#> 3278         3  90051.617    final        1200       1, 1   0.3 0.5         0.7
#> 3279         3  90258.477    final        1200       1, 1   0.3 0.5         0.7
#> 3280         3  90270.095    final        1200       1, 1   0.3 0.5         0.7
#> 3281         3  90432.464    final        1200       1, 1   0.3 0.5         0.7
#> 3282         3  90435.167    final        1200       1, 1   0.3 0.5         0.7
#> 3283         3  90481.485    final        1200       1, 1   0.3 0.5         0.7
#> 3284         3  90926.376    final        1200       1, 1   0.3 0.5         0.7
#> 3285         3  90980.638    final        1200       1, 1   0.3 0.5         0.7
#> 3286         3  91365.815    final        1200       1, 1   0.3 0.5         0.7
#> 3287         3  91416.290    final        1200       1, 1   0.3 0.5         0.7
#> 3288         3  91447.116    final        1200       1, 1   0.3 0.5         0.7
#> 3289         3  91516.028    final        1200       1, 1   0.3 0.5         0.7
#> 3290         3  91670.005    final        1200       1, 1   0.3 0.5         0.7
#> 3291         3  91913.331    final        1200       1, 1   0.3 0.5         0.7
#> 3292         3  91994.482    final        1200       1, 1   0.3 0.5         0.7
#> 3293         3  92047.497    final        1200       1, 1   0.3 0.5         0.7
#> 3294         3  92165.726    final        1200       1, 1   0.3 0.5         0.7
#> 3295         3  92307.812    final        1200       1, 1   0.3 0.5         0.7
#> 3296         3  92644.106    final        1200       1, 1   0.3 0.5         0.7
#> 3297         3  92702.908    final        1200       1, 1   0.3 0.5         0.7
#> 3298         3  92749.540    final        1200       1, 1   0.3 0.5         0.7
#> 3299         3  92811.706    final        1200       1, 1   0.3 0.5         0.7
#> 3300         3  92930.349    final        1200       1, 1   0.3 0.5         0.7
#> 3301         3  92975.327    final        1200       1, 1   0.3 0.5         0.7
#> 3302         3  92980.651    final        1200       1, 1   0.3 0.5         0.7
#> 3303         3  93027.568    final        1200       1, 1   0.3 0.5         0.7
#> 3304         3  93046.114    final        1200       1, 1   0.3 0.5         0.7
#> 3305         3  93087.645    final        1200       1, 1   0.3 0.5         0.7
#> 3306         3  93108.663    final        1200       1, 1   0.3 0.5         0.7
#> 3307         3  93109.014    final        1200       1, 1   0.3 0.5         0.7
#> 3308         3  93140.933    final        1200       1, 1   0.3 0.5         0.7
#> 3309         3  93317.929    final        1200       1, 1   0.3 0.5         0.7
#> 3310         3  93329.362    final        1200       1, 1   0.3 0.5         0.7
#> 3311         3  93334.225    final        1200       1, 1   0.3 0.5         0.7
#> 3312         3  93346.777    final        1200       1, 1   0.3 0.5         0.7
#> 3313         3  93429.656    final        1200       1, 1   0.3 0.5         0.7
#> 3314         3  93483.673    final        1200       1, 1   0.3 0.5         0.7
#> 3315         3  93685.391    final        1200       1, 1   0.3 0.5         0.7
#> 3316         3  93693.428    final        1200       1, 1   0.3 0.5         0.7
#> 3317         3  94088.569    final        1200       1, 1   0.3 0.5         0.7
#> 3318         3  94166.056    final        1200       1, 1   0.3 0.5         0.7
#> 3319         3  94240.473    final        1200       1, 1   0.3 0.5         0.7
#> 3320         3  94537.169    final        1200       1, 1   0.3 0.5         0.7
#> 3321         3  94541.551    final        1200       1, 1   0.3 0.5         0.7
#> 3322         3  94831.753    final        1200       1, 1   0.3 0.5         0.7
#> 3323         3  94920.197    final        1200       1, 1   0.3 0.5         0.7
#> 3324         3  95017.546    final        1200       1, 1   0.3 0.5         0.7
#> 3325         3  95199.569    final        1200       1, 1   0.3 0.5         0.7
#> 3326         3  95242.689    final        1200       1, 1   0.3 0.5         0.7
#> 3327         3  95313.607    final        1200       1, 1   0.3 0.5         0.7
#> 3328         3  95325.738    final        1200       1, 1   0.3 0.5         0.7
#> 3329         3  95358.355    final        1200       1, 1   0.3 0.5         0.7
#> 3330         3  95450.784    final        1200       1, 1   0.3 0.5         0.7
#> 3331         3  95484.881    final        1200       1, 1   0.3 0.5         0.7
#> 3332         3  95588.286    final        1200       1, 1   0.3 0.5         0.7
#> 3333         3  95617.319    final        1200       1, 1   0.3 0.5         0.7
#> 3334         3  95626.127    final        1200       1, 1   0.3 0.5         0.7
#> 3335         3  95633.095    final        1200       1, 1   0.3 0.5         0.7
#> 3336         3  95760.705    final        1200       1, 1   0.3 0.5         0.7
#> 3337         3  95862.569    final        1200       1, 1   0.3 0.5         0.7
#> 3338         3  95956.098    final        1200       1, 1   0.3 0.5         0.7
#> 3339         3  95967.239    final        1200       1, 1   0.3 0.5         0.7
#> 3340         3  96086.918    final        1200       1, 1   0.3 0.5         0.7
#> 3341         3  96435.242    final        1200       1, 1   0.3 0.5         0.7
#> 3342         3  96495.174    final        1200       1, 1   0.3 0.5         0.7
#> 3343         3  96504.391    final        1200       1, 1   0.3 0.5         0.7
#> 3344         3  97031.209    final        1200       1, 1   0.3 0.5         0.7
#> 3345         3  97076.123    final        1200       1, 1   0.3 0.5         0.7
#> 3346         3  97076.891    final        1200       1, 1   0.3 0.5         0.7
#> 3347         3  97199.070    final        1200       1, 1   0.3 0.5         0.7
#> 3348         3  97831.035    final        1200       1, 1   0.3 0.5         0.7
#> 3349         3  97955.468    final        1200       1, 1   0.3 0.5         0.7
#> 3350         3  98019.587    final        1200       1, 1   0.3 0.5         0.7
#> 3351         3  98036.574    final        1200       1, 1   0.3 0.5         0.7
#> 3352         3  98092.597    final        1200       1, 1   0.3 0.5         0.7
#> 3353         3  98140.431    final        1200       1, 1   0.3 0.5         0.7
#> 3354         3  98296.022    final        1200       1, 1   0.3 0.5         0.7
#> 3355         3  98298.260    final        1200       1, 1   0.3 0.5         0.7
#> 3356         3  98537.763    final        1200       1, 1   0.3 0.5         0.7
#> 3357         3  98600.601    final        1200       1, 1   0.3 0.5         0.7
#> 3358         3  98604.522    final        1200       1, 1   0.3 0.5         0.7
#> 3359         3  98702.386    final        1200       1, 1   0.3 0.5         0.7
#> 3360         3  98766.624    final        1200       1, 1   0.3 0.5         0.7
#> 3361         3  98881.075    final        1200       1, 1   0.3 0.5         0.7
#> 3362         3  98912.777    final        1200       1, 1   0.3 0.5         0.7
#> 3363         3  98994.040    final        1200       1, 1   0.3 0.5         0.7
#> 3364         3  99238.568    final        1200       1, 1   0.3 0.5         0.7
#> 3365         3  99455.444    final        1200       1, 1   0.3 0.5         0.7
#> 3366         3  99510.545    final        1200       1, 1   0.3 0.5         0.7
#> 3367         3  99538.672    final        1200       1, 1   0.3 0.5         0.7
#> 3368         3  99582.727    final        1200       1, 1   0.3 0.5         0.7
#> 3369         3  99629.092    final        1200       1, 1   0.3 0.5         0.7
#> 3370         3  99764.986    final        1200       1, 1   0.3 0.5         0.7
#> 3371         3  99789.474    final        1200       1, 1   0.3 0.5         0.7
#> 3372         3  99844.162    final        1200       1, 1   0.3 0.5         0.7
#> 3373         3  99886.054    final        1200       1, 1   0.3 0.5         0.7
#> 3374         3  99960.499    final        1200       1, 1   0.3 0.5         0.7
#> 3375         3 100000.057    final        1200       1, 1   0.3 0.5         0.7
#> 3376         3 100164.623    final        1200       1, 1   0.3 0.5         0.7
#> 3377         3 100273.596    final        1200       1, 1   0.3 0.5         0.7
#> 3378         3 100331.749    final        1200       1, 1   0.3 0.5         0.7
#> 3379         3 100343.950    final        1200       1, 1   0.3 0.5         0.7
#> 3380         3 100516.140    final        1200       1, 1   0.3 0.5         0.7
#> 3381         3 100677.997    final        1200       1, 1   0.3 0.5         0.7
#> 3382         3 100755.717    final        1200       1, 1   0.3 0.5         0.7
#> 3383         3 100903.728    final        1200       1, 1   0.3 0.5         0.7
#> 3384         3 101192.082    final        1200       1, 1   0.3 0.5         0.7
#> 3385         3 101231.272    final        1200       1, 1   0.3 0.5         0.7
#> 3386         3 101329.037    final        1200       1, 1   0.3 0.5         0.7
#> 3387         3 101335.136    final        1200       1, 1   0.3 0.5         0.7
#> 3388         3 101550.219    final        1200       1, 1   0.3 0.5         0.7
#> 3389         3 101602.960    final        1200       1, 1   0.3 0.5         0.7
#> 3390         3 101681.762    final        1200       1, 1   0.3 0.5         0.7
#> 3391         3 101773.375    final        1200       1, 1   0.3 0.5         0.7
#> 3392         3 101860.649    final        1200       1, 1   0.3 0.5         0.7
#> 3393         3 102001.473    final        1200       1, 1   0.3 0.5         0.7
#> 3394         3 102067.177    final        1200       1, 1   0.3 0.5         0.7
#> 3395         3 102097.763    final        1200       1, 1   0.3 0.5         0.7
#> 3396         3 102232.497    final        1200       1, 1   0.3 0.5         0.7
#> 3397         3 102265.628    final        1200       1, 1   0.3 0.5         0.7
#> 3398         3 102544.536    final        1200       1, 1   0.3 0.5         0.7
#> 3399         3 102547.714    final        1200       1, 1   0.3 0.5         0.7
#> 3400         3 102643.931    final        1200       1, 1   0.3 0.5         0.7
#> 3401         3 102713.397    final        1200       1, 1   0.3 0.5         0.7
#> 3402         3 102744.051    final        1200       1, 1   0.3 0.5         0.7
#> 3403         3 102744.391    final        1200       1, 1   0.3 0.5         0.7
#> 3404         3 102762.101    final        1200       1, 1   0.3 0.5         0.7
#> 3405         3 102985.908    final        1200       1, 1   0.3 0.5         0.7
#> 3406         3 103014.467    final        1200       1, 1   0.3 0.5         0.7
#> 3407         3 103040.420    final        1200       1, 1   0.3 0.5         0.7
#> 3408         3 103040.500    final        1200       1, 1   0.3 0.5         0.7
#> 3409         3 103098.504    final        1200       1, 1   0.3 0.5         0.7
#> 3410         3 103190.138    final        1200       1, 1   0.3 0.5         0.7
#> 3411         3 103195.721    final        1200       1, 1   0.3 0.5         0.7
#> 3412         3 103322.526    final        1200       1, 1   0.3 0.5         0.7
#> 3413         3 103348.309    final        1200       1, 1   0.3 0.5         0.7
#> 3414         3 103417.831    final        1200       1, 1   0.3 0.5         0.7
#> 3415         3 103418.463    final        1200       1, 1   0.3 0.5         0.7
#> 3416         3 103547.002    final        1200       1, 1   0.3 0.5         0.7
#> 3417         3 103785.107    final        1200       1, 1   0.3 0.5         0.7
#> 3418         3 103810.937    final        1200       1, 1   0.3 0.5         0.7
#> 3419         3 103909.101    final        1200       1, 1   0.3 0.5         0.7
#> 3420         3 104066.122    final        1200       1, 1   0.3 0.5         0.7
#> 3421         3 104068.542    final        1200       1, 1   0.3 0.5         0.7
#> 3422         3 104272.422    final        1200       1, 1   0.3 0.5         0.7
#> 3423         3 104385.150    final        1200       1, 1   0.3 0.5         0.7
#> 3424         3 104494.041    final        1200       1, 1   0.3 0.5         0.7
#> 3425         3 104628.213    final        1200       1, 1   0.3 0.5         0.7
#> 3426         3 104776.147    final        1200       1, 1   0.3 0.5         0.7
#> 3427         3 104808.235    final        1200       1, 1   0.3 0.5         0.7
#> 3428         3 104953.516    final        1200       1, 1   0.3 0.5         0.7
#> 3429         3 105064.343    final        1200       1, 1   0.3 0.5         0.7
#> 3430         3 105174.370    final        1200       1, 1   0.3 0.5         0.7
#> 3431         3 105210.719    final        1200       1, 1   0.3 0.5         0.7
#> 3432         3 105249.940    final        1200       1, 1   0.3 0.5         0.7
#> 3433         3 105260.874    final        1200       1, 1   0.3 0.5         0.7
#> 3434         3 105274.979    final        1200       1, 1   0.3 0.5         0.7
#> 3435         3 105328.081    final        1200       1, 1   0.3 0.5         0.7
#> 3436         3 105416.676    final        1200       1, 1   0.3 0.5         0.7
#> 3437         3 105541.729    final        1200       1, 1   0.3 0.5         0.7
#> 3438         3 105722.448    final        1200       1, 1   0.3 0.5         0.7
#> 3439         3 105798.999    final        1200       1, 1   0.3 0.5         0.7
#> 3440         3 105880.914    final        1200       1, 1   0.3 0.5         0.7
#> 3441         3 105904.262    final        1200       1, 1   0.3 0.5         0.7
#> 3442         3 105926.526    final        1200       1, 1   0.3 0.5         0.7
#> 3443         3 106042.134    final        1200       1, 1   0.3 0.5         0.7
#> 3444         3 106089.812    final        1200       1, 1   0.3 0.5         0.7
#> 3445         3 106351.247    final        1200       1, 1   0.3 0.5         0.7
#> 3446         3 106382.045    final        1200       1, 1   0.3 0.5         0.7
#> 3447         3 106388.001    final        1200       1, 1   0.3 0.5         0.7
#> 3448         3 106419.191    final        1200       1, 1   0.3 0.5         0.7
#> 3449         3 106497.775    final        1200       1, 1   0.3 0.5         0.7
#> 3450         3 106618.887    final        1200       1, 1   0.3 0.5         0.7
#> 3451         3 106642.118    final        1200       1, 1   0.3 0.5         0.7
#> 3452         3 106714.984    final        1200       1, 1   0.3 0.5         0.7
#> 3453         3 106833.252    final        1200       1, 1   0.3 0.5         0.7
#> 3454         3 107064.573    final        1200       1, 1   0.3 0.5         0.7
#> 3455         3 107094.673    final        1200       1, 1   0.3 0.5         0.7
#> 3456         3 107124.839    final        1200       1, 1   0.3 0.5         0.7
#> 3457         3 107147.852    final        1200       1, 1   0.3 0.5         0.7
#> 3458         3 107165.430    final        1200       1, 1   0.3 0.5         0.7
#> 3459         3 107166.014    final        1200       1, 1   0.3 0.5         0.7
#> 3460         3 107238.870    final        1200       1, 1   0.3 0.5         0.7
#> 3461         3 107249.319    final        1200       1, 1   0.3 0.5         0.7
#> 3462         3 107334.403    final        1200       1, 1   0.3 0.5         0.7
#> 3463         3 107582.673    final        1200       1, 1   0.3 0.5         0.7
#> 3464         3 107592.457    final        1200       1, 1   0.3 0.5         0.7
#> 3465         3 107599.915    final        1200       1, 1   0.3 0.5         0.7
#> 3466         3 107601.441    final        1200       1, 1   0.3 0.5         0.7
#> 3467         3 107685.484    final        1200       1, 1   0.3 0.5         0.7
#> 3468         3 107795.084    final        1200       1, 1   0.3 0.5         0.7
#> 3469         3 107832.455    final        1200       1, 1   0.3 0.5         0.7
#> 3470         3 107833.343    final        1200       1, 1   0.3 0.5         0.7
#> 3471         3 107956.312    final        1200       1, 1   0.3 0.5         0.7
#> 3472         3 108079.495    final        1200       1, 1   0.3 0.5         0.7
#> 3473         3 108091.326    final        1200       1, 1   0.3 0.5         0.7
#> 3474         3 108330.764    final        1200       1, 1   0.3 0.5         0.7
#> 3475         3 108552.800    final        1200       1, 1   0.3 0.5         0.7
#> 3476         3 108554.600    final        1200       1, 1   0.3 0.5         0.7
#> 3477         3 108640.494    final        1200       1, 1   0.3 0.5         0.7
#> 3478         3 108656.332    final        1200       1, 1   0.3 0.5         0.7
#> 3479         3 108732.919    final        1200       1, 1   0.3 0.5         0.7
#> 3480         3 108807.408    final        1200       1, 1   0.3 0.5         0.7
#> 3481         3 108859.498    final        1200       1, 1   0.3 0.5         0.7
#> 3482         3 109069.578    final        1200       1, 1   0.3 0.5         0.7
#> 3483         3 109111.636    final        1200       1, 1   0.3 0.5         0.7
#> 3484         3 109222.507    final        1200       1, 1   0.3 0.5         0.7
#> 3485         3 109297.690    final        1200       1, 1   0.3 0.5         0.7
#> 3486         3 109309.089    final        1200       1, 1   0.3 0.5         0.7
#> 3487         3 109457.691    final        1200       1, 1   0.3 0.5         0.7
#> 3488         3 109529.741    final        1200       1, 1   0.3 0.5         0.7
#> 3489         3 109564.581    final        1200       1, 1   0.3 0.5         0.7
#> 3490         3 109575.815    final        1200       1, 1   0.3 0.5         0.7
#> 3491         3 109731.805    final        1200       1, 1   0.3 0.5         0.7
#> 3492         3 109743.829    final        1200       1, 1   0.3 0.5         0.7
#> 3493         3 109838.744    final        1200       1, 1   0.3 0.5         0.7
#> 3494         3 109887.824    final        1200       1, 1   0.3 0.5         0.7
#> 3495         3 109979.626    final        1200       1, 1   0.3 0.5         0.7
#> 3496         3 110138.265    final        1200       1, 1   0.3 0.5         0.7
#> 3497         3 110235.425    final        1200       1, 1   0.3 0.5         0.7
#> 3498         3 110362.854    final        1200       1, 1   0.3 0.5         0.7
#> 3499         3 110371.428    final        1200       1, 1   0.3 0.5         0.7
#> 3500         3 110515.219    final        1200       1, 1   0.3 0.5         0.7
#> 3501         3 110524.429    final        1200       1, 1   0.3 0.5         0.7
#> 3502         3 111036.251    final        1200       1, 1   0.3 0.5         0.7
#> 3503         3 111094.592    final        1200       1, 1   0.3 0.5         0.7
#> 3504         3 111099.335    final        1200       1, 1   0.3 0.5         0.7
#> 3505         3 111110.642    final        1200       1, 1   0.3 0.5         0.7
#> 3506         3 111132.730    final        1200       1, 1   0.3 0.5         0.7
#> 3507         3 111134.482    final        1200       1, 1   0.3 0.5         0.7
#> 3508         3 111138.161    final        1200       1, 1   0.3 0.5         0.7
#> 3509         3 111440.573    final        1200       1, 1   0.3 0.5         0.7
#> 3510         3 111573.807    final        1200       1, 1   0.3 0.5         0.7
#> 3511         3 111602.993    final        1200       1, 1   0.3 0.5         0.7
#> 3512         3 111689.766    final        1200       1, 1   0.3 0.5         0.7
#> 3513         3 111790.384    final        1200       1, 1   0.3 0.5         0.7
#> 3514         3 112062.148    final        1200       1, 1   0.3 0.5         0.7
#> 3515         3 112248.162    final        1200       1, 1   0.3 0.5         0.7
#> 3516         3 112345.646    final        1200       1, 1   0.3 0.5         0.7
#> 3517         3 112350.091    final        1200       1, 1   0.3 0.5         0.7
#> 3518         3 112523.429    final        1200       1, 1   0.3 0.5         0.7
#> 3519         3 112680.218    final        1200       1, 1   0.3 0.5         0.7
#> 3520         3 112906.208    final        1200       1, 1   0.3 0.5         0.7
#> 3521         3 113271.963    final        1200       1, 1   0.3 0.5         0.7
#> 3522         3 113359.804    final        1200       1, 1   0.3 0.5         0.7
#> 3523         3 113428.469    final        1200       1, 1   0.3 0.5         0.7
#> 3524         3 113460.211    final        1200       1, 1   0.3 0.5         0.7
#> 3525         3 113501.796    final        1200       1, 1   0.3 0.5         0.7
#> 3526         3 113529.856    final        1200       1, 1   0.3 0.5         0.7
#> 3527         3 113631.131    final        1200       1, 1   0.3 0.5         0.7
#> 3528         3 113768.714    final        1200       1, 1   0.3 0.5         0.7
#> 3529         3 113769.526    final        1200       1, 1   0.3 0.5         0.7
#> 3530         3 113879.340    final        1200       1, 1   0.3 0.5         0.7
#> 3531         3 113885.976    final        1200       1, 1   0.3 0.5         0.7
#> 3532         3 113960.145    final        1200       1, 1   0.3 0.5         0.7
#> 3533         3 114056.451    final        1200       1, 1   0.3 0.5         0.7
#> 3534         3 114113.991    final        1200       1, 1   0.3 0.5         0.7
#> 3535         3 114255.967    final        1200       1, 1   0.3 0.5         0.7
#> 3536         3 114377.080    final        1200       1, 1   0.3 0.5         0.7
#> 3537         3 114462.210    final        1200       1, 1   0.3 0.5         0.7
#> 3538         3 114506.773    final        1200       1, 1   0.3 0.5         0.7
#> 3539         3 114721.009    final        1200       1, 1   0.3 0.5         0.7
#> 3540         3 115094.741    final        1200       1, 1   0.3 0.5         0.7
#> 3541         3 115111.425    final        1200       1, 1   0.3 0.5         0.7
#> 3542         3 115174.317    final        1200       1, 1   0.3 0.5         0.7
#> 3543         3 115188.353    final        1200       1, 1   0.3 0.5         0.7
#> 3544         3 115189.637    final        1200       1, 1   0.3 0.5         0.7
#> 3545         3 115193.713    final        1200       1, 1   0.3 0.5         0.7
#> 3546         3 115197.571    final        1200       1, 1   0.3 0.5         0.7
#> 3547         3 115200.820    final        1200       1, 1   0.3 0.5         0.7
#> 3548         3 115407.695    final        1200       1, 1   0.3 0.5         0.7
#> 3549         3 115652.848    final        1200       1, 1   0.3 0.5         0.7
#> 3550         3 115660.257    final        1200       1, 1   0.3 0.5         0.7
#> 3551         3 115753.275    final        1200       1, 1   0.3 0.5         0.7
#> 3552         3 115826.785    final        1200       1, 1   0.3 0.5         0.7
#> 3553         3 115870.556    final        1200       1, 1   0.3 0.5         0.7
#> 3554         3 116236.154    final        1200       1, 1   0.3 0.5         0.7
#> 3555         3 116240.217    final        1200       1, 1   0.3 0.5         0.7
#> 3556         3 116369.369    final        1200       1, 1   0.3 0.5         0.7
#> 3557         3 116638.657    final        1200       1, 1   0.3 0.5         0.7
#> 3558         3 116749.748    final        1200       1, 1   0.3 0.5         0.7
#> 3559         3 116806.896    final        1200       1, 1   0.3 0.5         0.7
#> 3560         3 116844.722    final        1200       1, 1   0.3 0.5         0.7
#> 3561         3 116845.168    final        1200       1, 1   0.3 0.5         0.7
#> 3562         3 116920.512    final        1200       1, 1   0.3 0.5         0.7
#> 3563         3 117063.213    final        1200       1, 1   0.3 0.5         0.7
#> 3564         3 117158.135    final        1200       1, 1   0.3 0.5         0.7
#> 3565         3 117198.723    final        1200       1, 1   0.3 0.5         0.7
#> 3566         3 117246.049    final        1200       1, 1   0.3 0.5         0.7
#> 3567         3 117312.898    final        1200       1, 1   0.3 0.5         0.7
#> 3568         3 117429.398    final        1200       1, 1   0.3 0.5         0.7
#> 3569         3 117436.095    final        1200       1, 1   0.3 0.5         0.7
#> 3570         3 117451.199    final        1200       1, 1   0.3 0.5         0.7
#> 3571         3 117488.425    final        1200       1, 1   0.3 0.5         0.7
#> 3572         3 117524.532    final        1200       1, 1   0.3 0.5         0.7
#> 3573         3 117775.389    final        1200       1, 1   0.3 0.5         0.7
#> 3574         3 117841.157    final        1200       1, 1   0.3 0.5         0.7
#> 3575         3 117842.271    final        1200       1, 1   0.3 0.5         0.7
#> 3576         3 118126.482    final        1200       1, 1   0.3 0.5         0.7
#>      n_total n_pbo n_trt    p_ttest_y    logrank_p    hr_cox  hr_ci_lo
#> 1       1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 2       1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 3       1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 4       1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 5       1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 6       1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 7       1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 8       1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 9       1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 10      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 11      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 12      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 13      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 14      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 15      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 16      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 17      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 18      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 19      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 20      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 21      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 22      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 23      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 24      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 25      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 26      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 27      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 28      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 29      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 30      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 31      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 32      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 33      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 34      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 35      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 36      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 37      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 38      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 39      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 40      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 41      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 42      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 43      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 44      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 45      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 46      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 47      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 48      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 49      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 50      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 51      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 52      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 53      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 54      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 55      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 56      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 57      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 58      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 59      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 60      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 61      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 62      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 63      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 64      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 65      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 66      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 67      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 68      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 69      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 70      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 71      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 72      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 73      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 74      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 75      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 76      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 77      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 78      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 79      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 80      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 81      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 82      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 83      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 84      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 85      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 86      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 87      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 88      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 89      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 90      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 91      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 92      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 93      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 94      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 95      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 96      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 97      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 98      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 99      1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 100     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 101     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 102     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 103     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 104     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 105     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 106     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 107     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 108     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 109     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 110     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 111     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 112     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 113     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 114     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 115     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 116     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 117     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 118     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 119     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 120     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 121     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 122     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 123     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 124     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 125     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 126     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 127     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 128     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 129     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 130     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 131     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 132     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 133     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 134     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 135     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 136     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 137     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 138     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 139     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 140     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 141     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 142     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 143     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 144     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 145     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 146     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 147     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 148     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 149     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 150     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 151     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 152     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 153     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 154     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 155     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 156     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 157     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 158     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 159     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 160     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 161     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 162     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 163     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 164     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 165     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 166     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 167     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 168     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 169     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 170     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 171     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 172     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 173     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 174     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 175     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 176     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 177     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 178     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 179     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 180     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 181     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 182     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 183     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 184     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 185     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 186     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 187     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 188     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 189     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 190     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 191     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 192     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 193     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 194     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 195     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 196     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 197     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 198     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 199     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 200     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 201     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 202     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 203     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 204     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 205     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 206     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 207     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 208     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 209     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 210     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 211     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 212     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 213     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 214     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 215     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 216     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 217     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 218     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 219     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 220     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 221     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 222     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 223     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 224     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 225     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 226     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 227     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 228     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 229     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 230     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 231     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 232     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 233     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 234     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 235     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 236     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 237     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 238     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 239     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 240     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 241     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 242     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 243     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 244     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 245     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 246     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 247     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 248     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 249     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 250     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 251     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 252     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 253     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 254     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 255     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 256     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 257     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 258     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 259     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 260     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 261     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 262     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 263     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 264     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 265     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 266     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 267     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 268     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 269     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 270     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 271     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 272     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 273     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 274     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 275     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 276     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 277     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 278     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 279     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 280     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 281     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 282     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 283     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 284     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 285     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 286     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 287     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 288     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 289     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 290     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 291     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 292     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 293     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 294     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 295     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 296     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 297     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 298     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 299     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 300     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 301     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 302     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 303     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 304     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 305     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 306     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 307     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 308     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 309     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 310     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 311     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 312     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 313     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 314     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 315     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 316     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 317     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 318     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 319     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 320     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 321     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 322     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 323     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 324     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 325     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 326     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 327     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 328     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 329     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 330     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 331     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 332     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 333     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 334     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 335     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 336     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 337     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 338     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 339     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 340     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 341     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 342     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 343     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 344     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 345     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 346     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 347     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 348     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 349     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 350     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 351     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 352     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 353     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 354     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 355     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 356     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 357     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 358     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 359     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 360     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 361     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 362     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 363     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 364     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 365     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 366     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 367     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 368     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 369     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 370     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 371     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 372     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 373     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 374     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 375     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 376     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 377     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 378     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 379     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 380     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 381     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 382     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 383     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 384     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 385     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 386     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 387     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 388     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 389     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 390     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 391     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 392     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 393     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 394     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 395     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 396     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 397     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 398     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 399     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 400     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 401     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 402     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 403     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 404     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 405     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 406     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 407     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 408     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 409     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 410     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 411     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 412     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 413     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 414     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 415     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 416     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 417     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 418     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 419     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 420     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 421     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 422     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 423     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 424     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 425     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 426     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 427     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 428     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 429     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 430     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 431     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 432     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 433     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 434     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 435     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 436     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 437     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 438     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 439     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 440     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 441     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 442     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 443     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 444     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 445     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 446     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 447     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 448     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 449     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 450     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 451     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 452     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 453     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 454     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 455     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 456     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 457     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 458     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 459     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 460     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 461     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 462     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 463     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 464     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 465     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 466     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 467     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 468     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 469     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 470     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 471     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 472     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 473     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 474     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 475     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 476     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 477     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 478     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 479     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 480     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 481     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 482     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 483     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 484     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 485     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 486     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 487     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 488     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 489     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 490     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 491     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 492     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 493     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 494     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 495     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 496     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 497     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 498     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 499     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 500     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 501     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 502     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 503     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 504     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 505     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 506     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 507     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 508     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 509     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 510     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 511     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 512     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 513     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 514     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 515     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 516     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 517     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 518     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 519     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 520     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 521     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 522     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 523     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 524     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 525     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 526     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 527     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 528     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 529     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 530     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 531     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 532     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 533     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 534     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 535     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 536     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 537     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 538     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 539     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 540     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 541     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 542     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 543     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 544     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 545     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 546     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 547     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 548     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 549     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 550     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 551     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 552     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 553     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 554     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 555     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 556     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 557     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 558     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 559     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 560     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 561     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 562     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 563     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 564     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 565     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 566     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 567     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 568     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 569     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 570     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 571     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 572     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 573     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 574     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 575     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 576     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 577     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 578     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 579     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 580     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 581     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 582     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 583     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 584     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 585     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 586     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 587     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 588     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 589     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 590     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 591     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 592     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 593     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 594     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 595     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 596     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 597     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 598     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 599     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 600     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 601     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 602     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 603     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 604     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 605     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 606     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 607     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 608     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 609     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 610     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 611     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 612     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 613     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 614     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 615     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 616     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 617     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 618     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 619     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 620     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 621     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 622     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 623     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 624     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 625     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 626     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 627     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 628     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 629     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 630     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 631     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 632     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 633     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 634     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 635     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 636     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 637     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 638     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 639     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 640     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 641     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 642     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 643     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 644     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 645     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 646     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 647     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 648     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 649     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 650     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 651     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 652     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 653     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 654     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 655     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 656     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 657     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 658     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 659     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 660     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 661     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 662     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 663     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 664     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 665     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 666     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 667     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 668     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 669     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 670     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 671     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 672     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 673     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 674     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 675     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 676     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 677     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 678     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 679     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 680     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 681     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 682     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 683     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 684     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 685     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 686     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 687     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 688     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 689     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 690     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 691     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 692     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 693     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 694     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 695     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 696     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 697     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 698     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 699     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 700     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 701     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 702     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 703     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 704     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 705     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 706     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 707     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 708     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 709     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 710     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 711     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 712     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 713     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 714     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 715     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 716     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 717     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 718     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 719     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 720     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 721     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 722     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 723     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 724     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 725     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 726     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 727     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 728     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 729     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 730     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 731     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 732     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 733     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 734     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 735     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 736     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 737     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 738     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 739     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 740     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 741     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 742     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 743     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 744     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 745     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 746     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 747     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 748     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 749     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 750     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 751     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 752     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 753     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 754     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 755     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 756     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 757     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 758     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 759     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 760     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 761     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 762     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 763     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 764     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 765     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 766     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 767     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 768     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 769     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 770     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 771     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 772     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 773     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 774     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 775     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 776     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 777     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 778     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 779     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 780     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 781     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 782     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 783     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 784     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 785     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 786     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 787     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 788     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 789     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 790     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 791     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 792     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 793     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 794     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 795     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 796     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 797     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 798     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 799     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 800     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 801     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 802     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 803     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 804     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 805     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 806     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 807     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 808     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 809     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 810     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 811     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 812     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 813     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 814     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 815     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 816     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 817     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 818     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 819     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 820     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 821     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 822     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 823     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 824     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 825     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 826     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 827     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 828     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 829     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 830     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 831     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 832     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 833     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 834     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 835     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 836     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 837     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 838     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 839     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 840     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 841     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 842     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 843     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 844     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 845     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 846     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 847     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 848     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 849     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 850     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 851     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 852     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 853     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 854     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 855     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 856     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 857     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 858     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 859     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 860     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 861     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 862     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 863     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 864     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 865     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 866     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 867     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 868     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 869     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 870     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 871     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 872     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 873     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 874     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 875     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 876     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 877     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 878     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 879     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 880     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 881     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 882     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 883     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 884     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 885     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 886     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 887     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 888     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 889     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 890     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 891     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 892     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 893     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 894     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 895     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 896     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 897     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 898     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 899     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 900     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 901     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 902     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 903     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 904     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 905     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 906     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 907     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 908     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 909     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 910     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 911     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 912     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 913     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 914     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 915     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 916     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 917     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 918     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 919     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 920     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 921     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 922     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 923     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 924     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 925     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 926     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 927     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 928     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 929     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 930     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 931     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 932     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 933     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 934     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 935     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 936     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 937     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 938     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 939     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 940     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 941     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 942     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 943     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 944     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 945     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 946     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 947     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 948     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 949     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 950     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 951     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 952     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 953     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 954     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 955     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 956     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 957     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 958     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 959     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 960     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 961     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 962     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 963     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 964     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 965     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 966     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 967     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 968     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 969     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 970     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 971     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 972     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 973     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 974     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 975     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 976     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 977     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 978     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 979     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 980     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 981     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 982     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 983     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 984     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 985     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 986     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 987     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 988     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 989     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 990     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 991     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 992     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 993     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 994     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 995     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 996     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 997     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 998     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 999     1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1000    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1001    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1002    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1003    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1004    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1005    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1006    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1007    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1008    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1009    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1010    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1011    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1012    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1013    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1014    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1015    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1016    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1017    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1018    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1019    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1020    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1021    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1022    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1023    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1024    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1025    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1026    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1027    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1028    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1029    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1030    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1031    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1032    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1033    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1034    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1035    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1036    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1037    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1038    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1039    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1040    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1041    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1042    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1043    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1044    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1045    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1046    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1047    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1048    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1049    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1050    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1051    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1052    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1053    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1054    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1055    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1056    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1057    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1058    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1059    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1060    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1061    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1062    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1063    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1064    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1065    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1066    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1067    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1068    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1069    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1070    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1071    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1072    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1073    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1074    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1075    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1076    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1077    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1078    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1079    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1080    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1081    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1082    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1083    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1084    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1085    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1086    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1087    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1088    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1089    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1090    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1091    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1092    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1093    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1094    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1095    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1096    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1097    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1098    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1099    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1100    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1101    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1102    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1103    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1104    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1105    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1106    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1107    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1108    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1109    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1110    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1111    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1112    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1113    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1114    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1115    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1116    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1117    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1118    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1119    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1120    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1121    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1122    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1123    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1124    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1125    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1126    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1127    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1128    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1129    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1130    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1131    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1132    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1133    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1134    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1135    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1136    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1137    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1138    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1139    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1140    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1141    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1142    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1143    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1144    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1145    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1146    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1147    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1148    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1149    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1150    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1151    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1152    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1153    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1154    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1155    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1156    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1157    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1158    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1159    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1160    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1161    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1162    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1163    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1164    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1165    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1166    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1167    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1168    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1169    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1170    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1171    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1172    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1173    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1174    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1175    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1176    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1177    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1178    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1179    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1180    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1181    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1182    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1183    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1184    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1185    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1186    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1187    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1188    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1189    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1190    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1191    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1192    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1193    1200   600   600 1.437333e-07 1.352989e-05 0.7768063 0.6930741
#> 1194    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1195    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1196    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1197    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1198    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1199    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1200    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1201    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1202    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1203    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1204    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1205    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1206    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1207    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1208    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1209    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1210    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1211    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1212    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1213    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1214    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1215    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1216    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1217    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1218    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1219    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1220    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1221    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1222    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1223    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1224    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1225    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1226    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1227    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1228    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1229    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1230    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1231    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1232    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1233    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1234    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1235    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1236    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1237    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1238    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1239    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1240    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1241    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1242    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1243    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1244    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1245    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1246    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1247    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1248    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1249    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1250    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1251    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1252    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1253    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1254    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1255    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1256    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1257    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1258    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1259    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1260    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1261    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1262    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1263    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1264    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1265    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1266    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1267    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1268    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1269    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1270    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1271    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1272    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1273    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1274    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1275    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1276    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1277    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1278    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1279    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1280    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1281    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1282    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1283    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1284    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1285    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1286    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1287    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1288    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1289    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1290    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1291    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1292    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1293    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1294    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1295    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1296    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1297    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1298    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1299    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1300    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1301    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1302    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1303    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1304    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1305    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1306    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1307    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1308    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1309    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1310    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1311    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1312    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1313    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1314    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1315    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1316    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1317    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1318    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1319    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1320    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1321    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1322    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1323    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1324    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1325    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1326    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1327    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1328    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1329    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1330    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1331    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1332    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1333    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1334    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1335    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1336    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1337    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1338    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1339    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1340    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1341    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1342    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1343    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1344    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1345    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1346    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1347    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1348    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1349    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1350    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1351    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1352    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1353    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1354    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1355    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1356    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1357    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1358    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1359    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1360    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1361    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1362    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1363    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1364    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1365    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1366    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1367    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1368    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1369    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1370    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1371    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1372    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1373    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1374    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1375    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1376    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1377    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1378    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1379    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1380    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1381    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1382    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1383    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1384    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1385    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1386    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1387    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1388    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1389    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1390    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1391    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1392    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1393    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1394    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1395    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1396    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1397    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1398    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1399    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1400    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1401    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1402    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1403    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1404    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1405    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1406    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1407    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1408    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1409    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1410    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1411    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1412    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1413    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1414    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1415    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1416    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1417    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1418    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1419    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1420    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1421    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1422    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1423    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1424    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1425    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1426    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1427    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1428    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1429    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1430    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1431    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1432    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1433    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1434    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1435    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1436    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1437    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1438    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1439    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1440    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1441    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1442    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1443    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1444    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1445    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1446    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1447    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1448    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1449    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1450    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1451    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1452    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1453    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1454    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1455    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1456    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1457    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1458    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1459    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1460    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1461    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1462    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1463    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1464    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1465    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1466    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1467    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1468    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1469    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1470    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1471    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1472    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1473    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1474    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1475    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1476    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1477    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1478    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1479    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1480    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1481    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1482    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1483    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1484    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1485    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1486    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1487    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1488    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1489    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1490    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1491    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1492    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1493    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1494    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1495    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1496    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1497    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1498    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1499    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1500    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1501    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1502    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1503    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1504    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1505    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1506    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1507    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1508    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1509    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1510    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1511    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1512    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1513    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1514    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1515    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1516    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1517    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1518    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1519    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1520    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1521    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1522    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1523    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1524    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1525    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1526    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1527    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1528    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1529    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1530    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1531    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1532    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1533    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1534    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1535    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1536    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1537    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1538    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1539    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1540    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1541    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1542    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1543    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1544    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1545    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1546    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1547    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1548    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1549    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1550    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1551    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1552    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1553    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1554    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1555    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1556    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1557    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1558    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1559    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1560    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1561    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1562    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1563    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1564    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1565    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1566    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1567    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1568    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1569    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1570    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1571    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1572    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1573    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1574    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1575    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1576    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1577    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1578    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1579    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1580    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1581    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1582    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1583    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1584    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1585    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1586    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1587    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1588    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1589    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1590    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1591    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1592    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1593    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1594    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1595    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1596    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1597    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1598    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1599    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1600    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1601    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1602    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1603    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1604    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1605    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1606    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1607    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1608    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1609    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1610    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1611    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1612    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1613    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1614    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1615    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1616    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1617    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1618    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1619    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1620    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1621    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1622    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1623    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1624    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1625    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1626    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1627    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1628    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1629    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1630    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1631    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1632    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1633    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1634    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1635    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1636    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1637    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1638    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1639    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1640    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1641    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1642    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1643    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1644    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1645    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1646    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1647    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1648    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1649    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1650    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1651    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1652    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1653    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1654    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1655    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1656    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1657    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1658    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1659    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1660    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1661    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1662    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1663    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1664    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1665    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1666    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1667    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1668    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1669    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1670    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1671    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1672    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1673    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1674    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1675    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1676    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1677    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1678    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1679    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1680    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1681    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1682    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1683    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1684    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1685    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1686    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1687    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1688    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1689    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1690    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1691    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1692    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1693    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1694    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1695    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1696    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1697    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1698    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1699    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1700    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1701    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1702    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1703    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1704    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1705    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1706    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1707    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1708    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1709    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1710    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1711    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1712    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1713    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1714    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1715    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1716    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1717    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1718    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1719    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1720    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1721    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1722    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1723    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1724    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1725    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1726    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1727    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1728    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1729    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1730    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1731    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1732    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1733    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1734    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1735    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1736    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1737    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1738    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1739    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1740    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1741    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1742    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1743    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1744    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1745    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1746    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1747    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1748    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1749    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1750    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1751    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1752    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1753    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1754    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1755    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1756    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1757    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1758    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1759    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1760    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1761    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1762    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1763    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1764    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1765    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1766    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1767    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1768    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1769    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1770    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1771    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1772    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1773    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1774    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1775    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1776    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1777    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1778    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1779    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1780    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1781    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1782    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1783    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1784    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1785    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1786    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1787    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1788    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1789    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1790    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1791    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1792    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1793    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1794    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1795    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1796    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1797    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1798    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1799    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1800    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1801    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1802    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1803    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1804    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1805    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1806    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1807    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1808    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1809    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1810    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1811    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1812    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1813    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1814    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1815    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1816    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1817    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1818    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1819    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1820    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1821    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1822    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1823    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1824    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1825    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1826    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1827    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1828    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1829    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1830    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1831    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1832    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1833    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1834    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1835    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1836    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1837    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1838    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1839    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1840    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1841    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1842    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1843    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1844    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1845    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1846    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1847    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1848    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1849    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1850    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1851    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1852    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1853    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1854    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1855    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1856    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1857    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1858    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1859    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1860    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1861    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1862    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1863    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1864    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1865    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1866    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1867    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1868    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1869    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1870    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1871    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1872    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1873    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1874    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1875    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1876    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1877    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1878    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1879    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1880    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1881    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1882    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1883    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1884    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1885    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1886    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1887    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1888    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1889    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1890    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1891    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1892    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1893    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1894    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1895    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1896    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1897    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1898    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1899    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1900    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1901    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1902    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1903    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1904    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1905    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1906    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1907    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1908    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1909    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1910    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1911    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1912    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1913    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1914    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1915    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1916    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1917    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1918    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1919    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1920    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1921    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1922    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1923    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1924    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1925    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1926    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1927    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1928    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1929    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1930    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1931    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1932    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1933    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1934    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1935    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1936    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1937    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1938    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1939    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1940    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1941    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1942    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1943    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1944    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1945    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1946    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1947    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1948    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1949    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1950    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1951    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1952    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1953    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1954    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1955    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1956    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1957    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1958    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1959    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1960    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1961    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1962    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1963    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1964    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1965    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1966    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1967    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1968    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1969    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1970    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1971    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1972    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1973    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1974    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1975    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1976    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1977    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1978    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1979    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1980    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1981    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1982    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1983    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1984    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1985    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1986    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1987    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1988    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1989    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1990    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1991    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1992    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1993    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1994    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1995    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1996    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1997    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1998    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 1999    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2000    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2001    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2002    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2003    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2004    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2005    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2006    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2007    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2008    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2009    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2010    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2011    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2012    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2013    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2014    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2015    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2016    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2017    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2018    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2019    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2020    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2021    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2022    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2023    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2024    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2025    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2026    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2027    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2028    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2029    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2030    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2031    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2032    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2033    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2034    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2035    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2036    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2037    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2038    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2039    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2040    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2041    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2042    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2043    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2044    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2045    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2046    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2047    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2048    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2049    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2050    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2051    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2052    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2053    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2054    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2055    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2056    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2057    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2058    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2059    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2060    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2061    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2062    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2063    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2064    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2065    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2066    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2067    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2068    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2069    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2070    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2071    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2072    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2073    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2074    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2075    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2076    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2077    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2078    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2079    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2080    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2081    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2082    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2083    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2084    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2085    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2086    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2087    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2088    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2089    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2090    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2091    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2092    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2093    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2094    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2095    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2096    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2097    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2098    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2099    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2100    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2101    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2102    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2103    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2104    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2105    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2106    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2107    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2108    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2109    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2110    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2111    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2112    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2113    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2114    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2115    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2116    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2117    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2118    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2119    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2120    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2121    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2122    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2123    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2124    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2125    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2126    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2127    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2128    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2129    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2130    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2131    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2132    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2133    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2134    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2135    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2136    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2137    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2138    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2139    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2140    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2141    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2142    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2143    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2144    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2145    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2146    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2147    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2148    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2149    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2150    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2151    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2152    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2153    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2154    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2155    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2156    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2157    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2158    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2159    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2160    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2161    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2162    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2163    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2164    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2165    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2166    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2167    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2168    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2169    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2170    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2171    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2172    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2173    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2174    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2175    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2176    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2177    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2178    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2179    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2180    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2181    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2182    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2183    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2184    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2185    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2186    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2187    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2188    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2189    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2190    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2191    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2192    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2193    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2194    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2195    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2196    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2197    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2198    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2199    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2200    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2201    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2202    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2203    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2204    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2205    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2206    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2207    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2208    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2209    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2210    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2211    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2212    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2213    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2214    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2215    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2216    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2217    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2218    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2219    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2220    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2221    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2222    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2223    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2224    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2225    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2226    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2227    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2228    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2229    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2230    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2231    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2232    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2233    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2234    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2235    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2236    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2237    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2238    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2239    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2240    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2241    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2242    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2243    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2244    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2245    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2246    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2247    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2248    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2249    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2250    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2251    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2252    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2253    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2254    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2255    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2256    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2257    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2258    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2259    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2260    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2261    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2262    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2263    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2264    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2265    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2266    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2267    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2268    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2269    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2270    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2271    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2272    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2273    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2274    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2275    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2276    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2277    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2278    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2279    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2280    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2281    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2282    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2283    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2284    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2285    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2286    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2287    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2288    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2289    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2290    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2291    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2292    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2293    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2294    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2295    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2296    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2297    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2298    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2299    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2300    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2301    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2302    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2303    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2304    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2305    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2306    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2307    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2308    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2309    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2310    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2311    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2312    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2313    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2314    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2315    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2316    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2317    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2318    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2319    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2320    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2321    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2322    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2323    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2324    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2325    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2326    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2327    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2328    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2329    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2330    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2331    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2332    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2333    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2334    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2335    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2336    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2337    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2338    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2339    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2340    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2341    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2342    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2343    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2344    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2345    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2346    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2347    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2348    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2349    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2350    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2351    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2352    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2353    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2354    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2355    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2356    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2357    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2358    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2359    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2360    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2361    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2362    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2363    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2364    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2365    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2366    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2367    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2368    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2369    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2370    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2371    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2372    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2373    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2374    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2375    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2376    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2377    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2378    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2379    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2380    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2381    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2382    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2383    1200   600   600 9.000125e-08 3.702888e-10 0.6950138 0.6198860
#> 2384    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2385    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2386    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2387    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2388    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2389    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2390    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2391    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2392    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2393    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2394    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2395    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2396    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2397    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2398    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2399    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2400    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2401    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2402    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2403    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2404    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2405    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2406    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2407    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2408    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2409    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2410    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2411    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2412    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2413    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2414    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2415    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2416    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2417    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2418    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2419    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2420    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2421    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2422    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2423    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2424    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2425    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2426    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2427    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2428    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2429    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2430    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2431    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2432    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2433    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2434    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2435    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2436    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2437    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2438    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2439    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2440    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2441    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2442    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2443    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2444    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2445    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2446    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2447    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2448    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2449    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2450    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2451    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2452    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2453    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2454    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2455    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2456    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2457    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2458    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2459    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2460    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2461    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2462    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2463    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2464    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2465    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2466    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2467    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2468    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2469    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2470    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2471    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2472    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2473    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2474    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2475    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2476    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2477    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2478    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2479    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2480    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2481    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2482    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2483    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2484    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2485    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2486    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2487    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2488    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2489    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2490    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2491    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2492    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2493    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2494    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2495    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2496    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2497    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2498    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2499    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2500    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2501    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2502    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2503    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2504    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2505    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2506    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2507    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2508    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2509    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2510    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2511    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2512    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2513    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2514    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2515    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2516    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2517    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2518    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2519    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2520    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2521    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2522    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2523    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2524    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2525    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2526    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2527    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2528    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2529    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2530    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2531    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2532    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2533    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2534    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2535    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2536    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2537    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2538    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2539    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2540    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2541    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2542    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2543    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2544    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2545    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2546    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2547    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2548    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2549    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2550    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2551    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2552    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2553    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2554    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2555    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2556    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2557    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2558    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2559    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2560    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2561    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2562    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2563    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2564    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2565    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2566    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2567    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2568    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2569    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2570    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2571    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2572    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2573    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2574    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2575    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2576    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2577    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2578    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2579    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2580    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2581    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2582    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2583    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2584    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2585    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2586    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2587    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2588    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2589    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2590    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2591    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2592    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2593    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2594    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2595    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2596    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2597    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2598    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2599    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2600    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2601    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2602    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2603    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2604    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2605    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2606    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2607    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2608    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2609    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2610    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2611    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2612    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2613    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2614    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2615    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2616    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2617    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2618    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2619    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2620    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2621    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2622    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2623    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2624    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2625    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2626    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2627    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2628    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2629    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2630    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2631    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2632    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2633    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2634    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2635    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2636    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2637    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2638    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2639    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2640    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2641    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2642    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2643    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2644    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2645    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2646    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2647    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2648    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2649    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2650    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2651    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2652    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2653    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2654    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2655    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2656    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2657    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2658    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2659    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2660    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2661    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2662    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2663    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2664    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2665    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2666    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2667    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2668    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2669    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2670    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2671    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2672    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2673    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2674    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2675    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2676    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2677    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2678    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2679    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2680    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2681    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2682    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2683    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2684    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2685    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2686    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2687    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2688    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2689    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2690    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2691    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2692    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2693    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2694    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2695    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2696    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2697    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2698    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2699    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2700    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2701    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2702    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2703    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2704    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2705    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2706    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2707    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2708    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2709    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2710    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2711    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2712    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2713    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2714    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2715    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2716    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2717    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2718    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2719    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2720    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2721    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2722    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2723    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2724    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2725    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2726    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2727    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2728    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2729    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2730    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2731    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2732    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2733    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2734    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2735    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2736    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2737    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2738    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2739    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2740    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2741    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2742    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2743    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2744    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2745    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2746    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2747    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2748    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2749    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2750    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2751    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2752    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2753    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2754    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2755    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2756    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2757    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2758    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2759    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2760    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2761    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2762    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2763    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2764    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2765    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2766    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2767    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2768    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2769    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2770    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2771    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2772    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2773    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2774    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2775    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2776    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2777    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2778    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2779    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2780    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2781    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2782    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2783    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2784    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2785    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2786    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2787    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2788    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2789    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2790    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2791    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2792    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2793    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2794    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2795    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2796    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2797    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2798    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2799    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2800    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2801    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2802    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2803    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2804    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2805    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2806    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2807    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2808    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2809    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2810    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2811    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2812    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2813    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2814    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2815    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2816    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2817    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2818    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2819    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2820    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2821    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2822    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2823    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2824    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2825    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2826    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2827    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2828    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2829    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2830    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2831    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2832    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2833    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2834    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2835    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2836    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2837    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2838    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2839    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2840    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2841    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2842    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2843    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2844    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2845    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2846    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2847    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2848    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2849    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2850    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2851    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2852    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2853    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2854    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2855    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2856    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2857    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2858    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2859    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2860    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2861    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2862    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2863    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2864    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2865    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2866    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2867    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2868    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2869    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2870    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2871    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2872    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2873    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2874    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2875    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2876    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2877    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2878    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2879    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2880    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2881    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2882    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2883    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2884    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2885    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2886    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2887    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2888    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2889    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2890    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2891    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2892    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2893    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2894    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2895    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2896    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2897    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2898    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2899    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2900    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2901    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2902    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2903    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2904    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2905    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2906    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2907    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2908    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2909    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2910    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2911    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2912    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2913    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2914    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2915    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2916    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2917    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2918    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2919    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2920    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2921    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2922    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2923    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2924    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2925    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2926    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2927    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2928    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2929    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2930    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2931    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2932    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2933    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2934    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2935    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2936    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2937    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2938    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2939    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2940    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2941    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2942    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2943    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2944    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2945    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2946    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2947    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2948    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2949    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2950    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2951    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2952    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2953    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2954    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2955    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2956    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2957    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2958    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2959    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2960    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2961    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2962    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2963    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2964    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2965    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2966    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2967    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2968    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2969    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2970    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2971    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2972    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2973    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2974    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2975    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2976    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2977    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2978    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2979    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2980    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2981    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2982    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2983    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2984    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2985    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2986    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2987    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2988    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2989    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2990    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2991    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2992    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2993    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2994    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2995    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2996    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2997    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2998    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 2999    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3000    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3001    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3002    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3003    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3004    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3005    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3006    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3007    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3008    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3009    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3010    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3011    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3012    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3013    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3014    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3015    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3016    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3017    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3018    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3019    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3020    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3021    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3022    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3023    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3024    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3025    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3026    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3027    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3028    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3029    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3030    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3031    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3032    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3033    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3034    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3035    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3036    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3037    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3038    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3039    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3040    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3041    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3042    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3043    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3044    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3045    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3046    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3047    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3048    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3049    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3050    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3051    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3052    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3053    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3054    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3055    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3056    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3057    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3058    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3059    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3060    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3061    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3062    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3063    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3064    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3065    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3066    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3067    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3068    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3069    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3070    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3071    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3072    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3073    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3074    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3075    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3076    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3077    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3078    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3079    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3080    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3081    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3082    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3083    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3084    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3085    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3086    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3087    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3088    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3089    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3090    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3091    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3092    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3093    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3094    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3095    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3096    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3097    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3098    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3099    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3100    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3101    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3102    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3103    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3104    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3105    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3106    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3107    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3108    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3109    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3110    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3111    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3112    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3113    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3114    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3115    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3116    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3117    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3118    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3119    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3120    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3121    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3122    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3123    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3124    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3125    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3126    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3127    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3128    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3129    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3130    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3131    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3132    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3133    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3134    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3135    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3136    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3137    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3138    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3139    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3140    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3141    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3142    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3143    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3144    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3145    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3146    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3147    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3148    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3149    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3150    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3151    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3152    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3153    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3154    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3155    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3156    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3157    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3158    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3159    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3160    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3161    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3162    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3163    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3164    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3165    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3166    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3167    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3168    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3169    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3170    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3171    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3172    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3173    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3174    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3175    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3176    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3177    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3178    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3179    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3180    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3181    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3182    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3183    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3184    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3185    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3186    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3187    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3188    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3189    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3190    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3191    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3192    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3193    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3194    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3195    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3196    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3197    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3198    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3199    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3200    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3201    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3202    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3203    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3204    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3205    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3206    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3207    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3208    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3209    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3210    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3211    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3212    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3213    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3214    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3215    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3216    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3217    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3218    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3219    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3220    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3221    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3222    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3223    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3224    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3225    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3226    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3227    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3228    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3229    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3230    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3231    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3232    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3233    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3234    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3235    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3236    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3237    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3238    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3239    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3240    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3241    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3242    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3243    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3244    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3245    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3246    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3247    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3248    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3249    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3250    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3251    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3252    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3253    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3254    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3255    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3256    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3257    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3258    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3259    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3260    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3261    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3262    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3263    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3264    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3265    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3266    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3267    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3268    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3269    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3270    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3271    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3272    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3273    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3274    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3275    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3276    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3277    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3278    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3279    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3280    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3281    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3282    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3283    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3284    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3285    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3286    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3287    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3288    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3289    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3290    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3291    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3292    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3293    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3294    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3295    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3296    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3297    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3298    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3299    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3300    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3301    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3302    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3303    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3304    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3305    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3306    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3307    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3308    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3309    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3310    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3311    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3312    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3313    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3314    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3315    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3316    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3317    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3318    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3319    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3320    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3321    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3322    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3323    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3324    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3325    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3326    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3327    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3328    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3329    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3330    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3331    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3332    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3333    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3334    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3335    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3336    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3337    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3338    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3339    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3340    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3341    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3342    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3343    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3344    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3345    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3346    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3347    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3348    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3349    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3350    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3351    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3352    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3353    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3354    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3355    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3356    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3357    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3358    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3359    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3360    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3361    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3362    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3363    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3364    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3365    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3366    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3367    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3368    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3369    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3370    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3371    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3372    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3373    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3374    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3375    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3376    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3377    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3378    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3379    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3380    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3381    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3382    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3383    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3384    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3385    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3386    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3387    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3388    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3389    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3390    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3391    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3392    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3393    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3394    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3395    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3396    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3397    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3398    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3399    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3400    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3401    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3402    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3403    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3404    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3405    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3406    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3407    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3408    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3409    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3410    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3411    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3412    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3413    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3414    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3415    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3416    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3417    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3418    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3419    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3420    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3421    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3422    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3423    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3424    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3425    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3426    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3427    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3428    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3429    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3430    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3431    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3432    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3433    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3434    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3435    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3436    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3437    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3438    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3439    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3440    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3441    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3442    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3443    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3444    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3445    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3446    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3447    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3448    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3449    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3450    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3451    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3452    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3453    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3454    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3455    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3456    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3457    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3458    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3459    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3460    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3461    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3462    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3463    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3464    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3465    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3466    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3467    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3468    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3469    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3470    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3471    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3472    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3473    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3474    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3475    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3476    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3477    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3478    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3479    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3480    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3481    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3482    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3483    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3484    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3485    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3486    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3487    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3488    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3489    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3490    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3491    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3492    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3493    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3494    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3495    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3496    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3497    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3498    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3499    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3500    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3501    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3502    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3503    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3504    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3505    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3506    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3507    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3508    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3509    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3510    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3511    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3512    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3513    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3514    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3515    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3516    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3517    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3518    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3519    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3520    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3521    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3522    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3523    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3524    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3525    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3526    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3527    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3528    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3529    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3530    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3531    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3532    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3533    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3534    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3535    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3536    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3537    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3538    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3539    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3540    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3541    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3542    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3543    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3544    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3545    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3546    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3547    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3548    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3549    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3550    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3551    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3552    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3553    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3554    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3555    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3556    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3557    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3558    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3559    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3560    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3561    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3562    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3563    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3564    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3565    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3566    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3567    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3568    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3569    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3570    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3571    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3572    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3573    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3574    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3575    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#> 3576    1200   600   600 9.572157e-06 3.198372e-07 0.7430600 0.6628171
#>       hr_ci_hi events
#> 1    0.8706544   1198
#> 2    0.8706544   1198
#> 3    0.8706544   1198
#> 4    0.8706544   1198
#> 5    0.8706544   1198
#> 6    0.8706544   1198
#> 7    0.8706544   1198
#> 8    0.8706544   1198
#> 9    0.8706544   1198
#> 10   0.8706544   1198
#> 11   0.8706544   1198
#> 12   0.8706544   1198
#> 13   0.8706544   1198
#> 14   0.8706544   1198
#> 15   0.8706544   1198
#> 16   0.8706544   1198
#> 17   0.8706544   1198
#> 18   0.8706544   1198
#> 19   0.8706544   1198
#> 20   0.8706544   1198
#> 21   0.8706544   1198
#> 22   0.8706544   1198
#> 23   0.8706544   1198
#> 24   0.8706544   1198
#> 25   0.8706544   1198
#> 26   0.8706544   1198
#> 27   0.8706544   1198
#> 28   0.8706544   1198
#> 29   0.8706544   1198
#> 30   0.8706544   1198
#> 31   0.8706544   1198
#> 32   0.8706544   1198
#> 33   0.8706544   1198
#> 34   0.8706544   1198
#> 35   0.8706544   1198
#> 36   0.8706544   1198
#> 37   0.8706544   1198
#> 38   0.8706544   1198
#> 39   0.8706544   1198
#> 40   0.8706544   1198
#> 41   0.8706544   1198
#> 42   0.8706544   1198
#> 43   0.8706544   1198
#> 44   0.8706544   1198
#> 45   0.8706544   1198
#> 46   0.8706544   1198
#> 47   0.8706544   1198
#> 48   0.8706544   1198
#> 49   0.8706544   1198
#> 50   0.8706544   1198
#> 51   0.8706544   1198
#> 52   0.8706544   1198
#> 53   0.8706544   1198
#> 54   0.8706544   1198
#> 55   0.8706544   1198
#> 56   0.8706544   1198
#> 57   0.8706544   1198
#> 58   0.8706544   1198
#> 59   0.8706544   1198
#> 60   0.8706544   1198
#> 61   0.8706544   1198
#> 62   0.8706544   1198
#> 63   0.8706544   1198
#> 64   0.8706544   1198
#> 65   0.8706544   1198
#> 66   0.8706544   1198
#> 67   0.8706544   1198
#> 68   0.8706544   1198
#> 69   0.8706544   1198
#> 70   0.8706544   1198
#> 71   0.8706544   1198
#> 72   0.8706544   1198
#> 73   0.8706544   1198
#> 74   0.8706544   1198
#> 75   0.8706544   1198
#> 76   0.8706544   1198
#> 77   0.8706544   1198
#> 78   0.8706544   1198
#> 79   0.8706544   1198
#> 80   0.8706544   1198
#> 81   0.8706544   1198
#> 82   0.8706544   1198
#> 83   0.8706544   1198
#> 84   0.8706544   1198
#> 85   0.8706544   1198
#> 86   0.8706544   1198
#> 87   0.8706544   1198
#> 88   0.8706544   1198
#> 89   0.8706544   1198
#> 90   0.8706544   1198
#> 91   0.8706544   1198
#> 92   0.8706544   1198
#> 93   0.8706544   1198
#> 94   0.8706544   1198
#> 95   0.8706544   1198
#> 96   0.8706544   1198
#> 97   0.8706544   1198
#> 98   0.8706544   1198
#> 99   0.8706544   1198
#> 100  0.8706544   1198
#> 101  0.8706544   1198
#> 102  0.8706544   1198
#> 103  0.8706544   1198
#> 104  0.8706544   1198
#> 105  0.8706544   1198
#> 106  0.8706544   1198
#> 107  0.8706544   1198
#> 108  0.8706544   1198
#> 109  0.8706544   1198
#> 110  0.8706544   1198
#> 111  0.8706544   1198
#> 112  0.8706544   1198
#> 113  0.8706544   1198
#> 114  0.8706544   1198
#> 115  0.8706544   1198
#> 116  0.8706544   1198
#> 117  0.8706544   1198
#> 118  0.8706544   1198
#> 119  0.8706544   1198
#> 120  0.8706544   1198
#> 121  0.8706544   1198
#> 122  0.8706544   1198
#> 123  0.8706544   1198
#> 124  0.8706544   1198
#> 125  0.8706544   1198
#> 126  0.8706544   1198
#> 127  0.8706544   1198
#> 128  0.8706544   1198
#> 129  0.8706544   1198
#> 130  0.8706544   1198
#> 131  0.8706544   1198
#> 132  0.8706544   1198
#> 133  0.8706544   1198
#> 134  0.8706544   1198
#> 135  0.8706544   1198
#> 136  0.8706544   1198
#> 137  0.8706544   1198
#> 138  0.8706544   1198
#> 139  0.8706544   1198
#> 140  0.8706544   1198
#> 141  0.8706544   1198
#> 142  0.8706544   1198
#> 143  0.8706544   1198
#> 144  0.8706544   1198
#> 145  0.8706544   1198
#> 146  0.8706544   1198
#> 147  0.8706544   1198
#> 148  0.8706544   1198
#> 149  0.8706544   1198
#> 150  0.8706544   1198
#> 151  0.8706544   1198
#> 152  0.8706544   1198
#> 153  0.8706544   1198
#> 154  0.8706544   1198
#> 155  0.8706544   1198
#> 156  0.8706544   1198
#> 157  0.8706544   1198
#> 158  0.8706544   1198
#> 159  0.8706544   1198
#> 160  0.8706544   1198
#> 161  0.8706544   1198
#> 162  0.8706544   1198
#> 163  0.8706544   1198
#> 164  0.8706544   1198
#> 165  0.8706544   1198
#> 166  0.8706544   1198
#> 167  0.8706544   1198
#> 168  0.8706544   1198
#> 169  0.8706544   1198
#> 170  0.8706544   1198
#> 171  0.8706544   1198
#> 172  0.8706544   1198
#> 173  0.8706544   1198
#> 174  0.8706544   1198
#> 175  0.8706544   1198
#> 176  0.8706544   1198
#> 177  0.8706544   1198
#> 178  0.8706544   1198
#> 179  0.8706544   1198
#> 180  0.8706544   1198
#> 181  0.8706544   1198
#> 182  0.8706544   1198
#> 183  0.8706544   1198
#> 184  0.8706544   1198
#> 185  0.8706544   1198
#> 186  0.8706544   1198
#> 187  0.8706544   1198
#> 188  0.8706544   1198
#> 189  0.8706544   1198
#> 190  0.8706544   1198
#> 191  0.8706544   1198
#> 192  0.8706544   1198
#> 193  0.8706544   1198
#> 194  0.8706544   1198
#> 195  0.8706544   1198
#> 196  0.8706544   1198
#> 197  0.8706544   1198
#> 198  0.8706544   1198
#> 199  0.8706544   1198
#> 200  0.8706544   1198
#> 201  0.8706544   1198
#> 202  0.8706544   1198
#> 203  0.8706544   1198
#> 204  0.8706544   1198
#> 205  0.8706544   1198
#> 206  0.8706544   1198
#> 207  0.8706544   1198
#> 208  0.8706544   1198
#> 209  0.8706544   1198
#> 210  0.8706544   1198
#> 211  0.8706544   1198
#> 212  0.8706544   1198
#> 213  0.8706544   1198
#> 214  0.8706544   1198
#> 215  0.8706544   1198
#> 216  0.8706544   1198
#> 217  0.8706544   1198
#> 218  0.8706544   1198
#> 219  0.8706544   1198
#> 220  0.8706544   1198
#> 221  0.8706544   1198
#> 222  0.8706544   1198
#> 223  0.8706544   1198
#> 224  0.8706544   1198
#> 225  0.8706544   1198
#> 226  0.8706544   1198
#> 227  0.8706544   1198
#> 228  0.8706544   1198
#> 229  0.8706544   1198
#> 230  0.8706544   1198
#> 231  0.8706544   1198
#> 232  0.8706544   1198
#> 233  0.8706544   1198
#> 234  0.8706544   1198
#> 235  0.8706544   1198
#> 236  0.8706544   1198
#> 237  0.8706544   1198
#> 238  0.8706544   1198
#> 239  0.8706544   1198
#> 240  0.8706544   1198
#> 241  0.8706544   1198
#> 242  0.8706544   1198
#> 243  0.8706544   1198
#> 244  0.8706544   1198
#> 245  0.8706544   1198
#> 246  0.8706544   1198
#> 247  0.8706544   1198
#> 248  0.8706544   1198
#> 249  0.8706544   1198
#> 250  0.8706544   1198
#> 251  0.8706544   1198
#> 252  0.8706544   1198
#> 253  0.8706544   1198
#> 254  0.8706544   1198
#> 255  0.8706544   1198
#> 256  0.8706544   1198
#> 257  0.8706544   1198
#> 258  0.8706544   1198
#> 259  0.8706544   1198
#> 260  0.8706544   1198
#> 261  0.8706544   1198
#> 262  0.8706544   1198
#> 263  0.8706544   1198
#> 264  0.8706544   1198
#> 265  0.8706544   1198
#> 266  0.8706544   1198
#> 267  0.8706544   1198
#> 268  0.8706544   1198
#> 269  0.8706544   1198
#> 270  0.8706544   1198
#> 271  0.8706544   1198
#> 272  0.8706544   1198
#> 273  0.8706544   1198
#> 274  0.8706544   1198
#> 275  0.8706544   1198
#> 276  0.8706544   1198
#> 277  0.8706544   1198
#> 278  0.8706544   1198
#> 279  0.8706544   1198
#> 280  0.8706544   1198
#> 281  0.8706544   1198
#> 282  0.8706544   1198
#> 283  0.8706544   1198
#> 284  0.8706544   1198
#> 285  0.8706544   1198
#> 286  0.8706544   1198
#> 287  0.8706544   1198
#> 288  0.8706544   1198
#> 289  0.8706544   1198
#> 290  0.8706544   1198
#> 291  0.8706544   1198
#> 292  0.8706544   1198
#> 293  0.8706544   1198
#> 294  0.8706544   1198
#> 295  0.8706544   1198
#> 296  0.8706544   1198
#> 297  0.8706544   1198
#> 298  0.8706544   1198
#> 299  0.8706544   1198
#> 300  0.8706544   1198
#> 301  0.8706544   1198
#> 302  0.8706544   1198
#> 303  0.8706544   1198
#> 304  0.8706544   1198
#> 305  0.8706544   1198
#> 306  0.8706544   1198
#> 307  0.8706544   1198
#> 308  0.8706544   1198
#> 309  0.8706544   1198
#> 310  0.8706544   1198
#> 311  0.8706544   1198
#> 312  0.8706544   1198
#> 313  0.8706544   1198
#> 314  0.8706544   1198
#> 315  0.8706544   1198
#> 316  0.8706544   1198
#> 317  0.8706544   1198
#> 318  0.8706544   1198
#> 319  0.8706544   1198
#> 320  0.8706544   1198
#> 321  0.8706544   1198
#> 322  0.8706544   1198
#> 323  0.8706544   1198
#> 324  0.8706544   1198
#> 325  0.8706544   1198
#> 326  0.8706544   1198
#> 327  0.8706544   1198
#> 328  0.8706544   1198
#> 329  0.8706544   1198
#> 330  0.8706544   1198
#> 331  0.8706544   1198
#> 332  0.8706544   1198
#> 333  0.8706544   1198
#> 334  0.8706544   1198
#> 335  0.8706544   1198
#> 336  0.8706544   1198
#> 337  0.8706544   1198
#> 338  0.8706544   1198
#> 339  0.8706544   1198
#> 340  0.8706544   1198
#> 341  0.8706544   1198
#> 342  0.8706544   1198
#> 343  0.8706544   1198
#> 344  0.8706544   1198
#> 345  0.8706544   1198
#> 346  0.8706544   1198
#> 347  0.8706544   1198
#> 348  0.8706544   1198
#> 349  0.8706544   1198
#> 350  0.8706544   1198
#> 351  0.8706544   1198
#> 352  0.8706544   1198
#> 353  0.8706544   1198
#> 354  0.8706544   1198
#> 355  0.8706544   1198
#> 356  0.8706544   1198
#> 357  0.8706544   1198
#> 358  0.8706544   1198
#> 359  0.8706544   1198
#> 360  0.8706544   1198
#> 361  0.8706544   1198
#> 362  0.8706544   1198
#> 363  0.8706544   1198
#> 364  0.8706544   1198
#> 365  0.8706544   1198
#> 366  0.8706544   1198
#> 367  0.8706544   1198
#> 368  0.8706544   1198
#> 369  0.8706544   1198
#> 370  0.8706544   1198
#> 371  0.8706544   1198
#> 372  0.8706544   1198
#> 373  0.8706544   1198
#> 374  0.8706544   1198
#> 375  0.8706544   1198
#> 376  0.8706544   1198
#> 377  0.8706544   1198
#> 378  0.8706544   1198
#> 379  0.8706544   1198
#> 380  0.8706544   1198
#> 381  0.8706544   1198
#> 382  0.8706544   1198
#> 383  0.8706544   1198
#> 384  0.8706544   1198
#> 385  0.8706544   1198
#> 386  0.8706544   1198
#> 387  0.8706544   1198
#> 388  0.8706544   1198
#> 389  0.8706544   1198
#> 390  0.8706544   1198
#> 391  0.8706544   1198
#> 392  0.8706544   1198
#> 393  0.8706544   1198
#> 394  0.8706544   1198
#> 395  0.8706544   1198
#> 396  0.8706544   1198
#> 397  0.8706544   1198
#> 398  0.8706544   1198
#> 399  0.8706544   1198
#> 400  0.8706544   1198
#> 401  0.8706544   1198
#> 402  0.8706544   1198
#> 403  0.8706544   1198
#> 404  0.8706544   1198
#> 405  0.8706544   1198
#> 406  0.8706544   1198
#> 407  0.8706544   1198
#> 408  0.8706544   1198
#> 409  0.8706544   1198
#> 410  0.8706544   1198
#> 411  0.8706544   1198
#> 412  0.8706544   1198
#> 413  0.8706544   1198
#> 414  0.8706544   1198
#> 415  0.8706544   1198
#> 416  0.8706544   1198
#> 417  0.8706544   1198
#> 418  0.8706544   1198
#> 419  0.8706544   1198
#> 420  0.8706544   1198
#> 421  0.8706544   1198
#> 422  0.8706544   1198
#> 423  0.8706544   1198
#> 424  0.8706544   1198
#> 425  0.8706544   1198
#> 426  0.8706544   1198
#> 427  0.8706544   1198
#> 428  0.8706544   1198
#> 429  0.8706544   1198
#> 430  0.8706544   1198
#> 431  0.8706544   1198
#> 432  0.8706544   1198
#> 433  0.8706544   1198
#> 434  0.8706544   1198
#> 435  0.8706544   1198
#> 436  0.8706544   1198
#> 437  0.8706544   1198
#> 438  0.8706544   1198
#> 439  0.8706544   1198
#> 440  0.8706544   1198
#> 441  0.8706544   1198
#> 442  0.8706544   1198
#> 443  0.8706544   1198
#> 444  0.8706544   1198
#> 445  0.8706544   1198
#> 446  0.8706544   1198
#> 447  0.8706544   1198
#> 448  0.8706544   1198
#> 449  0.8706544   1198
#> 450  0.8706544   1198
#> 451  0.8706544   1198
#> 452  0.8706544   1198
#> 453  0.8706544   1198
#> 454  0.8706544   1198
#> 455  0.8706544   1198
#> 456  0.8706544   1198
#> 457  0.8706544   1198
#> 458  0.8706544   1198
#> 459  0.8706544   1198
#> 460  0.8706544   1198
#> 461  0.8706544   1198
#> 462  0.8706544   1198
#> 463  0.8706544   1198
#> 464  0.8706544   1198
#> 465  0.8706544   1198
#> 466  0.8706544   1198
#> 467  0.8706544   1198
#> 468  0.8706544   1198
#> 469  0.8706544   1198
#> 470  0.8706544   1198
#> 471  0.8706544   1198
#> 472  0.8706544   1198
#> 473  0.8706544   1198
#> 474  0.8706544   1198
#> 475  0.8706544   1198
#> 476  0.8706544   1198
#> 477  0.8706544   1198
#> 478  0.8706544   1198
#> 479  0.8706544   1198
#> 480  0.8706544   1198
#> 481  0.8706544   1198
#> 482  0.8706544   1198
#> 483  0.8706544   1198
#> 484  0.8706544   1198
#> 485  0.8706544   1198
#> 486  0.8706544   1198
#> 487  0.8706544   1198
#> 488  0.8706544   1198
#> 489  0.8706544   1198
#> 490  0.8706544   1198
#> 491  0.8706544   1198
#> 492  0.8706544   1198
#> 493  0.8706544   1198
#> 494  0.8706544   1198
#> 495  0.8706544   1198
#> 496  0.8706544   1198
#> 497  0.8706544   1198
#> 498  0.8706544   1198
#> 499  0.8706544   1198
#> 500  0.8706544   1198
#> 501  0.8706544   1198
#> 502  0.8706544   1198
#> 503  0.8706544   1198
#> 504  0.8706544   1198
#> 505  0.8706544   1198
#> 506  0.8706544   1198
#> 507  0.8706544   1198
#> 508  0.8706544   1198
#> 509  0.8706544   1198
#> 510  0.8706544   1198
#> 511  0.8706544   1198
#> 512  0.8706544   1198
#> 513  0.8706544   1198
#> 514  0.8706544   1198
#> 515  0.8706544   1198
#> 516  0.8706544   1198
#> 517  0.8706544   1198
#> 518  0.8706544   1198
#> 519  0.8706544   1198
#> 520  0.8706544   1198
#> 521  0.8706544   1198
#> 522  0.8706544   1198
#> 523  0.8706544   1198
#> 524  0.8706544   1198
#> 525  0.8706544   1198
#> 526  0.8706544   1198
#> 527  0.8706544   1198
#> 528  0.8706544   1198
#> 529  0.8706544   1198
#> 530  0.8706544   1198
#> 531  0.8706544   1198
#> 532  0.8706544   1198
#> 533  0.8706544   1198
#> 534  0.8706544   1198
#> 535  0.8706544   1198
#> 536  0.8706544   1198
#> 537  0.8706544   1198
#> 538  0.8706544   1198
#> 539  0.8706544   1198
#> 540  0.8706544   1198
#> 541  0.8706544   1198
#> 542  0.8706544   1198
#> 543  0.8706544   1198
#> 544  0.8706544   1198
#> 545  0.8706544   1198
#> 546  0.8706544   1198
#> 547  0.8706544   1198
#> 548  0.8706544   1198
#> 549  0.8706544   1198
#> 550  0.8706544   1198
#> 551  0.8706544   1198
#> 552  0.8706544   1198
#> 553  0.8706544   1198
#> 554  0.8706544   1198
#> 555  0.8706544   1198
#> 556  0.8706544   1198
#> 557  0.8706544   1198
#> 558  0.8706544   1198
#> 559  0.8706544   1198
#> 560  0.8706544   1198
#> 561  0.8706544   1198
#> 562  0.8706544   1198
#> 563  0.8706544   1198
#> 564  0.8706544   1198
#> 565  0.8706544   1198
#> 566  0.8706544   1198
#> 567  0.8706544   1198
#> 568  0.8706544   1198
#> 569  0.8706544   1198
#> 570  0.8706544   1198
#> 571  0.8706544   1198
#> 572  0.8706544   1198
#> 573  0.8706544   1198
#> 574  0.8706544   1198
#> 575  0.8706544   1198
#> 576  0.8706544   1198
#> 577  0.8706544   1198
#> 578  0.8706544   1198
#> 579  0.8706544   1198
#> 580  0.8706544   1198
#> 581  0.8706544   1198
#> 582  0.8706544   1198
#> 583  0.8706544   1198
#> 584  0.8706544   1198
#> 585  0.8706544   1198
#> 586  0.8706544   1198
#> 587  0.8706544   1198
#> 588  0.8706544   1198
#> 589  0.8706544   1198
#> 590  0.8706544   1198
#> 591  0.8706544   1198
#> 592  0.8706544   1198
#> 593  0.8706544   1198
#> 594  0.8706544   1198
#> 595  0.8706544   1198
#> 596  0.8706544   1198
#> 597  0.8706544   1198
#> 598  0.8706544   1198
#> 599  0.8706544   1198
#> 600  0.8706544   1198
#> 601  0.8706544   1198
#> 602  0.8706544   1198
#> 603  0.8706544   1198
#> 604  0.8706544   1198
#> 605  0.8706544   1198
#> 606  0.8706544   1198
#> 607  0.8706544   1198
#> 608  0.8706544   1198
#> 609  0.8706544   1198
#> 610  0.8706544   1198
#> 611  0.8706544   1198
#> 612  0.8706544   1198
#> 613  0.8706544   1198
#> 614  0.8706544   1198
#> 615  0.8706544   1198
#> 616  0.8706544   1198
#> 617  0.8706544   1198
#> 618  0.8706544   1198
#> 619  0.8706544   1198
#> 620  0.8706544   1198
#> 621  0.8706544   1198
#> 622  0.8706544   1198
#> 623  0.8706544   1198
#> 624  0.8706544   1198
#> 625  0.8706544   1198
#> 626  0.8706544   1198
#> 627  0.8706544   1198
#> 628  0.8706544   1198
#> 629  0.8706544   1198
#> 630  0.8706544   1198
#> 631  0.8706544   1198
#> 632  0.8706544   1198
#> 633  0.8706544   1198
#> 634  0.8706544   1198
#> 635  0.8706544   1198
#> 636  0.8706544   1198
#> 637  0.8706544   1198
#> 638  0.8706544   1198
#> 639  0.8706544   1198
#> 640  0.8706544   1198
#> 641  0.8706544   1198
#> 642  0.8706544   1198
#> 643  0.8706544   1198
#> 644  0.8706544   1198
#> 645  0.8706544   1198
#> 646  0.8706544   1198
#> 647  0.8706544   1198
#> 648  0.8706544   1198
#> 649  0.8706544   1198
#> 650  0.8706544   1198
#> 651  0.8706544   1198
#> 652  0.8706544   1198
#> 653  0.8706544   1198
#> 654  0.8706544   1198
#> 655  0.8706544   1198
#> 656  0.8706544   1198
#> 657  0.8706544   1198
#> 658  0.8706544   1198
#> 659  0.8706544   1198
#> 660  0.8706544   1198
#> 661  0.8706544   1198
#> 662  0.8706544   1198
#> 663  0.8706544   1198
#> 664  0.8706544   1198
#> 665  0.8706544   1198
#> 666  0.8706544   1198
#> 667  0.8706544   1198
#> 668  0.8706544   1198
#> 669  0.8706544   1198
#> 670  0.8706544   1198
#> 671  0.8706544   1198
#> 672  0.8706544   1198
#> 673  0.8706544   1198
#> 674  0.8706544   1198
#> 675  0.8706544   1198
#> 676  0.8706544   1198
#> 677  0.8706544   1198
#> 678  0.8706544   1198
#> 679  0.8706544   1198
#> 680  0.8706544   1198
#> 681  0.8706544   1198
#> 682  0.8706544   1198
#> 683  0.8706544   1198
#> 684  0.8706544   1198
#> 685  0.8706544   1198
#> 686  0.8706544   1198
#> 687  0.8706544   1198
#> 688  0.8706544   1198
#> 689  0.8706544   1198
#> 690  0.8706544   1198
#> 691  0.8706544   1198
#> 692  0.8706544   1198
#> 693  0.8706544   1198
#> 694  0.8706544   1198
#> 695  0.8706544   1198
#> 696  0.8706544   1198
#> 697  0.8706544   1198
#> 698  0.8706544   1198
#> 699  0.8706544   1198
#> 700  0.8706544   1198
#> 701  0.8706544   1198
#> 702  0.8706544   1198
#> 703  0.8706544   1198
#> 704  0.8706544   1198
#> 705  0.8706544   1198
#> 706  0.8706544   1198
#> 707  0.8706544   1198
#> 708  0.8706544   1198
#> 709  0.8706544   1198
#> 710  0.8706544   1198
#> 711  0.8706544   1198
#> 712  0.8706544   1198
#> 713  0.8706544   1198
#> 714  0.8706544   1198
#> 715  0.8706544   1198
#> 716  0.8706544   1198
#> 717  0.8706544   1198
#> 718  0.8706544   1198
#> 719  0.8706544   1198
#> 720  0.8706544   1198
#> 721  0.8706544   1198
#> 722  0.8706544   1198
#> 723  0.8706544   1198
#> 724  0.8706544   1198
#> 725  0.8706544   1198
#> 726  0.8706544   1198
#> 727  0.8706544   1198
#> 728  0.8706544   1198
#> 729  0.8706544   1198
#> 730  0.8706544   1198
#> 731  0.8706544   1198
#> 732  0.8706544   1198
#> 733  0.8706544   1198
#> 734  0.8706544   1198
#> 735  0.8706544   1198
#> 736  0.8706544   1198
#> 737  0.8706544   1198
#> 738  0.8706544   1198
#> 739  0.8706544   1198
#> 740  0.8706544   1198
#> 741  0.8706544   1198
#> 742  0.8706544   1198
#> 743  0.8706544   1198
#> 744  0.8706544   1198
#> 745  0.8706544   1198
#> 746  0.8706544   1198
#> 747  0.8706544   1198
#> 748  0.8706544   1198
#> 749  0.8706544   1198
#> 750  0.8706544   1198
#> 751  0.8706544   1198
#> 752  0.8706544   1198
#> 753  0.8706544   1198
#> 754  0.8706544   1198
#> 755  0.8706544   1198
#> 756  0.8706544   1198
#> 757  0.8706544   1198
#> 758  0.8706544   1198
#> 759  0.8706544   1198
#> 760  0.8706544   1198
#> 761  0.8706544   1198
#> 762  0.8706544   1198
#> 763  0.8706544   1198
#> 764  0.8706544   1198
#> 765  0.8706544   1198
#> 766  0.8706544   1198
#> 767  0.8706544   1198
#> 768  0.8706544   1198
#> 769  0.8706544   1198
#> 770  0.8706544   1198
#> 771  0.8706544   1198
#> 772  0.8706544   1198
#> 773  0.8706544   1198
#> 774  0.8706544   1198
#> 775  0.8706544   1198
#> 776  0.8706544   1198
#> 777  0.8706544   1198
#> 778  0.8706544   1198
#> 779  0.8706544   1198
#> 780  0.8706544   1198
#> 781  0.8706544   1198
#> 782  0.8706544   1198
#> 783  0.8706544   1198
#> 784  0.8706544   1198
#> 785  0.8706544   1198
#> 786  0.8706544   1198
#> 787  0.8706544   1198
#> 788  0.8706544   1198
#> 789  0.8706544   1198
#> 790  0.8706544   1198
#> 791  0.8706544   1198
#> 792  0.8706544   1198
#> 793  0.8706544   1198
#> 794  0.8706544   1198
#> 795  0.8706544   1198
#> 796  0.8706544   1198
#> 797  0.8706544   1198
#> 798  0.8706544   1198
#> 799  0.8706544   1198
#> 800  0.8706544   1198
#> 801  0.8706544   1198
#> 802  0.8706544   1198
#> 803  0.8706544   1198
#> 804  0.8706544   1198
#> 805  0.8706544   1198
#> 806  0.8706544   1198
#> 807  0.8706544   1198
#> 808  0.8706544   1198
#> 809  0.8706544   1198
#> 810  0.8706544   1198
#> 811  0.8706544   1198
#> 812  0.8706544   1198
#> 813  0.8706544   1198
#> 814  0.8706544   1198
#> 815  0.8706544   1198
#> 816  0.8706544   1198
#> 817  0.8706544   1198
#> 818  0.8706544   1198
#> 819  0.8706544   1198
#> 820  0.8706544   1198
#> 821  0.8706544   1198
#> 822  0.8706544   1198
#> 823  0.8706544   1198
#> 824  0.8706544   1198
#> 825  0.8706544   1198
#> 826  0.8706544   1198
#> 827  0.8706544   1198
#> 828  0.8706544   1198
#> 829  0.8706544   1198
#> 830  0.8706544   1198
#> 831  0.8706544   1198
#> 832  0.8706544   1198
#> 833  0.8706544   1198
#> 834  0.8706544   1198
#> 835  0.8706544   1198
#> 836  0.8706544   1198
#> 837  0.8706544   1198
#> 838  0.8706544   1198
#> 839  0.8706544   1198
#> 840  0.8706544   1198
#> 841  0.8706544   1198
#> 842  0.8706544   1198
#> 843  0.8706544   1198
#> 844  0.8706544   1198
#> 845  0.8706544   1198
#> 846  0.8706544   1198
#> 847  0.8706544   1198
#> 848  0.8706544   1198
#> 849  0.8706544   1198
#> 850  0.8706544   1198
#> 851  0.8706544   1198
#> 852  0.8706544   1198
#> 853  0.8706544   1198
#> 854  0.8706544   1198
#> 855  0.8706544   1198
#> 856  0.8706544   1198
#> 857  0.8706544   1198
#> 858  0.8706544   1198
#> 859  0.8706544   1198
#> 860  0.8706544   1198
#> 861  0.8706544   1198
#> 862  0.8706544   1198
#> 863  0.8706544   1198
#> 864  0.8706544   1198
#> 865  0.8706544   1198
#> 866  0.8706544   1198
#> 867  0.8706544   1198
#> 868  0.8706544   1198
#> 869  0.8706544   1198
#> 870  0.8706544   1198
#> 871  0.8706544   1198
#> 872  0.8706544   1198
#> 873  0.8706544   1198
#> 874  0.8706544   1198
#> 875  0.8706544   1198
#> 876  0.8706544   1198
#> 877  0.8706544   1198
#> 878  0.8706544   1198
#> 879  0.8706544   1198
#> 880  0.8706544   1198
#> 881  0.8706544   1198
#> 882  0.8706544   1198
#> 883  0.8706544   1198
#> 884  0.8706544   1198
#> 885  0.8706544   1198
#> 886  0.8706544   1198
#> 887  0.8706544   1198
#> 888  0.8706544   1198
#> 889  0.8706544   1198
#> 890  0.8706544   1198
#> 891  0.8706544   1198
#> 892  0.8706544   1198
#> 893  0.8706544   1198
#> 894  0.8706544   1198
#> 895  0.8706544   1198
#> 896  0.8706544   1198
#> 897  0.8706544   1198
#> 898  0.8706544   1198
#> 899  0.8706544   1198
#> 900  0.8706544   1198
#> 901  0.8706544   1198
#> 902  0.8706544   1198
#> 903  0.8706544   1198
#> 904  0.8706544   1198
#> 905  0.8706544   1198
#> 906  0.8706544   1198
#> 907  0.8706544   1198
#> 908  0.8706544   1198
#> 909  0.8706544   1198
#> 910  0.8706544   1198
#> 911  0.8706544   1198
#> 912  0.8706544   1198
#> 913  0.8706544   1198
#> 914  0.8706544   1198
#> 915  0.8706544   1198
#> 916  0.8706544   1198
#> 917  0.8706544   1198
#> 918  0.8706544   1198
#> 919  0.8706544   1198
#> 920  0.8706544   1198
#> 921  0.8706544   1198
#> 922  0.8706544   1198
#> 923  0.8706544   1198
#> 924  0.8706544   1198
#> 925  0.8706544   1198
#> 926  0.8706544   1198
#> 927  0.8706544   1198
#> 928  0.8706544   1198
#> 929  0.8706544   1198
#> 930  0.8706544   1198
#> 931  0.8706544   1198
#> 932  0.8706544   1198
#> 933  0.8706544   1198
#> 934  0.8706544   1198
#> 935  0.8706544   1198
#> 936  0.8706544   1198
#> 937  0.8706544   1198
#> 938  0.8706544   1198
#> 939  0.8706544   1198
#> 940  0.8706544   1198
#> 941  0.8706544   1198
#> 942  0.8706544   1198
#> 943  0.8706544   1198
#> 944  0.8706544   1198
#> 945  0.8706544   1198
#> 946  0.8706544   1198
#> 947  0.8706544   1198
#> 948  0.8706544   1198
#> 949  0.8706544   1198
#> 950  0.8706544   1198
#> 951  0.8706544   1198
#> 952  0.8706544   1198
#> 953  0.8706544   1198
#> 954  0.8706544   1198
#> 955  0.8706544   1198
#> 956  0.8706544   1198
#> 957  0.8706544   1198
#> 958  0.8706544   1198
#> 959  0.8706544   1198
#> 960  0.8706544   1198
#> 961  0.8706544   1198
#> 962  0.8706544   1198
#> 963  0.8706544   1198
#> 964  0.8706544   1198
#> 965  0.8706544   1198
#> 966  0.8706544   1198
#> 967  0.8706544   1198
#> 968  0.8706544   1198
#> 969  0.8706544   1198
#> 970  0.8706544   1198
#> 971  0.8706544   1198
#> 972  0.8706544   1198
#> 973  0.8706544   1198
#> 974  0.8706544   1198
#> 975  0.8706544   1198
#> 976  0.8706544   1198
#> 977  0.8706544   1198
#> 978  0.8706544   1198
#> 979  0.8706544   1198
#> 980  0.8706544   1198
#> 981  0.8706544   1198
#> 982  0.8706544   1198
#> 983  0.8706544   1198
#> 984  0.8706544   1198
#> 985  0.8706544   1198
#> 986  0.8706544   1198
#> 987  0.8706544   1198
#> 988  0.8706544   1198
#> 989  0.8706544   1198
#> 990  0.8706544   1198
#> 991  0.8706544   1198
#> 992  0.8706544   1198
#> 993  0.8706544   1198
#> 994  0.8706544   1198
#> 995  0.8706544   1198
#> 996  0.8706544   1198
#> 997  0.8706544   1198
#> 998  0.8706544   1198
#> 999  0.8706544   1198
#> 1000 0.8706544   1198
#> 1001 0.8706544   1198
#> 1002 0.8706544   1198
#> 1003 0.8706544   1198
#> 1004 0.8706544   1198
#> 1005 0.8706544   1198
#> 1006 0.8706544   1198
#> 1007 0.8706544   1198
#> 1008 0.8706544   1198
#> 1009 0.8706544   1198
#> 1010 0.8706544   1198
#> 1011 0.8706544   1198
#> 1012 0.8706544   1198
#> 1013 0.8706544   1198
#> 1014 0.8706544   1198
#> 1015 0.8706544   1198
#> 1016 0.8706544   1198
#> 1017 0.8706544   1198
#> 1018 0.8706544   1198
#> 1019 0.8706544   1198
#> 1020 0.8706544   1198
#> 1021 0.8706544   1198
#> 1022 0.8706544   1198
#> 1023 0.8706544   1198
#> 1024 0.8706544   1198
#> 1025 0.8706544   1198
#> 1026 0.8706544   1198
#> 1027 0.8706544   1198
#> 1028 0.8706544   1198
#> 1029 0.8706544   1198
#> 1030 0.8706544   1198
#> 1031 0.8706544   1198
#> 1032 0.8706544   1198
#> 1033 0.8706544   1198
#> 1034 0.8706544   1198
#> 1035 0.8706544   1198
#> 1036 0.8706544   1198
#> 1037 0.8706544   1198
#> 1038 0.8706544   1198
#> 1039 0.8706544   1198
#> 1040 0.8706544   1198
#> 1041 0.8706544   1198
#> 1042 0.8706544   1198
#> 1043 0.8706544   1198
#> 1044 0.8706544   1198
#> 1045 0.8706544   1198
#> 1046 0.8706544   1198
#> 1047 0.8706544   1198
#> 1048 0.8706544   1198
#> 1049 0.8706544   1198
#> 1050 0.8706544   1198
#> 1051 0.8706544   1198
#> 1052 0.8706544   1198
#> 1053 0.8706544   1198
#> 1054 0.8706544   1198
#> 1055 0.8706544   1198
#> 1056 0.8706544   1198
#> 1057 0.8706544   1198
#> 1058 0.8706544   1198
#> 1059 0.8706544   1198
#> 1060 0.8706544   1198
#> 1061 0.8706544   1198
#> 1062 0.8706544   1198
#> 1063 0.8706544   1198
#> 1064 0.8706544   1198
#> 1065 0.8706544   1198
#> 1066 0.8706544   1198
#> 1067 0.8706544   1198
#> 1068 0.8706544   1198
#> 1069 0.8706544   1198
#> 1070 0.8706544   1198
#> 1071 0.8706544   1198
#> 1072 0.8706544   1198
#> 1073 0.8706544   1198
#> 1074 0.8706544   1198
#> 1075 0.8706544   1198
#> 1076 0.8706544   1198
#> 1077 0.8706544   1198
#> 1078 0.8706544   1198
#> 1079 0.8706544   1198
#> 1080 0.8706544   1198
#> 1081 0.8706544   1198
#> 1082 0.8706544   1198
#> 1083 0.8706544   1198
#> 1084 0.8706544   1198
#> 1085 0.8706544   1198
#> 1086 0.8706544   1198
#> 1087 0.8706544   1198
#> 1088 0.8706544   1198
#> 1089 0.8706544   1198
#> 1090 0.8706544   1198
#> 1091 0.8706544   1198
#> 1092 0.8706544   1198
#> 1093 0.8706544   1198
#> 1094 0.8706544   1198
#> 1095 0.8706544   1198
#> 1096 0.8706544   1198
#> 1097 0.8706544   1198
#> 1098 0.8706544   1198
#> 1099 0.8706544   1198
#> 1100 0.8706544   1198
#> 1101 0.8706544   1198
#> 1102 0.8706544   1198
#> 1103 0.8706544   1198
#> 1104 0.8706544   1198
#> 1105 0.8706544   1198
#> 1106 0.8706544   1198
#> 1107 0.8706544   1198
#> 1108 0.8706544   1198
#> 1109 0.8706544   1198
#> 1110 0.8706544   1198
#> 1111 0.8706544   1198
#> 1112 0.8706544   1198
#> 1113 0.8706544   1198
#> 1114 0.8706544   1198
#> 1115 0.8706544   1198
#> 1116 0.8706544   1198
#> 1117 0.8706544   1198
#> 1118 0.8706544   1198
#> 1119 0.8706544   1198
#> 1120 0.8706544   1198
#> 1121 0.8706544   1198
#> 1122 0.8706544   1198
#> 1123 0.8706544   1198
#> 1124 0.8706544   1198
#> 1125 0.8706544   1198
#> 1126 0.8706544   1198
#> 1127 0.8706544   1198
#> 1128 0.8706544   1198
#> 1129 0.8706544   1198
#> 1130 0.8706544   1198
#> 1131 0.8706544   1198
#> 1132 0.8706544   1198
#> 1133 0.8706544   1198
#> 1134 0.8706544   1198
#> 1135 0.8706544   1198
#> 1136 0.8706544   1198
#> 1137 0.8706544   1198
#> 1138 0.8706544   1198
#> 1139 0.8706544   1198
#> 1140 0.8706544   1198
#> 1141 0.8706544   1198
#> 1142 0.8706544   1198
#> 1143 0.8706544   1198
#> 1144 0.8706544   1198
#> 1145 0.8706544   1198
#> 1146 0.8706544   1198
#> 1147 0.8706544   1198
#> 1148 0.8706544   1198
#> 1149 0.8706544   1198
#> 1150 0.8706544   1198
#> 1151 0.8706544   1198
#> 1152 0.8706544   1198
#> 1153 0.8706544   1198
#> 1154 0.8706544   1198
#> 1155 0.8706544   1198
#> 1156 0.8706544   1198
#> 1157 0.8706544   1198
#> 1158 0.8706544   1198
#> 1159 0.8706544   1198
#> 1160 0.8706544   1198
#> 1161 0.8706544   1198
#> 1162 0.8706544   1198
#> 1163 0.8706544   1198
#> 1164 0.8706544   1198
#> 1165 0.8706544   1198
#> 1166 0.8706544   1198
#> 1167 0.8706544   1198
#> 1168 0.8706544   1198
#> 1169 0.8706544   1198
#> 1170 0.8706544   1198
#> 1171 0.8706544   1198
#> 1172 0.8706544   1198
#> 1173 0.8706544   1198
#> 1174 0.8706544   1198
#> 1175 0.8706544   1198
#> 1176 0.8706544   1198
#> 1177 0.8706544   1198
#> 1178 0.8706544   1198
#> 1179 0.8706544   1198
#> 1180 0.8706544   1198
#> 1181 0.8706544   1198
#> 1182 0.8706544   1198
#> 1183 0.8706544   1198
#> 1184 0.8706544   1198
#> 1185 0.8706544   1198
#> 1186 0.8706544   1198
#> 1187 0.8706544   1198
#> 1188 0.8706544   1198
#> 1189 0.8706544   1198
#> 1190 0.8706544   1198
#> 1191 0.8706544   1198
#> 1192 0.8706544   1198
#> 1193 0.8706544   1198
#> 1194 0.7792467   1200
#> 1195 0.7792467   1200
#> 1196 0.7792467   1200
#> 1197 0.7792467   1200
#> 1198 0.7792467   1200
#> 1199 0.7792467   1200
#> 1200 0.7792467   1200
#> 1201 0.7792467   1200
#> 1202 0.7792467   1200
#> 1203 0.7792467   1200
#> 1204 0.7792467   1200
#> 1205 0.7792467   1200
#> 1206 0.7792467   1200
#> 1207 0.7792467   1200
#> 1208 0.7792467   1200
#> 1209 0.7792467   1200
#> 1210 0.7792467   1200
#> 1211 0.7792467   1200
#> 1212 0.7792467   1200
#> 1213 0.7792467   1200
#> 1214 0.7792467   1200
#> 1215 0.7792467   1200
#> 1216 0.7792467   1200
#> 1217 0.7792467   1200
#> 1218 0.7792467   1200
#> 1219 0.7792467   1200
#> 1220 0.7792467   1200
#> 1221 0.7792467   1200
#> 1222 0.7792467   1200
#> 1223 0.7792467   1200
#> 1224 0.7792467   1200
#> 1225 0.7792467   1200
#> 1226 0.7792467   1200
#> 1227 0.7792467   1200
#> 1228 0.7792467   1200
#> 1229 0.7792467   1200
#> 1230 0.7792467   1200
#> 1231 0.7792467   1200
#> 1232 0.7792467   1200
#> 1233 0.7792467   1200
#> 1234 0.7792467   1200
#> 1235 0.7792467   1200
#> 1236 0.7792467   1200
#> 1237 0.7792467   1200
#> 1238 0.7792467   1200
#> 1239 0.7792467   1200
#> 1240 0.7792467   1200
#> 1241 0.7792467   1200
#> 1242 0.7792467   1200
#> 1243 0.7792467   1200
#> 1244 0.7792467   1200
#> 1245 0.7792467   1200
#> 1246 0.7792467   1200
#> 1247 0.7792467   1200
#> 1248 0.7792467   1200
#> 1249 0.7792467   1200
#> 1250 0.7792467   1200
#> 1251 0.7792467   1200
#> 1252 0.7792467   1200
#> 1253 0.7792467   1200
#> 1254 0.7792467   1200
#> 1255 0.7792467   1200
#> 1256 0.7792467   1200
#> 1257 0.7792467   1200
#> 1258 0.7792467   1200
#> 1259 0.7792467   1200
#> 1260 0.7792467   1200
#> 1261 0.7792467   1200
#> 1262 0.7792467   1200
#> 1263 0.7792467   1200
#> 1264 0.7792467   1200
#> 1265 0.7792467   1200
#> 1266 0.7792467   1200
#> 1267 0.7792467   1200
#> 1268 0.7792467   1200
#> 1269 0.7792467   1200
#> 1270 0.7792467   1200
#> 1271 0.7792467   1200
#> 1272 0.7792467   1200
#> 1273 0.7792467   1200
#> 1274 0.7792467   1200
#> 1275 0.7792467   1200
#> 1276 0.7792467   1200
#> 1277 0.7792467   1200
#> 1278 0.7792467   1200
#> 1279 0.7792467   1200
#> 1280 0.7792467   1200
#> 1281 0.7792467   1200
#> 1282 0.7792467   1200
#> 1283 0.7792467   1200
#> 1284 0.7792467   1200
#> 1285 0.7792467   1200
#> 1286 0.7792467   1200
#> 1287 0.7792467   1200
#> 1288 0.7792467   1200
#> 1289 0.7792467   1200
#> 1290 0.7792467   1200
#> 1291 0.7792467   1200
#> 1292 0.7792467   1200
#> 1293 0.7792467   1200
#> 1294 0.7792467   1200
#> 1295 0.7792467   1200
#> 1296 0.7792467   1200
#> 1297 0.7792467   1200
#> 1298 0.7792467   1200
#> 1299 0.7792467   1200
#> 1300 0.7792467   1200
#> 1301 0.7792467   1200
#> 1302 0.7792467   1200
#> 1303 0.7792467   1200
#> 1304 0.7792467   1200
#> 1305 0.7792467   1200
#> 1306 0.7792467   1200
#> 1307 0.7792467   1200
#> 1308 0.7792467   1200
#> 1309 0.7792467   1200
#> 1310 0.7792467   1200
#> 1311 0.7792467   1200
#> 1312 0.7792467   1200
#> 1313 0.7792467   1200
#> 1314 0.7792467   1200
#> 1315 0.7792467   1200
#> 1316 0.7792467   1200
#> 1317 0.7792467   1200
#> 1318 0.7792467   1200
#> 1319 0.7792467   1200
#> 1320 0.7792467   1200
#> 1321 0.7792467   1200
#> 1322 0.7792467   1200
#> 1323 0.7792467   1200
#> 1324 0.7792467   1200
#> 1325 0.7792467   1200
#> 1326 0.7792467   1200
#> 1327 0.7792467   1200
#> 1328 0.7792467   1200
#> 1329 0.7792467   1200
#> 1330 0.7792467   1200
#> 1331 0.7792467   1200
#> 1332 0.7792467   1200
#> 1333 0.7792467   1200
#> 1334 0.7792467   1200
#> 1335 0.7792467   1200
#> 1336 0.7792467   1200
#> 1337 0.7792467   1200
#> 1338 0.7792467   1200
#> 1339 0.7792467   1200
#> 1340 0.7792467   1200
#> 1341 0.7792467   1200
#> 1342 0.7792467   1200
#> 1343 0.7792467   1200
#> 1344 0.7792467   1200
#> 1345 0.7792467   1200
#> 1346 0.7792467   1200
#> 1347 0.7792467   1200
#> 1348 0.7792467   1200
#> 1349 0.7792467   1200
#> 1350 0.7792467   1200
#> 1351 0.7792467   1200
#> 1352 0.7792467   1200
#> 1353 0.7792467   1200
#> 1354 0.7792467   1200
#> 1355 0.7792467   1200
#> 1356 0.7792467   1200
#> 1357 0.7792467   1200
#> 1358 0.7792467   1200
#> 1359 0.7792467   1200
#> 1360 0.7792467   1200
#> 1361 0.7792467   1200
#> 1362 0.7792467   1200
#> 1363 0.7792467   1200
#> 1364 0.7792467   1200
#> 1365 0.7792467   1200
#> 1366 0.7792467   1200
#> 1367 0.7792467   1200
#> 1368 0.7792467   1200
#> 1369 0.7792467   1200
#> 1370 0.7792467   1200
#> 1371 0.7792467   1200
#> 1372 0.7792467   1200
#> 1373 0.7792467   1200
#> 1374 0.7792467   1200
#> 1375 0.7792467   1200
#> 1376 0.7792467   1200
#> 1377 0.7792467   1200
#> 1378 0.7792467   1200
#> 1379 0.7792467   1200
#> 1380 0.7792467   1200
#> 1381 0.7792467   1200
#> 1382 0.7792467   1200
#> 1383 0.7792467   1200
#> 1384 0.7792467   1200
#> 1385 0.7792467   1200
#> 1386 0.7792467   1200
#> 1387 0.7792467   1200
#> 1388 0.7792467   1200
#> 1389 0.7792467   1200
#> 1390 0.7792467   1200
#> 1391 0.7792467   1200
#> 1392 0.7792467   1200
#> 1393 0.7792467   1200
#> 1394 0.7792467   1200
#> 1395 0.7792467   1200
#> 1396 0.7792467   1200
#> 1397 0.7792467   1200
#> 1398 0.7792467   1200
#> 1399 0.7792467   1200
#> 1400 0.7792467   1200
#> 1401 0.7792467   1200
#> 1402 0.7792467   1200
#> 1403 0.7792467   1200
#> 1404 0.7792467   1200
#> 1405 0.7792467   1200
#> 1406 0.7792467   1200
#> 1407 0.7792467   1200
#> 1408 0.7792467   1200
#> 1409 0.7792467   1200
#> 1410 0.7792467   1200
#> 1411 0.7792467   1200
#> 1412 0.7792467   1200
#> 1413 0.7792467   1200
#> 1414 0.7792467   1200
#> 1415 0.7792467   1200
#> 1416 0.7792467   1200
#> 1417 0.7792467   1200
#> 1418 0.7792467   1200
#> 1419 0.7792467   1200
#> 1420 0.7792467   1200
#> 1421 0.7792467   1200
#> 1422 0.7792467   1200
#> 1423 0.7792467   1200
#> 1424 0.7792467   1200
#> 1425 0.7792467   1200
#> 1426 0.7792467   1200
#> 1427 0.7792467   1200
#> 1428 0.7792467   1200
#> 1429 0.7792467   1200
#> 1430 0.7792467   1200
#> 1431 0.7792467   1200
#> 1432 0.7792467   1200
#> 1433 0.7792467   1200
#> 1434 0.7792467   1200
#> 1435 0.7792467   1200
#> 1436 0.7792467   1200
#> 1437 0.7792467   1200
#> 1438 0.7792467   1200
#> 1439 0.7792467   1200
#> 1440 0.7792467   1200
#> 1441 0.7792467   1200
#> 1442 0.7792467   1200
#> 1443 0.7792467   1200
#> 1444 0.7792467   1200
#> 1445 0.7792467   1200
#> 1446 0.7792467   1200
#> 1447 0.7792467   1200
#> 1448 0.7792467   1200
#> 1449 0.7792467   1200
#> 1450 0.7792467   1200
#> 1451 0.7792467   1200
#> 1452 0.7792467   1200
#> 1453 0.7792467   1200
#> 1454 0.7792467   1200
#> 1455 0.7792467   1200
#> 1456 0.7792467   1200
#> 1457 0.7792467   1200
#> 1458 0.7792467   1200
#> 1459 0.7792467   1200
#> 1460 0.7792467   1200
#> 1461 0.7792467   1200
#> 1462 0.7792467   1200
#> 1463 0.7792467   1200
#> 1464 0.7792467   1200
#> 1465 0.7792467   1200
#> 1466 0.7792467   1200
#> 1467 0.7792467   1200
#> 1468 0.7792467   1200
#> 1469 0.7792467   1200
#> 1470 0.7792467   1200
#> 1471 0.7792467   1200
#> 1472 0.7792467   1200
#> 1473 0.7792467   1200
#> 1474 0.7792467   1200
#> 1475 0.7792467   1200
#> 1476 0.7792467   1200
#> 1477 0.7792467   1200
#> 1478 0.7792467   1200
#> 1479 0.7792467   1200
#> 1480 0.7792467   1200
#> 1481 0.7792467   1200
#> 1482 0.7792467   1200
#> 1483 0.7792467   1200
#> 1484 0.7792467   1200
#> 1485 0.7792467   1200
#> 1486 0.7792467   1200
#> 1487 0.7792467   1200
#> 1488 0.7792467   1200
#> 1489 0.7792467   1200
#> 1490 0.7792467   1200
#> 1491 0.7792467   1200
#> 1492 0.7792467   1200
#> 1493 0.7792467   1200
#> 1494 0.7792467   1200
#> 1495 0.7792467   1200
#> 1496 0.7792467   1200
#> 1497 0.7792467   1200
#> 1498 0.7792467   1200
#> 1499 0.7792467   1200
#> 1500 0.7792467   1200
#> 1501 0.7792467   1200
#> 1502 0.7792467   1200
#> 1503 0.7792467   1200
#> 1504 0.7792467   1200
#> 1505 0.7792467   1200
#> 1506 0.7792467   1200
#> 1507 0.7792467   1200
#> 1508 0.7792467   1200
#> 1509 0.7792467   1200
#> 1510 0.7792467   1200
#> 1511 0.7792467   1200
#> 1512 0.7792467   1200
#> 1513 0.7792467   1200
#> 1514 0.7792467   1200
#> 1515 0.7792467   1200
#> 1516 0.7792467   1200
#> 1517 0.7792467   1200
#> 1518 0.7792467   1200
#> 1519 0.7792467   1200
#> 1520 0.7792467   1200
#> 1521 0.7792467   1200
#> 1522 0.7792467   1200
#> 1523 0.7792467   1200
#> 1524 0.7792467   1200
#> 1525 0.7792467   1200
#> 1526 0.7792467   1200
#> 1527 0.7792467   1200
#> 1528 0.7792467   1200
#> 1529 0.7792467   1200
#> 1530 0.7792467   1200
#> 1531 0.7792467   1200
#> 1532 0.7792467   1200
#> 1533 0.7792467   1200
#> 1534 0.7792467   1200
#> 1535 0.7792467   1200
#> 1536 0.7792467   1200
#> 1537 0.7792467   1200
#> 1538 0.7792467   1200
#> 1539 0.7792467   1200
#> 1540 0.7792467   1200
#> 1541 0.7792467   1200
#> 1542 0.7792467   1200
#> 1543 0.7792467   1200
#> 1544 0.7792467   1200
#> 1545 0.7792467   1200
#> 1546 0.7792467   1200
#> 1547 0.7792467   1200
#> 1548 0.7792467   1200
#> 1549 0.7792467   1200
#> 1550 0.7792467   1200
#> 1551 0.7792467   1200
#> 1552 0.7792467   1200
#> 1553 0.7792467   1200
#> 1554 0.7792467   1200
#> 1555 0.7792467   1200
#> 1556 0.7792467   1200
#> 1557 0.7792467   1200
#> 1558 0.7792467   1200
#> 1559 0.7792467   1200
#> 1560 0.7792467   1200
#> 1561 0.7792467   1200
#> 1562 0.7792467   1200
#> 1563 0.7792467   1200
#> 1564 0.7792467   1200
#> 1565 0.7792467   1200
#> 1566 0.7792467   1200
#> 1567 0.7792467   1200
#> 1568 0.7792467   1200
#> 1569 0.7792467   1200
#> 1570 0.7792467   1200
#> 1571 0.7792467   1200
#> 1572 0.7792467   1200
#> 1573 0.7792467   1200
#> 1574 0.7792467   1200
#> 1575 0.7792467   1200
#> 1576 0.7792467   1200
#> 1577 0.7792467   1200
#> 1578 0.7792467   1200
#> 1579 0.7792467   1200
#> 1580 0.7792467   1200
#> 1581 0.7792467   1200
#> 1582 0.7792467   1200
#> 1583 0.7792467   1200
#> 1584 0.7792467   1200
#> 1585 0.7792467   1200
#> 1586 0.7792467   1200
#> 1587 0.7792467   1200
#> 1588 0.7792467   1200
#> 1589 0.7792467   1200
#> 1590 0.7792467   1200
#> 1591 0.7792467   1200
#> 1592 0.7792467   1200
#> 1593 0.7792467   1200
#> 1594 0.7792467   1200
#> 1595 0.7792467   1200
#> 1596 0.7792467   1200
#> 1597 0.7792467   1200
#> 1598 0.7792467   1200
#> 1599 0.7792467   1200
#> 1600 0.7792467   1200
#> 1601 0.7792467   1200
#> 1602 0.7792467   1200
#> 1603 0.7792467   1200
#> 1604 0.7792467   1200
#> 1605 0.7792467   1200
#> 1606 0.7792467   1200
#> 1607 0.7792467   1200
#> 1608 0.7792467   1200
#> 1609 0.7792467   1200
#> 1610 0.7792467   1200
#> 1611 0.7792467   1200
#> 1612 0.7792467   1200
#> 1613 0.7792467   1200
#> 1614 0.7792467   1200
#> 1615 0.7792467   1200
#> 1616 0.7792467   1200
#> 1617 0.7792467   1200
#> 1618 0.7792467   1200
#> 1619 0.7792467   1200
#> 1620 0.7792467   1200
#> 1621 0.7792467   1200
#> 1622 0.7792467   1200
#> 1623 0.7792467   1200
#> 1624 0.7792467   1200
#> 1625 0.7792467   1200
#> 1626 0.7792467   1200
#> 1627 0.7792467   1200
#> 1628 0.7792467   1200
#> 1629 0.7792467   1200
#> 1630 0.7792467   1200
#> 1631 0.7792467   1200
#> 1632 0.7792467   1200
#> 1633 0.7792467   1200
#> 1634 0.7792467   1200
#> 1635 0.7792467   1200
#> 1636 0.7792467   1200
#> 1637 0.7792467   1200
#> 1638 0.7792467   1200
#> 1639 0.7792467   1200
#> 1640 0.7792467   1200
#> 1641 0.7792467   1200
#> 1642 0.7792467   1200
#> 1643 0.7792467   1200
#> 1644 0.7792467   1200
#> 1645 0.7792467   1200
#> 1646 0.7792467   1200
#> 1647 0.7792467   1200
#> 1648 0.7792467   1200
#> 1649 0.7792467   1200
#> 1650 0.7792467   1200
#> 1651 0.7792467   1200
#> 1652 0.7792467   1200
#> 1653 0.7792467   1200
#> 1654 0.7792467   1200
#> 1655 0.7792467   1200
#> 1656 0.7792467   1200
#> 1657 0.7792467   1200
#> 1658 0.7792467   1200
#> 1659 0.7792467   1200
#> 1660 0.7792467   1200
#> 1661 0.7792467   1200
#> 1662 0.7792467   1200
#> 1663 0.7792467   1200
#> 1664 0.7792467   1200
#> 1665 0.7792467   1200
#> 1666 0.7792467   1200
#> 1667 0.7792467   1200
#> 1668 0.7792467   1200
#> 1669 0.7792467   1200
#> 1670 0.7792467   1200
#> 1671 0.7792467   1200
#> 1672 0.7792467   1200
#> 1673 0.7792467   1200
#> 1674 0.7792467   1200
#> 1675 0.7792467   1200
#> 1676 0.7792467   1200
#> 1677 0.7792467   1200
#> 1678 0.7792467   1200
#> 1679 0.7792467   1200
#> 1680 0.7792467   1200
#> 1681 0.7792467   1200
#> 1682 0.7792467   1200
#> 1683 0.7792467   1200
#> 1684 0.7792467   1200
#> 1685 0.7792467   1200
#> 1686 0.7792467   1200
#> 1687 0.7792467   1200
#> 1688 0.7792467   1200
#> 1689 0.7792467   1200
#> 1690 0.7792467   1200
#> 1691 0.7792467   1200
#> 1692 0.7792467   1200
#> 1693 0.7792467   1200
#> 1694 0.7792467   1200
#> 1695 0.7792467   1200
#> 1696 0.7792467   1200
#> 1697 0.7792467   1200
#> 1698 0.7792467   1200
#> 1699 0.7792467   1200
#> 1700 0.7792467   1200
#> 1701 0.7792467   1200
#> 1702 0.7792467   1200
#> 1703 0.7792467   1200
#> 1704 0.7792467   1200
#> 1705 0.7792467   1200
#> 1706 0.7792467   1200
#> 1707 0.7792467   1200
#> 1708 0.7792467   1200
#> 1709 0.7792467   1200
#> 1710 0.7792467   1200
#> 1711 0.7792467   1200
#> 1712 0.7792467   1200
#> 1713 0.7792467   1200
#> 1714 0.7792467   1200
#> 1715 0.7792467   1200
#> 1716 0.7792467   1200
#> 1717 0.7792467   1200
#> 1718 0.7792467   1200
#> 1719 0.7792467   1200
#> 1720 0.7792467   1200
#> 1721 0.7792467   1200
#> 1722 0.7792467   1200
#> 1723 0.7792467   1200
#> 1724 0.7792467   1200
#> 1725 0.7792467   1200
#> 1726 0.7792467   1200
#> 1727 0.7792467   1200
#> 1728 0.7792467   1200
#> 1729 0.7792467   1200
#> 1730 0.7792467   1200
#> 1731 0.7792467   1200
#> 1732 0.7792467   1200
#> 1733 0.7792467   1200
#> 1734 0.7792467   1200
#> 1735 0.7792467   1200
#> 1736 0.7792467   1200
#> 1737 0.7792467   1200
#> 1738 0.7792467   1200
#> 1739 0.7792467   1200
#> 1740 0.7792467   1200
#> 1741 0.7792467   1200
#> 1742 0.7792467   1200
#> 1743 0.7792467   1200
#> 1744 0.7792467   1200
#> 1745 0.7792467   1200
#> 1746 0.7792467   1200
#> 1747 0.7792467   1200
#> 1748 0.7792467   1200
#> 1749 0.7792467   1200
#> 1750 0.7792467   1200
#> 1751 0.7792467   1200
#> 1752 0.7792467   1200
#> 1753 0.7792467   1200
#> 1754 0.7792467   1200
#> 1755 0.7792467   1200
#> 1756 0.7792467   1200
#> 1757 0.7792467   1200
#> 1758 0.7792467   1200
#> 1759 0.7792467   1200
#> 1760 0.7792467   1200
#> 1761 0.7792467   1200
#> 1762 0.7792467   1200
#> 1763 0.7792467   1200
#> 1764 0.7792467   1200
#> 1765 0.7792467   1200
#> 1766 0.7792467   1200
#> 1767 0.7792467   1200
#> 1768 0.7792467   1200
#> 1769 0.7792467   1200
#> 1770 0.7792467   1200
#> 1771 0.7792467   1200
#> 1772 0.7792467   1200
#> 1773 0.7792467   1200
#> 1774 0.7792467   1200
#> 1775 0.7792467   1200
#> 1776 0.7792467   1200
#> 1777 0.7792467   1200
#> 1778 0.7792467   1200
#> 1779 0.7792467   1200
#> 1780 0.7792467   1200
#> 1781 0.7792467   1200
#> 1782 0.7792467   1200
#> 1783 0.7792467   1200
#> 1784 0.7792467   1200
#> 1785 0.7792467   1200
#> 1786 0.7792467   1200
#> 1787 0.7792467   1200
#> 1788 0.7792467   1200
#> 1789 0.7792467   1200
#> 1790 0.7792467   1200
#> 1791 0.7792467   1200
#> 1792 0.7792467   1200
#> 1793 0.7792467   1200
#> 1794 0.7792467   1200
#> 1795 0.7792467   1200
#> 1796 0.7792467   1200
#> 1797 0.7792467   1200
#> 1798 0.7792467   1200
#> 1799 0.7792467   1200
#> 1800 0.7792467   1200
#> 1801 0.7792467   1200
#> 1802 0.7792467   1200
#> 1803 0.7792467   1200
#> 1804 0.7792467   1200
#> 1805 0.7792467   1200
#> 1806 0.7792467   1200
#> 1807 0.7792467   1200
#> 1808 0.7792467   1200
#> 1809 0.7792467   1200
#> 1810 0.7792467   1200
#> 1811 0.7792467   1200
#> 1812 0.7792467   1200
#> 1813 0.7792467   1200
#> 1814 0.7792467   1200
#> 1815 0.7792467   1200
#> 1816 0.7792467   1200
#> 1817 0.7792467   1200
#> 1818 0.7792467   1200
#> 1819 0.7792467   1200
#> 1820 0.7792467   1200
#> 1821 0.7792467   1200
#> 1822 0.7792467   1200
#> 1823 0.7792467   1200
#> 1824 0.7792467   1200
#> 1825 0.7792467   1200
#> 1826 0.7792467   1200
#> 1827 0.7792467   1200
#> 1828 0.7792467   1200
#> 1829 0.7792467   1200
#> 1830 0.7792467   1200
#> 1831 0.7792467   1200
#> 1832 0.7792467   1200
#> 1833 0.7792467   1200
#> 1834 0.7792467   1200
#> 1835 0.7792467   1200
#> 1836 0.7792467   1200
#> 1837 0.7792467   1200
#> 1838 0.7792467   1200
#> 1839 0.7792467   1200
#> 1840 0.7792467   1200
#> 1841 0.7792467   1200
#> 1842 0.7792467   1200
#> 1843 0.7792467   1200
#> 1844 0.7792467   1200
#> 1845 0.7792467   1200
#> 1846 0.7792467   1200
#> 1847 0.7792467   1200
#> 1848 0.7792467   1200
#> 1849 0.7792467   1200
#> 1850 0.7792467   1200
#> 1851 0.7792467   1200
#> 1852 0.7792467   1200
#> 1853 0.7792467   1200
#> 1854 0.7792467   1200
#> 1855 0.7792467   1200
#> 1856 0.7792467   1200
#> 1857 0.7792467   1200
#> 1858 0.7792467   1200
#> 1859 0.7792467   1200
#> 1860 0.7792467   1200
#> 1861 0.7792467   1200
#> 1862 0.7792467   1200
#> 1863 0.7792467   1200
#> 1864 0.7792467   1200
#> 1865 0.7792467   1200
#> 1866 0.7792467   1200
#> 1867 0.7792467   1200
#> 1868 0.7792467   1200
#> 1869 0.7792467   1200
#> 1870 0.7792467   1200
#> 1871 0.7792467   1200
#> 1872 0.7792467   1200
#> 1873 0.7792467   1200
#> 1874 0.7792467   1200
#> 1875 0.7792467   1200
#> 1876 0.7792467   1200
#> 1877 0.7792467   1200
#> 1878 0.7792467   1200
#> 1879 0.7792467   1200
#> 1880 0.7792467   1200
#> 1881 0.7792467   1200
#> 1882 0.7792467   1200
#> 1883 0.7792467   1200
#> 1884 0.7792467   1200
#> 1885 0.7792467   1200
#> 1886 0.7792467   1200
#> 1887 0.7792467   1200
#> 1888 0.7792467   1200
#> 1889 0.7792467   1200
#> 1890 0.7792467   1200
#> 1891 0.7792467   1200
#> 1892 0.7792467   1200
#> 1893 0.7792467   1200
#> 1894 0.7792467   1200
#> 1895 0.7792467   1200
#> 1896 0.7792467   1200
#> 1897 0.7792467   1200
#> 1898 0.7792467   1200
#> 1899 0.7792467   1200
#> 1900 0.7792467   1200
#> 1901 0.7792467   1200
#> 1902 0.7792467   1200
#> 1903 0.7792467   1200
#> 1904 0.7792467   1200
#> 1905 0.7792467   1200
#> 1906 0.7792467   1200
#> 1907 0.7792467   1200
#> 1908 0.7792467   1200
#> 1909 0.7792467   1200
#> 1910 0.7792467   1200
#> 1911 0.7792467   1200
#> 1912 0.7792467   1200
#> 1913 0.7792467   1200
#> 1914 0.7792467   1200
#> 1915 0.7792467   1200
#> 1916 0.7792467   1200
#> 1917 0.7792467   1200
#> 1918 0.7792467   1200
#> 1919 0.7792467   1200
#> 1920 0.7792467   1200
#> 1921 0.7792467   1200
#> 1922 0.7792467   1200
#> 1923 0.7792467   1200
#> 1924 0.7792467   1200
#> 1925 0.7792467   1200
#> 1926 0.7792467   1200
#> 1927 0.7792467   1200
#> 1928 0.7792467   1200
#> 1929 0.7792467   1200
#> 1930 0.7792467   1200
#> 1931 0.7792467   1200
#> 1932 0.7792467   1200
#> 1933 0.7792467   1200
#> 1934 0.7792467   1200
#> 1935 0.7792467   1200
#> 1936 0.7792467   1200
#> 1937 0.7792467   1200
#> 1938 0.7792467   1200
#> 1939 0.7792467   1200
#> 1940 0.7792467   1200
#> 1941 0.7792467   1200
#> 1942 0.7792467   1200
#> 1943 0.7792467   1200
#> 1944 0.7792467   1200
#> 1945 0.7792467   1200
#> 1946 0.7792467   1200
#> 1947 0.7792467   1200
#> 1948 0.7792467   1200
#> 1949 0.7792467   1200
#> 1950 0.7792467   1200
#> 1951 0.7792467   1200
#> 1952 0.7792467   1200
#> 1953 0.7792467   1200
#> 1954 0.7792467   1200
#> 1955 0.7792467   1200
#> 1956 0.7792467   1200
#> 1957 0.7792467   1200
#> 1958 0.7792467   1200
#> 1959 0.7792467   1200
#> 1960 0.7792467   1200
#> 1961 0.7792467   1200
#> 1962 0.7792467   1200
#> 1963 0.7792467   1200
#> 1964 0.7792467   1200
#> 1965 0.7792467   1200
#> 1966 0.7792467   1200
#> 1967 0.7792467   1200
#> 1968 0.7792467   1200
#> 1969 0.7792467   1200
#> 1970 0.7792467   1200
#> 1971 0.7792467   1200
#> 1972 0.7792467   1200
#> 1973 0.7792467   1200
#> 1974 0.7792467   1200
#> 1975 0.7792467   1200
#> 1976 0.7792467   1200
#> 1977 0.7792467   1200
#> 1978 0.7792467   1200
#> 1979 0.7792467   1200
#> 1980 0.7792467   1200
#> 1981 0.7792467   1200
#> 1982 0.7792467   1200
#> 1983 0.7792467   1200
#> 1984 0.7792467   1200
#> 1985 0.7792467   1200
#> 1986 0.7792467   1200
#> 1987 0.7792467   1200
#> 1988 0.7792467   1200
#> 1989 0.7792467   1200
#> 1990 0.7792467   1200
#> 1991 0.7792467   1200
#> 1992 0.7792467   1200
#> 1993 0.7792467   1200
#> 1994 0.7792467   1200
#> 1995 0.7792467   1200
#> 1996 0.7792467   1200
#> 1997 0.7792467   1200
#> 1998 0.7792467   1200
#> 1999 0.7792467   1200
#> 2000 0.7792467   1200
#> 2001 0.7792467   1200
#> 2002 0.7792467   1200
#> 2003 0.7792467   1200
#> 2004 0.7792467   1200
#> 2005 0.7792467   1200
#> 2006 0.7792467   1200
#> 2007 0.7792467   1200
#> 2008 0.7792467   1200
#> 2009 0.7792467   1200
#> 2010 0.7792467   1200
#> 2011 0.7792467   1200
#> 2012 0.7792467   1200
#> 2013 0.7792467   1200
#> 2014 0.7792467   1200
#> 2015 0.7792467   1200
#> 2016 0.7792467   1200
#> 2017 0.7792467   1200
#> 2018 0.7792467   1200
#> 2019 0.7792467   1200
#> 2020 0.7792467   1200
#> 2021 0.7792467   1200
#> 2022 0.7792467   1200
#> 2023 0.7792467   1200
#> 2024 0.7792467   1200
#> 2025 0.7792467   1200
#> 2026 0.7792467   1200
#> 2027 0.7792467   1200
#> 2028 0.7792467   1200
#> 2029 0.7792467   1200
#> 2030 0.7792467   1200
#> 2031 0.7792467   1200
#> 2032 0.7792467   1200
#> 2033 0.7792467   1200
#> 2034 0.7792467   1200
#> 2035 0.7792467   1200
#> 2036 0.7792467   1200
#> 2037 0.7792467   1200
#> 2038 0.7792467   1200
#> 2039 0.7792467   1200
#> 2040 0.7792467   1200
#> 2041 0.7792467   1200
#> 2042 0.7792467   1200
#> 2043 0.7792467   1200
#> 2044 0.7792467   1200
#> 2045 0.7792467   1200
#> 2046 0.7792467   1200
#> 2047 0.7792467   1200
#> 2048 0.7792467   1200
#> 2049 0.7792467   1200
#> 2050 0.7792467   1200
#> 2051 0.7792467   1200
#> 2052 0.7792467   1200
#> 2053 0.7792467   1200
#> 2054 0.7792467   1200
#> 2055 0.7792467   1200
#> 2056 0.7792467   1200
#> 2057 0.7792467   1200
#> 2058 0.7792467   1200
#> 2059 0.7792467   1200
#> 2060 0.7792467   1200
#> 2061 0.7792467   1200
#> 2062 0.7792467   1200
#> 2063 0.7792467   1200
#> 2064 0.7792467   1200
#> 2065 0.7792467   1200
#> 2066 0.7792467   1200
#> 2067 0.7792467   1200
#> 2068 0.7792467   1200
#> 2069 0.7792467   1200
#> 2070 0.7792467   1200
#> 2071 0.7792467   1200
#> 2072 0.7792467   1200
#> 2073 0.7792467   1200
#> 2074 0.7792467   1200
#> 2075 0.7792467   1200
#> 2076 0.7792467   1200
#> 2077 0.7792467   1200
#> 2078 0.7792467   1200
#> 2079 0.7792467   1200
#> 2080 0.7792467   1200
#> 2081 0.7792467   1200
#> 2082 0.7792467   1200
#> 2083 0.7792467   1200
#> 2084 0.7792467   1200
#> 2085 0.7792467   1200
#> 2086 0.7792467   1200
#> 2087 0.7792467   1200
#> 2088 0.7792467   1200
#> 2089 0.7792467   1200
#> 2090 0.7792467   1200
#> 2091 0.7792467   1200
#> 2092 0.7792467   1200
#> 2093 0.7792467   1200
#> 2094 0.7792467   1200
#> 2095 0.7792467   1200
#> 2096 0.7792467   1200
#> 2097 0.7792467   1200
#> 2098 0.7792467   1200
#> 2099 0.7792467   1200
#> 2100 0.7792467   1200
#> 2101 0.7792467   1200
#> 2102 0.7792467   1200
#> 2103 0.7792467   1200
#> 2104 0.7792467   1200
#> 2105 0.7792467   1200
#> 2106 0.7792467   1200
#> 2107 0.7792467   1200
#> 2108 0.7792467   1200
#> 2109 0.7792467   1200
#> 2110 0.7792467   1200
#> 2111 0.7792467   1200
#> 2112 0.7792467   1200
#> 2113 0.7792467   1200
#> 2114 0.7792467   1200
#> 2115 0.7792467   1200
#> 2116 0.7792467   1200
#> 2117 0.7792467   1200
#> 2118 0.7792467   1200
#> 2119 0.7792467   1200
#> 2120 0.7792467   1200
#> 2121 0.7792467   1200
#> 2122 0.7792467   1200
#> 2123 0.7792467   1200
#> 2124 0.7792467   1200
#> 2125 0.7792467   1200
#> 2126 0.7792467   1200
#> 2127 0.7792467   1200
#> 2128 0.7792467   1200
#> 2129 0.7792467   1200
#> 2130 0.7792467   1200
#> 2131 0.7792467   1200
#> 2132 0.7792467   1200
#> 2133 0.7792467   1200
#> 2134 0.7792467   1200
#> 2135 0.7792467   1200
#> 2136 0.7792467   1200
#> 2137 0.7792467   1200
#> 2138 0.7792467   1200
#> 2139 0.7792467   1200
#> 2140 0.7792467   1200
#> 2141 0.7792467   1200
#> 2142 0.7792467   1200
#> 2143 0.7792467   1200
#> 2144 0.7792467   1200
#> 2145 0.7792467   1200
#> 2146 0.7792467   1200
#> 2147 0.7792467   1200
#> 2148 0.7792467   1200
#> 2149 0.7792467   1200
#> 2150 0.7792467   1200
#> 2151 0.7792467   1200
#> 2152 0.7792467   1200
#> 2153 0.7792467   1200
#> 2154 0.7792467   1200
#> 2155 0.7792467   1200
#> 2156 0.7792467   1200
#> 2157 0.7792467   1200
#> 2158 0.7792467   1200
#> 2159 0.7792467   1200
#> 2160 0.7792467   1200
#> 2161 0.7792467   1200
#> 2162 0.7792467   1200
#> 2163 0.7792467   1200
#> 2164 0.7792467   1200
#> 2165 0.7792467   1200
#> 2166 0.7792467   1200
#> 2167 0.7792467   1200
#> 2168 0.7792467   1200
#> 2169 0.7792467   1200
#> 2170 0.7792467   1200
#> 2171 0.7792467   1200
#> 2172 0.7792467   1200
#> 2173 0.7792467   1200
#> 2174 0.7792467   1200
#> 2175 0.7792467   1200
#> 2176 0.7792467   1200
#> 2177 0.7792467   1200
#> 2178 0.7792467   1200
#> 2179 0.7792467   1200
#> 2180 0.7792467   1200
#> 2181 0.7792467   1200
#> 2182 0.7792467   1200
#> 2183 0.7792467   1200
#> 2184 0.7792467   1200
#> 2185 0.7792467   1200
#> 2186 0.7792467   1200
#> 2187 0.7792467   1200
#> 2188 0.7792467   1200
#> 2189 0.7792467   1200
#> 2190 0.7792467   1200
#> 2191 0.7792467   1200
#> 2192 0.7792467   1200
#> 2193 0.7792467   1200
#> 2194 0.7792467   1200
#> 2195 0.7792467   1200
#> 2196 0.7792467   1200
#> 2197 0.7792467   1200
#> 2198 0.7792467   1200
#> 2199 0.7792467   1200
#> 2200 0.7792467   1200
#> 2201 0.7792467   1200
#> 2202 0.7792467   1200
#> 2203 0.7792467   1200
#> 2204 0.7792467   1200
#> 2205 0.7792467   1200
#> 2206 0.7792467   1200
#> 2207 0.7792467   1200
#> 2208 0.7792467   1200
#> 2209 0.7792467   1200
#> 2210 0.7792467   1200
#> 2211 0.7792467   1200
#> 2212 0.7792467   1200
#> 2213 0.7792467   1200
#> 2214 0.7792467   1200
#> 2215 0.7792467   1200
#> 2216 0.7792467   1200
#> 2217 0.7792467   1200
#> 2218 0.7792467   1200
#> 2219 0.7792467   1200
#> 2220 0.7792467   1200
#> 2221 0.7792467   1200
#> 2222 0.7792467   1200
#> 2223 0.7792467   1200
#> 2224 0.7792467   1200
#> 2225 0.7792467   1200
#> 2226 0.7792467   1200
#> 2227 0.7792467   1200
#> 2228 0.7792467   1200
#> 2229 0.7792467   1200
#> 2230 0.7792467   1200
#> 2231 0.7792467   1200
#> 2232 0.7792467   1200
#> 2233 0.7792467   1200
#> 2234 0.7792467   1200
#> 2235 0.7792467   1200
#> 2236 0.7792467   1200
#> 2237 0.7792467   1200
#> 2238 0.7792467   1200
#> 2239 0.7792467   1200
#> 2240 0.7792467   1200
#> 2241 0.7792467   1200
#> 2242 0.7792467   1200
#> 2243 0.7792467   1200
#> 2244 0.7792467   1200
#> 2245 0.7792467   1200
#> 2246 0.7792467   1200
#> 2247 0.7792467   1200
#> 2248 0.7792467   1200
#> 2249 0.7792467   1200
#> 2250 0.7792467   1200
#> 2251 0.7792467   1200
#> 2252 0.7792467   1200
#> 2253 0.7792467   1200
#> 2254 0.7792467   1200
#> 2255 0.7792467   1200
#> 2256 0.7792467   1200
#> 2257 0.7792467   1200
#> 2258 0.7792467   1200
#> 2259 0.7792467   1200
#> 2260 0.7792467   1200
#> 2261 0.7792467   1200
#> 2262 0.7792467   1200
#> 2263 0.7792467   1200
#> 2264 0.7792467   1200
#> 2265 0.7792467   1200
#> 2266 0.7792467   1200
#> 2267 0.7792467   1200
#> 2268 0.7792467   1200
#> 2269 0.7792467   1200
#> 2270 0.7792467   1200
#> 2271 0.7792467   1200
#> 2272 0.7792467   1200
#> 2273 0.7792467   1200
#> 2274 0.7792467   1200
#> 2275 0.7792467   1200
#> 2276 0.7792467   1200
#> 2277 0.7792467   1200
#> 2278 0.7792467   1200
#> 2279 0.7792467   1200
#> 2280 0.7792467   1200
#> 2281 0.7792467   1200
#> 2282 0.7792467   1200
#> 2283 0.7792467   1200
#> 2284 0.7792467   1200
#> 2285 0.7792467   1200
#> 2286 0.7792467   1200
#> 2287 0.7792467   1200
#> 2288 0.7792467   1200
#> 2289 0.7792467   1200
#> 2290 0.7792467   1200
#> 2291 0.7792467   1200
#> 2292 0.7792467   1200
#> 2293 0.7792467   1200
#> 2294 0.7792467   1200
#> 2295 0.7792467   1200
#> 2296 0.7792467   1200
#> 2297 0.7792467   1200
#> 2298 0.7792467   1200
#> 2299 0.7792467   1200
#> 2300 0.7792467   1200
#> 2301 0.7792467   1200
#> 2302 0.7792467   1200
#> 2303 0.7792467   1200
#> 2304 0.7792467   1200
#> 2305 0.7792467   1200
#> 2306 0.7792467   1200
#> 2307 0.7792467   1200
#> 2308 0.7792467   1200
#> 2309 0.7792467   1200
#> 2310 0.7792467   1200
#> 2311 0.7792467   1200
#> 2312 0.7792467   1200
#> 2313 0.7792467   1200
#> 2314 0.7792467   1200
#> 2315 0.7792467   1200
#> 2316 0.7792467   1200
#> 2317 0.7792467   1200
#> 2318 0.7792467   1200
#> 2319 0.7792467   1200
#> 2320 0.7792467   1200
#> 2321 0.7792467   1200
#> 2322 0.7792467   1200
#> 2323 0.7792467   1200
#> 2324 0.7792467   1200
#> 2325 0.7792467   1200
#> 2326 0.7792467   1200
#> 2327 0.7792467   1200
#> 2328 0.7792467   1200
#> 2329 0.7792467   1200
#> 2330 0.7792467   1200
#> 2331 0.7792467   1200
#> 2332 0.7792467   1200
#> 2333 0.7792467   1200
#> 2334 0.7792467   1200
#> 2335 0.7792467   1200
#> 2336 0.7792467   1200
#> 2337 0.7792467   1200
#> 2338 0.7792467   1200
#> 2339 0.7792467   1200
#> 2340 0.7792467   1200
#> 2341 0.7792467   1200
#> 2342 0.7792467   1200
#> 2343 0.7792467   1200
#> 2344 0.7792467   1200
#> 2345 0.7792467   1200
#> 2346 0.7792467   1200
#> 2347 0.7792467   1200
#> 2348 0.7792467   1200
#> 2349 0.7792467   1200
#> 2350 0.7792467   1200
#> 2351 0.7792467   1200
#> 2352 0.7792467   1200
#> 2353 0.7792467   1200
#> 2354 0.7792467   1200
#> 2355 0.7792467   1200
#> 2356 0.7792467   1200
#> 2357 0.7792467   1200
#> 2358 0.7792467   1200
#> 2359 0.7792467   1200
#> 2360 0.7792467   1200
#> 2361 0.7792467   1200
#> 2362 0.7792467   1200
#> 2363 0.7792467   1200
#> 2364 0.7792467   1200
#> 2365 0.7792467   1200
#> 2366 0.7792467   1200
#> 2367 0.7792467   1200
#> 2368 0.7792467   1200
#> 2369 0.7792467   1200
#> 2370 0.7792467   1200
#> 2371 0.7792467   1200
#> 2372 0.7792467   1200
#> 2373 0.7792467   1200
#> 2374 0.7792467   1200
#> 2375 0.7792467   1200
#> 2376 0.7792467   1200
#> 2377 0.7792467   1200
#> 2378 0.7792467   1200
#> 2379 0.7792467   1200
#> 2380 0.7792467   1200
#> 2381 0.7792467   1200
#> 2382 0.7792467   1200
#> 2383 0.7792467   1200
#> 2384 0.8330175   1199
#> 2385 0.8330175   1199
#> 2386 0.8330175   1199
#> 2387 0.8330175   1199
#> 2388 0.8330175   1199
#> 2389 0.8330175   1199
#> 2390 0.8330175   1199
#> 2391 0.8330175   1199
#> 2392 0.8330175   1199
#> 2393 0.8330175   1199
#> 2394 0.8330175   1199
#> 2395 0.8330175   1199
#> 2396 0.8330175   1199
#> 2397 0.8330175   1199
#> 2398 0.8330175   1199
#> 2399 0.8330175   1199
#> 2400 0.8330175   1199
#> 2401 0.8330175   1199
#> 2402 0.8330175   1199
#> 2403 0.8330175   1199
#> 2404 0.8330175   1199
#> 2405 0.8330175   1199
#> 2406 0.8330175   1199
#> 2407 0.8330175   1199
#> 2408 0.8330175   1199
#> 2409 0.8330175   1199
#> 2410 0.8330175   1199
#> 2411 0.8330175   1199
#> 2412 0.8330175   1199
#> 2413 0.8330175   1199
#> 2414 0.8330175   1199
#> 2415 0.8330175   1199
#> 2416 0.8330175   1199
#> 2417 0.8330175   1199
#> 2418 0.8330175   1199
#> 2419 0.8330175   1199
#> 2420 0.8330175   1199
#> 2421 0.8330175   1199
#> 2422 0.8330175   1199
#> 2423 0.8330175   1199
#> 2424 0.8330175   1199
#> 2425 0.8330175   1199
#> 2426 0.8330175   1199
#> 2427 0.8330175   1199
#> 2428 0.8330175   1199
#> 2429 0.8330175   1199
#> 2430 0.8330175   1199
#> 2431 0.8330175   1199
#> 2432 0.8330175   1199
#> 2433 0.8330175   1199
#> 2434 0.8330175   1199
#> 2435 0.8330175   1199
#> 2436 0.8330175   1199
#> 2437 0.8330175   1199
#> 2438 0.8330175   1199
#> 2439 0.8330175   1199
#> 2440 0.8330175   1199
#> 2441 0.8330175   1199
#> 2442 0.8330175   1199
#> 2443 0.8330175   1199
#> 2444 0.8330175   1199
#> 2445 0.8330175   1199
#> 2446 0.8330175   1199
#> 2447 0.8330175   1199
#> 2448 0.8330175   1199
#> 2449 0.8330175   1199
#> 2450 0.8330175   1199
#> 2451 0.8330175   1199
#> 2452 0.8330175   1199
#> 2453 0.8330175   1199
#> 2454 0.8330175   1199
#> 2455 0.8330175   1199
#> 2456 0.8330175   1199
#> 2457 0.8330175   1199
#> 2458 0.8330175   1199
#> 2459 0.8330175   1199
#> 2460 0.8330175   1199
#> 2461 0.8330175   1199
#> 2462 0.8330175   1199
#> 2463 0.8330175   1199
#> 2464 0.8330175   1199
#> 2465 0.8330175   1199
#> 2466 0.8330175   1199
#> 2467 0.8330175   1199
#> 2468 0.8330175   1199
#> 2469 0.8330175   1199
#> 2470 0.8330175   1199
#> 2471 0.8330175   1199
#> 2472 0.8330175   1199
#> 2473 0.8330175   1199
#> 2474 0.8330175   1199
#> 2475 0.8330175   1199
#> 2476 0.8330175   1199
#> 2477 0.8330175   1199
#> 2478 0.8330175   1199
#> 2479 0.8330175   1199
#> 2480 0.8330175   1199
#> 2481 0.8330175   1199
#> 2482 0.8330175   1199
#> 2483 0.8330175   1199
#> 2484 0.8330175   1199
#> 2485 0.8330175   1199
#> 2486 0.8330175   1199
#> 2487 0.8330175   1199
#> 2488 0.8330175   1199
#> 2489 0.8330175   1199
#> 2490 0.8330175   1199
#> 2491 0.8330175   1199
#> 2492 0.8330175   1199
#> 2493 0.8330175   1199
#> 2494 0.8330175   1199
#> 2495 0.8330175   1199
#> 2496 0.8330175   1199
#> 2497 0.8330175   1199
#> 2498 0.8330175   1199
#> 2499 0.8330175   1199
#> 2500 0.8330175   1199
#> 2501 0.8330175   1199
#> 2502 0.8330175   1199
#> 2503 0.8330175   1199
#> 2504 0.8330175   1199
#> 2505 0.8330175   1199
#> 2506 0.8330175   1199
#> 2507 0.8330175   1199
#> 2508 0.8330175   1199
#> 2509 0.8330175   1199
#> 2510 0.8330175   1199
#> 2511 0.8330175   1199
#> 2512 0.8330175   1199
#> 2513 0.8330175   1199
#> 2514 0.8330175   1199
#> 2515 0.8330175   1199
#> 2516 0.8330175   1199
#> 2517 0.8330175   1199
#> 2518 0.8330175   1199
#> 2519 0.8330175   1199
#> 2520 0.8330175   1199
#> 2521 0.8330175   1199
#> 2522 0.8330175   1199
#> 2523 0.8330175   1199
#> 2524 0.8330175   1199
#> 2525 0.8330175   1199
#> 2526 0.8330175   1199
#> 2527 0.8330175   1199
#> 2528 0.8330175   1199
#> 2529 0.8330175   1199
#> 2530 0.8330175   1199
#> 2531 0.8330175   1199
#> 2532 0.8330175   1199
#> 2533 0.8330175   1199
#> 2534 0.8330175   1199
#> 2535 0.8330175   1199
#> 2536 0.8330175   1199
#> 2537 0.8330175   1199
#> 2538 0.8330175   1199
#> 2539 0.8330175   1199
#> 2540 0.8330175   1199
#> 2541 0.8330175   1199
#> 2542 0.8330175   1199
#> 2543 0.8330175   1199
#> 2544 0.8330175   1199
#> 2545 0.8330175   1199
#> 2546 0.8330175   1199
#> 2547 0.8330175   1199
#> 2548 0.8330175   1199
#> 2549 0.8330175   1199
#> 2550 0.8330175   1199
#> 2551 0.8330175   1199
#> 2552 0.8330175   1199
#> 2553 0.8330175   1199
#> 2554 0.8330175   1199
#> 2555 0.8330175   1199
#> 2556 0.8330175   1199
#> 2557 0.8330175   1199
#> 2558 0.8330175   1199
#> 2559 0.8330175   1199
#> 2560 0.8330175   1199
#> 2561 0.8330175   1199
#> 2562 0.8330175   1199
#> 2563 0.8330175   1199
#> 2564 0.8330175   1199
#> 2565 0.8330175   1199
#> 2566 0.8330175   1199
#> 2567 0.8330175   1199
#> 2568 0.8330175   1199
#> 2569 0.8330175   1199
#> 2570 0.8330175   1199
#> 2571 0.8330175   1199
#> 2572 0.8330175   1199
#> 2573 0.8330175   1199
#> 2574 0.8330175   1199
#> 2575 0.8330175   1199
#> 2576 0.8330175   1199
#> 2577 0.8330175   1199
#> 2578 0.8330175   1199
#> 2579 0.8330175   1199
#> 2580 0.8330175   1199
#> 2581 0.8330175   1199
#> 2582 0.8330175   1199
#> 2583 0.8330175   1199
#> 2584 0.8330175   1199
#> 2585 0.8330175   1199
#> 2586 0.8330175   1199
#> 2587 0.8330175   1199
#> 2588 0.8330175   1199
#> 2589 0.8330175   1199
#> 2590 0.8330175   1199
#> 2591 0.8330175   1199
#> 2592 0.8330175   1199
#> 2593 0.8330175   1199
#> 2594 0.8330175   1199
#> 2595 0.8330175   1199
#> 2596 0.8330175   1199
#> 2597 0.8330175   1199
#> 2598 0.8330175   1199
#> 2599 0.8330175   1199
#> 2600 0.8330175   1199
#> 2601 0.8330175   1199
#> 2602 0.8330175   1199
#> 2603 0.8330175   1199
#> 2604 0.8330175   1199
#> 2605 0.8330175   1199
#> 2606 0.8330175   1199
#> 2607 0.8330175   1199
#> 2608 0.8330175   1199
#> 2609 0.8330175   1199
#> 2610 0.8330175   1199
#> 2611 0.8330175   1199
#> 2612 0.8330175   1199
#> 2613 0.8330175   1199
#> 2614 0.8330175   1199
#> 2615 0.8330175   1199
#> 2616 0.8330175   1199
#> 2617 0.8330175   1199
#> 2618 0.8330175   1199
#> 2619 0.8330175   1199
#> 2620 0.8330175   1199
#> 2621 0.8330175   1199
#> 2622 0.8330175   1199
#> 2623 0.8330175   1199
#> 2624 0.8330175   1199
#> 2625 0.8330175   1199
#> 2626 0.8330175   1199
#> 2627 0.8330175   1199
#> 2628 0.8330175   1199
#> 2629 0.8330175   1199
#> 2630 0.8330175   1199
#> 2631 0.8330175   1199
#> 2632 0.8330175   1199
#> 2633 0.8330175   1199
#> 2634 0.8330175   1199
#> 2635 0.8330175   1199
#> 2636 0.8330175   1199
#> 2637 0.8330175   1199
#> 2638 0.8330175   1199
#> 2639 0.8330175   1199
#> 2640 0.8330175   1199
#> 2641 0.8330175   1199
#> 2642 0.8330175   1199
#> 2643 0.8330175   1199
#> 2644 0.8330175   1199
#> 2645 0.8330175   1199
#> 2646 0.8330175   1199
#> 2647 0.8330175   1199
#> 2648 0.8330175   1199
#> 2649 0.8330175   1199
#> 2650 0.8330175   1199
#> 2651 0.8330175   1199
#> 2652 0.8330175   1199
#> 2653 0.8330175   1199
#> 2654 0.8330175   1199
#> 2655 0.8330175   1199
#> 2656 0.8330175   1199
#> 2657 0.8330175   1199
#> 2658 0.8330175   1199
#> 2659 0.8330175   1199
#> 2660 0.8330175   1199
#> 2661 0.8330175   1199
#> 2662 0.8330175   1199
#> 2663 0.8330175   1199
#> 2664 0.8330175   1199
#> 2665 0.8330175   1199
#> 2666 0.8330175   1199
#> 2667 0.8330175   1199
#> 2668 0.8330175   1199
#> 2669 0.8330175   1199
#> 2670 0.8330175   1199
#> 2671 0.8330175   1199
#> 2672 0.8330175   1199
#> 2673 0.8330175   1199
#> 2674 0.8330175   1199
#> 2675 0.8330175   1199
#> 2676 0.8330175   1199
#> 2677 0.8330175   1199
#> 2678 0.8330175   1199
#> 2679 0.8330175   1199
#> 2680 0.8330175   1199
#> 2681 0.8330175   1199
#> 2682 0.8330175   1199
#> 2683 0.8330175   1199
#> 2684 0.8330175   1199
#> 2685 0.8330175   1199
#> 2686 0.8330175   1199
#> 2687 0.8330175   1199
#> 2688 0.8330175   1199
#> 2689 0.8330175   1199
#> 2690 0.8330175   1199
#> 2691 0.8330175   1199
#> 2692 0.8330175   1199
#> 2693 0.8330175   1199
#> 2694 0.8330175   1199
#> 2695 0.8330175   1199
#> 2696 0.8330175   1199
#> 2697 0.8330175   1199
#> 2698 0.8330175   1199
#> 2699 0.8330175   1199
#> 2700 0.8330175   1199
#> 2701 0.8330175   1199
#> 2702 0.8330175   1199
#> 2703 0.8330175   1199
#> 2704 0.8330175   1199
#> 2705 0.8330175   1199
#> 2706 0.8330175   1199
#> 2707 0.8330175   1199
#> 2708 0.8330175   1199
#> 2709 0.8330175   1199
#> 2710 0.8330175   1199
#> 2711 0.8330175   1199
#> 2712 0.8330175   1199
#> 2713 0.8330175   1199
#> 2714 0.8330175   1199
#> 2715 0.8330175   1199
#> 2716 0.8330175   1199
#> 2717 0.8330175   1199
#> 2718 0.8330175   1199
#> 2719 0.8330175   1199
#> 2720 0.8330175   1199
#> 2721 0.8330175   1199
#> 2722 0.8330175   1199
#> 2723 0.8330175   1199
#> 2724 0.8330175   1199
#> 2725 0.8330175   1199
#> 2726 0.8330175   1199
#> 2727 0.8330175   1199
#> 2728 0.8330175   1199
#> 2729 0.8330175   1199
#> 2730 0.8330175   1199
#> 2731 0.8330175   1199
#> 2732 0.8330175   1199
#> 2733 0.8330175   1199
#> 2734 0.8330175   1199
#> 2735 0.8330175   1199
#> 2736 0.8330175   1199
#> 2737 0.8330175   1199
#> 2738 0.8330175   1199
#> 2739 0.8330175   1199
#> 2740 0.8330175   1199
#> 2741 0.8330175   1199
#> 2742 0.8330175   1199
#> 2743 0.8330175   1199
#> 2744 0.8330175   1199
#> 2745 0.8330175   1199
#> 2746 0.8330175   1199
#> 2747 0.8330175   1199
#> 2748 0.8330175   1199
#> 2749 0.8330175   1199
#> 2750 0.8330175   1199
#> 2751 0.8330175   1199
#> 2752 0.8330175   1199
#> 2753 0.8330175   1199
#> 2754 0.8330175   1199
#> 2755 0.8330175   1199
#> 2756 0.8330175   1199
#> 2757 0.8330175   1199
#> 2758 0.8330175   1199
#> 2759 0.8330175   1199
#> 2760 0.8330175   1199
#> 2761 0.8330175   1199
#> 2762 0.8330175   1199
#> 2763 0.8330175   1199
#> 2764 0.8330175   1199
#> 2765 0.8330175   1199
#> 2766 0.8330175   1199
#> 2767 0.8330175   1199
#> 2768 0.8330175   1199
#> 2769 0.8330175   1199
#> 2770 0.8330175   1199
#> 2771 0.8330175   1199
#> 2772 0.8330175   1199
#> 2773 0.8330175   1199
#> 2774 0.8330175   1199
#> 2775 0.8330175   1199
#> 2776 0.8330175   1199
#> 2777 0.8330175   1199
#> 2778 0.8330175   1199
#> 2779 0.8330175   1199
#> 2780 0.8330175   1199
#> 2781 0.8330175   1199
#> 2782 0.8330175   1199
#> 2783 0.8330175   1199
#> 2784 0.8330175   1199
#> 2785 0.8330175   1199
#> 2786 0.8330175   1199
#> 2787 0.8330175   1199
#> 2788 0.8330175   1199
#> 2789 0.8330175   1199
#> 2790 0.8330175   1199
#> 2791 0.8330175   1199
#> 2792 0.8330175   1199
#> 2793 0.8330175   1199
#> 2794 0.8330175   1199
#> 2795 0.8330175   1199
#> 2796 0.8330175   1199
#> 2797 0.8330175   1199
#> 2798 0.8330175   1199
#> 2799 0.8330175   1199
#> 2800 0.8330175   1199
#> 2801 0.8330175   1199
#> 2802 0.8330175   1199
#> 2803 0.8330175   1199
#> 2804 0.8330175   1199
#> 2805 0.8330175   1199
#> 2806 0.8330175   1199
#> 2807 0.8330175   1199
#> 2808 0.8330175   1199
#> 2809 0.8330175   1199
#> 2810 0.8330175   1199
#> 2811 0.8330175   1199
#> 2812 0.8330175   1199
#> 2813 0.8330175   1199
#> 2814 0.8330175   1199
#> 2815 0.8330175   1199
#> 2816 0.8330175   1199
#> 2817 0.8330175   1199
#> 2818 0.8330175   1199
#> 2819 0.8330175   1199
#> 2820 0.8330175   1199
#> 2821 0.8330175   1199
#> 2822 0.8330175   1199
#> 2823 0.8330175   1199
#> 2824 0.8330175   1199
#> 2825 0.8330175   1199
#> 2826 0.8330175   1199
#> 2827 0.8330175   1199
#> 2828 0.8330175   1199
#> 2829 0.8330175   1199
#> 2830 0.8330175   1199
#> 2831 0.8330175   1199
#> 2832 0.8330175   1199
#> 2833 0.8330175   1199
#> 2834 0.8330175   1199
#> 2835 0.8330175   1199
#> 2836 0.8330175   1199
#> 2837 0.8330175   1199
#> 2838 0.8330175   1199
#> 2839 0.8330175   1199
#> 2840 0.8330175   1199
#> 2841 0.8330175   1199
#> 2842 0.8330175   1199
#> 2843 0.8330175   1199
#> 2844 0.8330175   1199
#> 2845 0.8330175   1199
#> 2846 0.8330175   1199
#> 2847 0.8330175   1199
#> 2848 0.8330175   1199
#> 2849 0.8330175   1199
#> 2850 0.8330175   1199
#> 2851 0.8330175   1199
#> 2852 0.8330175   1199
#> 2853 0.8330175   1199
#> 2854 0.8330175   1199
#> 2855 0.8330175   1199
#> 2856 0.8330175   1199
#> 2857 0.8330175   1199
#> 2858 0.8330175   1199
#> 2859 0.8330175   1199
#> 2860 0.8330175   1199
#> 2861 0.8330175   1199
#> 2862 0.8330175   1199
#> 2863 0.8330175   1199
#> 2864 0.8330175   1199
#> 2865 0.8330175   1199
#> 2866 0.8330175   1199
#> 2867 0.8330175   1199
#> 2868 0.8330175   1199
#> 2869 0.8330175   1199
#> 2870 0.8330175   1199
#> 2871 0.8330175   1199
#> 2872 0.8330175   1199
#> 2873 0.8330175   1199
#> 2874 0.8330175   1199
#> 2875 0.8330175   1199
#> 2876 0.8330175   1199
#> 2877 0.8330175   1199
#> 2878 0.8330175   1199
#> 2879 0.8330175   1199
#> 2880 0.8330175   1199
#> 2881 0.8330175   1199
#> 2882 0.8330175   1199
#> 2883 0.8330175   1199
#> 2884 0.8330175   1199
#> 2885 0.8330175   1199
#> 2886 0.8330175   1199
#> 2887 0.8330175   1199
#> 2888 0.8330175   1199
#> 2889 0.8330175   1199
#> 2890 0.8330175   1199
#> 2891 0.8330175   1199
#> 2892 0.8330175   1199
#> 2893 0.8330175   1199
#> 2894 0.8330175   1199
#> 2895 0.8330175   1199
#> 2896 0.8330175   1199
#> 2897 0.8330175   1199
#> 2898 0.8330175   1199
#> 2899 0.8330175   1199
#> 2900 0.8330175   1199
#> 2901 0.8330175   1199
#> 2902 0.8330175   1199
#> 2903 0.8330175   1199
#> 2904 0.8330175   1199
#> 2905 0.8330175   1199
#> 2906 0.8330175   1199
#> 2907 0.8330175   1199
#> 2908 0.8330175   1199
#> 2909 0.8330175   1199
#> 2910 0.8330175   1199
#> 2911 0.8330175   1199
#> 2912 0.8330175   1199
#> 2913 0.8330175   1199
#> 2914 0.8330175   1199
#> 2915 0.8330175   1199
#> 2916 0.8330175   1199
#> 2917 0.8330175   1199
#> 2918 0.8330175   1199
#> 2919 0.8330175   1199
#> 2920 0.8330175   1199
#> 2921 0.8330175   1199
#> 2922 0.8330175   1199
#> 2923 0.8330175   1199
#> 2924 0.8330175   1199
#> 2925 0.8330175   1199
#> 2926 0.8330175   1199
#> 2927 0.8330175   1199
#> 2928 0.8330175   1199
#> 2929 0.8330175   1199
#> 2930 0.8330175   1199
#> 2931 0.8330175   1199
#> 2932 0.8330175   1199
#> 2933 0.8330175   1199
#> 2934 0.8330175   1199
#> 2935 0.8330175   1199
#> 2936 0.8330175   1199
#> 2937 0.8330175   1199
#> 2938 0.8330175   1199
#> 2939 0.8330175   1199
#> 2940 0.8330175   1199
#> 2941 0.8330175   1199
#> 2942 0.8330175   1199
#> 2943 0.8330175   1199
#> 2944 0.8330175   1199
#> 2945 0.8330175   1199
#> 2946 0.8330175   1199
#> 2947 0.8330175   1199
#> 2948 0.8330175   1199
#> 2949 0.8330175   1199
#> 2950 0.8330175   1199
#> 2951 0.8330175   1199
#> 2952 0.8330175   1199
#> 2953 0.8330175   1199
#> 2954 0.8330175   1199
#> 2955 0.8330175   1199
#> 2956 0.8330175   1199
#> 2957 0.8330175   1199
#> 2958 0.8330175   1199
#> 2959 0.8330175   1199
#> 2960 0.8330175   1199
#> 2961 0.8330175   1199
#> 2962 0.8330175   1199
#> 2963 0.8330175   1199
#> 2964 0.8330175   1199
#> 2965 0.8330175   1199
#> 2966 0.8330175   1199
#> 2967 0.8330175   1199
#> 2968 0.8330175   1199
#> 2969 0.8330175   1199
#> 2970 0.8330175   1199
#> 2971 0.8330175   1199
#> 2972 0.8330175   1199
#> 2973 0.8330175   1199
#> 2974 0.8330175   1199
#> 2975 0.8330175   1199
#> 2976 0.8330175   1199
#> 2977 0.8330175   1199
#> 2978 0.8330175   1199
#> 2979 0.8330175   1199
#> 2980 0.8330175   1199
#> 2981 0.8330175   1199
#> 2982 0.8330175   1199
#> 2983 0.8330175   1199
#> 2984 0.8330175   1199
#> 2985 0.8330175   1199
#> 2986 0.8330175   1199
#> 2987 0.8330175   1199
#> 2988 0.8330175   1199
#> 2989 0.8330175   1199
#> 2990 0.8330175   1199
#> 2991 0.8330175   1199
#> 2992 0.8330175   1199
#> 2993 0.8330175   1199
#> 2994 0.8330175   1199
#> 2995 0.8330175   1199
#> 2996 0.8330175   1199
#> 2997 0.8330175   1199
#> 2998 0.8330175   1199
#> 2999 0.8330175   1199
#> 3000 0.8330175   1199
#> 3001 0.8330175   1199
#> 3002 0.8330175   1199
#> 3003 0.8330175   1199
#> 3004 0.8330175   1199
#> 3005 0.8330175   1199
#> 3006 0.8330175   1199
#> 3007 0.8330175   1199
#> 3008 0.8330175   1199
#> 3009 0.8330175   1199
#> 3010 0.8330175   1199
#> 3011 0.8330175   1199
#> 3012 0.8330175   1199
#> 3013 0.8330175   1199
#> 3014 0.8330175   1199
#> 3015 0.8330175   1199
#> 3016 0.8330175   1199
#> 3017 0.8330175   1199
#> 3018 0.8330175   1199
#> 3019 0.8330175   1199
#> 3020 0.8330175   1199
#> 3021 0.8330175   1199
#> 3022 0.8330175   1199
#> 3023 0.8330175   1199
#> 3024 0.8330175   1199
#> 3025 0.8330175   1199
#> 3026 0.8330175   1199
#> 3027 0.8330175   1199
#> 3028 0.8330175   1199
#> 3029 0.8330175   1199
#> 3030 0.8330175   1199
#> 3031 0.8330175   1199
#> 3032 0.8330175   1199
#> 3033 0.8330175   1199
#> 3034 0.8330175   1199
#> 3035 0.8330175   1199
#> 3036 0.8330175   1199
#> 3037 0.8330175   1199
#> 3038 0.8330175   1199
#> 3039 0.8330175   1199
#> 3040 0.8330175   1199
#> 3041 0.8330175   1199
#> 3042 0.8330175   1199
#> 3043 0.8330175   1199
#> 3044 0.8330175   1199
#> 3045 0.8330175   1199
#> 3046 0.8330175   1199
#> 3047 0.8330175   1199
#> 3048 0.8330175   1199
#> 3049 0.8330175   1199
#> 3050 0.8330175   1199
#> 3051 0.8330175   1199
#> 3052 0.8330175   1199
#> 3053 0.8330175   1199
#> 3054 0.8330175   1199
#> 3055 0.8330175   1199
#> 3056 0.8330175   1199
#> 3057 0.8330175   1199
#> 3058 0.8330175   1199
#> 3059 0.8330175   1199
#> 3060 0.8330175   1199
#> 3061 0.8330175   1199
#> 3062 0.8330175   1199
#> 3063 0.8330175   1199
#> 3064 0.8330175   1199
#> 3065 0.8330175   1199
#> 3066 0.8330175   1199
#> 3067 0.8330175   1199
#> 3068 0.8330175   1199
#> 3069 0.8330175   1199
#> 3070 0.8330175   1199
#> 3071 0.8330175   1199
#> 3072 0.8330175   1199
#> 3073 0.8330175   1199
#> 3074 0.8330175   1199
#> 3075 0.8330175   1199
#> 3076 0.8330175   1199
#> 3077 0.8330175   1199
#> 3078 0.8330175   1199
#> 3079 0.8330175   1199
#> 3080 0.8330175   1199
#> 3081 0.8330175   1199
#> 3082 0.8330175   1199
#> 3083 0.8330175   1199
#> 3084 0.8330175   1199
#> 3085 0.8330175   1199
#> 3086 0.8330175   1199
#> 3087 0.8330175   1199
#> 3088 0.8330175   1199
#> 3089 0.8330175   1199
#> 3090 0.8330175   1199
#> 3091 0.8330175   1199
#> 3092 0.8330175   1199
#> 3093 0.8330175   1199
#> 3094 0.8330175   1199
#> 3095 0.8330175   1199
#> 3096 0.8330175   1199
#> 3097 0.8330175   1199
#> 3098 0.8330175   1199
#> 3099 0.8330175   1199
#> 3100 0.8330175   1199
#> 3101 0.8330175   1199
#> 3102 0.8330175   1199
#> 3103 0.8330175   1199
#> 3104 0.8330175   1199
#> 3105 0.8330175   1199
#> 3106 0.8330175   1199
#> 3107 0.8330175   1199
#> 3108 0.8330175   1199
#> 3109 0.8330175   1199
#> 3110 0.8330175   1199
#> 3111 0.8330175   1199
#> 3112 0.8330175   1199
#> 3113 0.8330175   1199
#> 3114 0.8330175   1199
#> 3115 0.8330175   1199
#> 3116 0.8330175   1199
#> 3117 0.8330175   1199
#> 3118 0.8330175   1199
#> 3119 0.8330175   1199
#> 3120 0.8330175   1199
#> 3121 0.8330175   1199
#> 3122 0.8330175   1199
#> 3123 0.8330175   1199
#> 3124 0.8330175   1199
#> 3125 0.8330175   1199
#> 3126 0.8330175   1199
#> 3127 0.8330175   1199
#> 3128 0.8330175   1199
#> 3129 0.8330175   1199
#> 3130 0.8330175   1199
#> 3131 0.8330175   1199
#> 3132 0.8330175   1199
#> 3133 0.8330175   1199
#> 3134 0.8330175   1199
#> 3135 0.8330175   1199
#> 3136 0.8330175   1199
#> 3137 0.8330175   1199
#> 3138 0.8330175   1199
#> 3139 0.8330175   1199
#> 3140 0.8330175   1199
#> 3141 0.8330175   1199
#> 3142 0.8330175   1199
#> 3143 0.8330175   1199
#> 3144 0.8330175   1199
#> 3145 0.8330175   1199
#> 3146 0.8330175   1199
#> 3147 0.8330175   1199
#> 3148 0.8330175   1199
#> 3149 0.8330175   1199
#> 3150 0.8330175   1199
#> 3151 0.8330175   1199
#> 3152 0.8330175   1199
#> 3153 0.8330175   1199
#> 3154 0.8330175   1199
#> 3155 0.8330175   1199
#> 3156 0.8330175   1199
#> 3157 0.8330175   1199
#> 3158 0.8330175   1199
#> 3159 0.8330175   1199
#> 3160 0.8330175   1199
#> 3161 0.8330175   1199
#> 3162 0.8330175   1199
#> 3163 0.8330175   1199
#> 3164 0.8330175   1199
#> 3165 0.8330175   1199
#> 3166 0.8330175   1199
#> 3167 0.8330175   1199
#> 3168 0.8330175   1199
#> 3169 0.8330175   1199
#> 3170 0.8330175   1199
#> 3171 0.8330175   1199
#> 3172 0.8330175   1199
#> 3173 0.8330175   1199
#> 3174 0.8330175   1199
#> 3175 0.8330175   1199
#> 3176 0.8330175   1199
#> 3177 0.8330175   1199
#> 3178 0.8330175   1199
#> 3179 0.8330175   1199
#> 3180 0.8330175   1199
#> 3181 0.8330175   1199
#> 3182 0.8330175   1199
#> 3183 0.8330175   1199
#> 3184 0.8330175   1199
#> 3185 0.8330175   1199
#> 3186 0.8330175   1199
#> 3187 0.8330175   1199
#> 3188 0.8330175   1199
#> 3189 0.8330175   1199
#> 3190 0.8330175   1199
#> 3191 0.8330175   1199
#> 3192 0.8330175   1199
#> 3193 0.8330175   1199
#> 3194 0.8330175   1199
#> 3195 0.8330175   1199
#> 3196 0.8330175   1199
#> 3197 0.8330175   1199
#> 3198 0.8330175   1199
#> 3199 0.8330175   1199
#> 3200 0.8330175   1199
#> 3201 0.8330175   1199
#> 3202 0.8330175   1199
#> 3203 0.8330175   1199
#> 3204 0.8330175   1199
#> 3205 0.8330175   1199
#> 3206 0.8330175   1199
#> 3207 0.8330175   1199
#> 3208 0.8330175   1199
#> 3209 0.8330175   1199
#> 3210 0.8330175   1199
#> 3211 0.8330175   1199
#> 3212 0.8330175   1199
#> 3213 0.8330175   1199
#> 3214 0.8330175   1199
#> 3215 0.8330175   1199
#> 3216 0.8330175   1199
#> 3217 0.8330175   1199
#> 3218 0.8330175   1199
#> 3219 0.8330175   1199
#> 3220 0.8330175   1199
#> 3221 0.8330175   1199
#> 3222 0.8330175   1199
#> 3223 0.8330175   1199
#> 3224 0.8330175   1199
#> 3225 0.8330175   1199
#> 3226 0.8330175   1199
#> 3227 0.8330175   1199
#> 3228 0.8330175   1199
#> 3229 0.8330175   1199
#> 3230 0.8330175   1199
#> 3231 0.8330175   1199
#> 3232 0.8330175   1199
#> 3233 0.8330175   1199
#> 3234 0.8330175   1199
#> 3235 0.8330175   1199
#> 3236 0.8330175   1199
#> 3237 0.8330175   1199
#> 3238 0.8330175   1199
#> 3239 0.8330175   1199
#> 3240 0.8330175   1199
#> 3241 0.8330175   1199
#> 3242 0.8330175   1199
#> 3243 0.8330175   1199
#> 3244 0.8330175   1199
#> 3245 0.8330175   1199
#> 3246 0.8330175   1199
#> 3247 0.8330175   1199
#> 3248 0.8330175   1199
#> 3249 0.8330175   1199
#> 3250 0.8330175   1199
#> 3251 0.8330175   1199
#> 3252 0.8330175   1199
#> 3253 0.8330175   1199
#> 3254 0.8330175   1199
#> 3255 0.8330175   1199
#> 3256 0.8330175   1199
#> 3257 0.8330175   1199
#> 3258 0.8330175   1199
#> 3259 0.8330175   1199
#> 3260 0.8330175   1199
#> 3261 0.8330175   1199
#> 3262 0.8330175   1199
#> 3263 0.8330175   1199
#> 3264 0.8330175   1199
#> 3265 0.8330175   1199
#> 3266 0.8330175   1199
#> 3267 0.8330175   1199
#> 3268 0.8330175   1199
#> 3269 0.8330175   1199
#> 3270 0.8330175   1199
#> 3271 0.8330175   1199
#> 3272 0.8330175   1199
#> 3273 0.8330175   1199
#> 3274 0.8330175   1199
#> 3275 0.8330175   1199
#> 3276 0.8330175   1199
#> 3277 0.8330175   1199
#> 3278 0.8330175   1199
#> 3279 0.8330175   1199
#> 3280 0.8330175   1199
#> 3281 0.8330175   1199
#> 3282 0.8330175   1199
#> 3283 0.8330175   1199
#> 3284 0.8330175   1199
#> 3285 0.8330175   1199
#> 3286 0.8330175   1199
#> 3287 0.8330175   1199
#> 3288 0.8330175   1199
#> 3289 0.8330175   1199
#> 3290 0.8330175   1199
#> 3291 0.8330175   1199
#> 3292 0.8330175   1199
#> 3293 0.8330175   1199
#> 3294 0.8330175   1199
#> 3295 0.8330175   1199
#> 3296 0.8330175   1199
#> 3297 0.8330175   1199
#> 3298 0.8330175   1199
#> 3299 0.8330175   1199
#> 3300 0.8330175   1199
#> 3301 0.8330175   1199
#> 3302 0.8330175   1199
#> 3303 0.8330175   1199
#> 3304 0.8330175   1199
#> 3305 0.8330175   1199
#> 3306 0.8330175   1199
#> 3307 0.8330175   1199
#> 3308 0.8330175   1199
#> 3309 0.8330175   1199
#> 3310 0.8330175   1199
#> 3311 0.8330175   1199
#> 3312 0.8330175   1199
#> 3313 0.8330175   1199
#> 3314 0.8330175   1199
#> 3315 0.8330175   1199
#> 3316 0.8330175   1199
#> 3317 0.8330175   1199
#> 3318 0.8330175   1199
#> 3319 0.8330175   1199
#> 3320 0.8330175   1199
#> 3321 0.8330175   1199
#> 3322 0.8330175   1199
#> 3323 0.8330175   1199
#> 3324 0.8330175   1199
#> 3325 0.8330175   1199
#> 3326 0.8330175   1199
#> 3327 0.8330175   1199
#> 3328 0.8330175   1199
#> 3329 0.8330175   1199
#> 3330 0.8330175   1199
#> 3331 0.8330175   1199
#> 3332 0.8330175   1199
#> 3333 0.8330175   1199
#> 3334 0.8330175   1199
#> 3335 0.8330175   1199
#> 3336 0.8330175   1199
#> 3337 0.8330175   1199
#> 3338 0.8330175   1199
#> 3339 0.8330175   1199
#> 3340 0.8330175   1199
#> 3341 0.8330175   1199
#> 3342 0.8330175   1199
#> 3343 0.8330175   1199
#> 3344 0.8330175   1199
#> 3345 0.8330175   1199
#> 3346 0.8330175   1199
#> 3347 0.8330175   1199
#> 3348 0.8330175   1199
#> 3349 0.8330175   1199
#> 3350 0.8330175   1199
#> 3351 0.8330175   1199
#> 3352 0.8330175   1199
#> 3353 0.8330175   1199
#> 3354 0.8330175   1199
#> 3355 0.8330175   1199
#> 3356 0.8330175   1199
#> 3357 0.8330175   1199
#> 3358 0.8330175   1199
#> 3359 0.8330175   1199
#> 3360 0.8330175   1199
#> 3361 0.8330175   1199
#> 3362 0.8330175   1199
#> 3363 0.8330175   1199
#> 3364 0.8330175   1199
#> 3365 0.8330175   1199
#> 3366 0.8330175   1199
#> 3367 0.8330175   1199
#> 3368 0.8330175   1199
#> 3369 0.8330175   1199
#> 3370 0.8330175   1199
#> 3371 0.8330175   1199
#> 3372 0.8330175   1199
#> 3373 0.8330175   1199
#> 3374 0.8330175   1199
#> 3375 0.8330175   1199
#> 3376 0.8330175   1199
#> 3377 0.8330175   1199
#> 3378 0.8330175   1199
#> 3379 0.8330175   1199
#> 3380 0.8330175   1199
#> 3381 0.8330175   1199
#> 3382 0.8330175   1199
#> 3383 0.8330175   1199
#> 3384 0.8330175   1199
#> 3385 0.8330175   1199
#> 3386 0.8330175   1199
#> 3387 0.8330175   1199
#> 3388 0.8330175   1199
#> 3389 0.8330175   1199
#> 3390 0.8330175   1199
#> 3391 0.8330175   1199
#> 3392 0.8330175   1199
#> 3393 0.8330175   1199
#> 3394 0.8330175   1199
#> 3395 0.8330175   1199
#> 3396 0.8330175   1199
#> 3397 0.8330175   1199
#> 3398 0.8330175   1199
#> 3399 0.8330175   1199
#> 3400 0.8330175   1199
#> 3401 0.8330175   1199
#> 3402 0.8330175   1199
#> 3403 0.8330175   1199
#> 3404 0.8330175   1199
#> 3405 0.8330175   1199
#> 3406 0.8330175   1199
#> 3407 0.8330175   1199
#> 3408 0.8330175   1199
#> 3409 0.8330175   1199
#> 3410 0.8330175   1199
#> 3411 0.8330175   1199
#> 3412 0.8330175   1199
#> 3413 0.8330175   1199
#> 3414 0.8330175   1199
#> 3415 0.8330175   1199
#> 3416 0.8330175   1199
#> 3417 0.8330175   1199
#> 3418 0.8330175   1199
#> 3419 0.8330175   1199
#> 3420 0.8330175   1199
#> 3421 0.8330175   1199
#> 3422 0.8330175   1199
#> 3423 0.8330175   1199
#> 3424 0.8330175   1199
#> 3425 0.8330175   1199
#> 3426 0.8330175   1199
#> 3427 0.8330175   1199
#> 3428 0.8330175   1199
#> 3429 0.8330175   1199
#> 3430 0.8330175   1199
#> 3431 0.8330175   1199
#> 3432 0.8330175   1199
#> 3433 0.8330175   1199
#> 3434 0.8330175   1199
#> 3435 0.8330175   1199
#> 3436 0.8330175   1199
#> 3437 0.8330175   1199
#> 3438 0.8330175   1199
#> 3439 0.8330175   1199
#> 3440 0.8330175   1199
#> 3441 0.8330175   1199
#> 3442 0.8330175   1199
#> 3443 0.8330175   1199
#> 3444 0.8330175   1199
#> 3445 0.8330175   1199
#> 3446 0.8330175   1199
#> 3447 0.8330175   1199
#> 3448 0.8330175   1199
#> 3449 0.8330175   1199
#> 3450 0.8330175   1199
#> 3451 0.8330175   1199
#> 3452 0.8330175   1199
#> 3453 0.8330175   1199
#> 3454 0.8330175   1199
#> 3455 0.8330175   1199
#> 3456 0.8330175   1199
#> 3457 0.8330175   1199
#> 3458 0.8330175   1199
#> 3459 0.8330175   1199
#> 3460 0.8330175   1199
#> 3461 0.8330175   1199
#> 3462 0.8330175   1199
#> 3463 0.8330175   1199
#> 3464 0.8330175   1199
#> 3465 0.8330175   1199
#> 3466 0.8330175   1199
#> 3467 0.8330175   1199
#> 3468 0.8330175   1199
#> 3469 0.8330175   1199
#> 3470 0.8330175   1199
#> 3471 0.8330175   1199
#> 3472 0.8330175   1199
#> 3473 0.8330175   1199
#> 3474 0.8330175   1199
#> 3475 0.8330175   1199
#> 3476 0.8330175   1199
#> 3477 0.8330175   1199
#> 3478 0.8330175   1199
#> 3479 0.8330175   1199
#> 3480 0.8330175   1199
#> 3481 0.8330175   1199
#> 3482 0.8330175   1199
#> 3483 0.8330175   1199
#> 3484 0.8330175   1199
#> 3485 0.8330175   1199
#> 3486 0.8330175   1199
#> 3487 0.8330175   1199
#> 3488 0.8330175   1199
#> 3489 0.8330175   1199
#> 3490 0.8330175   1199
#> 3491 0.8330175   1199
#> 3492 0.8330175   1199
#> 3493 0.8330175   1199
#> 3494 0.8330175   1199
#> 3495 0.8330175   1199
#> 3496 0.8330175   1199
#> 3497 0.8330175   1199
#> 3498 0.8330175   1199
#> 3499 0.8330175   1199
#> 3500 0.8330175   1199
#> 3501 0.8330175   1199
#> 3502 0.8330175   1199
#> 3503 0.8330175   1199
#> 3504 0.8330175   1199
#> 3505 0.8330175   1199
#> 3506 0.8330175   1199
#> 3507 0.8330175   1199
#> 3508 0.8330175   1199
#> 3509 0.8330175   1199
#> 3510 0.8330175   1199
#> 3511 0.8330175   1199
#> 3512 0.8330175   1199
#> 3513 0.8330175   1199
#> 3514 0.8330175   1199
#> 3515 0.8330175   1199
#> 3516 0.8330175   1199
#> 3517 0.8330175   1199
#> 3518 0.8330175   1199
#> 3519 0.8330175   1199
#> 3520 0.8330175   1199
#> 3521 0.8330175   1199
#> 3522 0.8330175   1199
#> 3523 0.8330175   1199
#> 3524 0.8330175   1199
#> 3525 0.8330175   1199
#> 3526 0.8330175   1199
#> 3527 0.8330175   1199
#> 3528 0.8330175   1199
#> 3529 0.8330175   1199
#> 3530 0.8330175   1199
#> 3531 0.8330175   1199
#> 3532 0.8330175   1199
#> 3533 0.8330175   1199
#> 3534 0.8330175   1199
#> 3535 0.8330175   1199
#> 3536 0.8330175   1199
#> 3537 0.8330175   1199
#> 3538 0.8330175   1199
#> 3539 0.8330175   1199
#> 3540 0.8330175   1199
#> 3541 0.8330175   1199
#> 3542 0.8330175   1199
#> 3543 0.8330175   1199
#> 3544 0.8330175   1199
#> 3545 0.8330175   1199
#> 3546 0.8330175   1199
#> 3547 0.8330175   1199
#> 3548 0.8330175   1199
#> 3549 0.8330175   1199
#> 3550 0.8330175   1199
#> 3551 0.8330175   1199
#> 3552 0.8330175   1199
#> 3553 0.8330175   1199
#> 3554 0.8330175   1199
#> 3555 0.8330175   1199
#> 3556 0.8330175   1199
#> 3557 0.8330175   1199
#> 3558 0.8330175   1199
#> 3559 0.8330175   1199
#> 3560 0.8330175   1199
#> 3561 0.8330175   1199
#> 3562 0.8330175   1199
#> 3563 0.8330175   1199
#> 3564 0.8330175   1199
#> 3565 0.8330175   1199
#> 3566 0.8330175   1199
#> 3567 0.8330175   1199
#> 3568 0.8330175   1199
#> 3569 0.8330175   1199
#> 3570 0.8330175   1199
#> 3571 0.8330175   1199
#> 3572 0.8330175   1199
#> 3573 0.8330175   1199
#> 3574 0.8330175   1199
#> 3575 0.8330175   1199
#> 3576 0.8330175   1199
```

[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
prepends `replicate`, `timepoint`, and `analysis` columns to the result
data. The `events` column shows how many TTE events were observed — with
n=1200 and `lambda0 = 0.08`, this should be large. `logrank_p` tests the
null of equal survival curves; `hr_cox` is the estimated hazard ratio
from the Cox model (values \< 1 favour treatment), and
`hr_ci_lo`/`hr_ci_hi` are its 95% CI bounds. Because the continuous
endpoint Y and TTE share the same latent drivers, replicates with a
small `p_ttest_y` tend to also have a small `logrank_p`, though the
correspondence is imperfect due to the different statistical tests and
the censoring mechanism.

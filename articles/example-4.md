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

## Scenario

We consider a **two-arm fixed design** with a **continuous endpoint
(Y)** and a **time-to-event (TTE)** endpoint per subject.  
The endpoints are **correlated** via a shared Gaussian latent structure.
Enrollment is piece-wise linear; analysis happens **once at the final
time**.

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
mu_trt <- 0.3      # treatment mean shift
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
  allocation = list(allocation),
  rho = rho,
  hr_trt_true = HR_trt
)
```

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

## Arms (Populations)

Per arm we simulate subject-level latent bivariate normals `(Z1, Z2)`
with correlation `rho`.  
- The **continuous** endpoint is `Y = mu_arm + sd_y * Z1`. - The **TTE**
endpoint uses a proportional-hazards exponential model with
subject-specific frailty:  
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
  trt = mk_population_generator(mu_y = mu_trt, lambda_arm = lambda_trt)
)
```

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
replicate_results <- do.call(rbind, lapply(seq_along(trials), function(i) {
  data.frame(
    replicate = i,
    trials[[i]]$results[[1]]$final,
    stringsAsFactors = FALSE
  )
}))

replicate_results
#>   replicate sample_size allocation rho hr_trt_true n_total n_pbo n_trt
#> 1         1        1200       1, 1 0.5         0.7    1200   600   600
#> 2         2        1200       1, 1 0.5         0.7    1200   600   600
#> 3         3        1200       1, 1 0.5         0.7    1200   600   600
#>      p_ttest_y    logrank_p    hr_cox  hr_ci_lo  hr_ci_hi events
#> 1 1.437333e-07 1.352989e-05 0.7768063 0.6930741 0.8706544   1198
#> 2 9.000125e-08 3.702888e-10 0.6950138 0.6198860 0.7792467   1200
#> 3 9.572157e-06 3.198372e-07 0.7430600 0.6628171 0.8330175   1199
```

# Example 7: Two-arm \| Fixed design \| Single continuous endpoint \| Bayesian Go/No-Go with placebo borrowing

``` r
# Core simulation framework
library(rxsim)

# Bayesian borrowing utilities
library(RBesT)
#> This is RBesT version 1.8.2 (released 2025-04-25, git-sha b9dab00)

# Helpers
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

set.seed(777)
```

## Scenario

A **two-arm, fixed design** trial (placebo vs treatment) with a **single
continuous endpoint** and a **single final analysis**.

- Placebo arm uses **historical borrowing** through an informative prior
  derived from historical placebo data and robustified with a
  non-informative component.
- Treatment arm uses a **non-informative prior**.
- Go/No-Go rule at final analysis:

$${\text{Go if}\mspace{6mu}}P\left( \Delta > \delta \mid \text{data} \right) \geq \gamma,$$

where $\Delta = \mu_{T} - \mu_{P}$, with $\delta = 0.1$ and
$\gamma = 0.8$.

``` r
# Trial design
sample_size <- 80
allocation  <- c(1, 1)
arms        <- c("placebo", "treatment")

# Data-generating truth for this simulation example
mu_placebo_true <- 0.00
mu_treat_true   <- 0.30
sigma_known     <- 1.0

# Enrollment/dropout profile (fixed design timeline)
enrollment <- list(
  end_time = c(4, 8),
  rate     = c(6, 10)
)
dropout <- list(
  end_time = c(4, 8),
  rate     = c(0, 0)
)

# Bayesian decision thresholds
delta <- 0.1
gamma <- 0.8

scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation = list(allocation),
  delta = delta,
  gamma = gamma
)
```

## Historical placebo dataset and priors

We include a small historical placebo dataset as study-level means and
sample sizes.

``` r
hist_placebo <- data.frame(
  study = c("H1", "H2", "H3"),
  n     = c(24, 18, 20),
  mean  = c(-0.05, 0.02, 0.00)
)

hist_placebo
#>   study  n  mean
#> 1    H1 24 -0.05
#> 2    H2 18  0.02
#> 3    H3 20  0.00

# Pooled historical summary used to derive an informative placebo prior
n_hist_total <- sum(hist_placebo$n)
m_hist_pooled <- weighted.mean(hist_placebo$mean, hist_placebo$n)

# Non-informative prior (very small pseudo-sample size in mn parametrization)
# Used as treatment prior and as base for historical update.
prior_noninf <- mixnorm(
  c(1, 0, 1e-6),
  sigma = sigma_known,
  param = "mn"
)

# Informative placebo prior from historical placebo summary
prior_placebo_inf <- postmix(
  prior_noninf,
  n = n_hist_total,
  m = m_hist_pooled
)
#> Using default prior reference scale 1

# Robustified placebo prior to protect against prior-data conflict
prior_placebo <- robustify(
  prior_placebo_inf,
  weight = 0.2,
  mean   = 0,
  n      = 1,
  sigma  = sigma_known
)

# Treatment prior remains non-informative
prior_treatment <- prior_noninf

summary(prior_placebo)
#>        mean          sd        2.5%       50.0%       97.5% 
#> -0.01083871  0.46144620 -1.15034971 -0.01312827  1.15034975
```

## Populations

``` r
population_generators <- list(
  placebo = function(n) {
    data.frame(
      id = 1:n,
      y = rnorm(n, mean = mu_placebo_true, sd = sigma_known),
      readout_time = 1
    )
  },
  treatment = function(n) {
    data.frame(
      id = 1:n,
      y = rnorm(n, mean = mu_treat_true, sd = sigma_known),
      readout_time = 1
    )
  }
)
```

## Trial Parameters

``` r
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn <- function(n) rexp(n, rate = 0.01)
```

## Final analysis trigger (Bayesian Go/No-Go)

Only one trigger is added at final analysis.

``` r
analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer) {
      dat <- df |>
        dplyr::filter(!is.na(enroll_time)) |>
        dplyr::mutate(arm = factor(arm, levels = c("placebo", "treatment")))

      y_p <- dat$y[dat$arm == "placebo"]
      y_t <- dat$y[dat$arm == "treatment"]

      n_p <- length(y_p)
      n_t <- length(y_t)

      m_p <- mean(y_p)
      m_t <- mean(y_t)

      # Posterior for each arm mean
      post_placebo <- RBesT::postmix(prior_placebo, n = n_p, m = m_p)
      post_treat   <- RBesT::postmix(prior_treatment, n = n_t, m = m_t)

      # Posterior probability for treatment effect exceeding delta
      prob_delta <- RBesT::pmixdiff(post_treat, post_placebo, delta, lower.tail = FALSE)

      # Posterior summaries for Delta = mu_T - mu_P
      delta_ci <- RBesT::qmixdiff(post_treat, post_placebo, c(0.025, 0.5, 0.975))

      data.frame(
        scenario,
        n_placebo = n_p,
        n_treatment = n_t,
        mean_placebo = m_p,
        mean_treatment = m_t,
        post_prob_delta = prob_delta,
        delta_q025 = delta_ci[1],
        delta_q500 = delta_ci[2],
        delta_q975 = delta_ci[3],
        decision_go = as.integer(prob_delta >= gamma),
        stringsAsFactors = FALSE
      )
    }
  )
)
```

## Trial

``` r
trials <- replicate_trial(
  trial_name = "example_7_bayes_two_arm_fixed",
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
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: example_7_bayes_two_arm_fixed_1
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
#>     name: example_7_bayes_two_arm_fixed_2
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
#>     name: example_7_bayes_two_arm_fixed_3
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
```

## Results

``` r
replicate_results <- do.call(rbind, lapply(seq_along(trials), function(i) {
  data.frame(
    replicate = i,
    trials[[i]]$results[[1]]$final,
    stringsAsFactors = FALSE
  )
}))

replicate_results
#>   replicate sample_size allocation delta gamma n_placebo n_treatment
#> 1         1          80       1, 1   0.1   0.8        40          40
#> 2         2          80       1, 1   0.1   0.8        40          40
#> 3         3          80       1, 1   0.1   0.8        40          40
#>   mean_placebo mean_treatment post_prob_delta  delta_q025 delta_q500 delta_q975
#> 1    0.1581546      0.3864720       0.8826182 -0.05183437  0.3272105  0.6968002
#> 2    0.1669418      0.1654918       0.5048670 -0.27784403  0.1023216  0.4720979
#> 3   -0.4034483      0.3106430       0.9825414  0.12871903  0.5184149  0.9892527
#>   decision_go
#> 1           1
#> 2           0
#> 3           1
```

`decision_go = 1` indicates **Go**, and `decision_go = 0` indicates
**No-Go** under the criterion
$P\left( \Delta > 0.1 \mid \text{data} \right) \geq 0.8$.

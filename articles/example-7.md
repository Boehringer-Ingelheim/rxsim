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

Early development trials often face small sample sizes where borrowing
information from historical placebo data can increase efficiency and
reduce required sample size. This example demonstrates a Bayesian
Go/No-Go decision rule for a Phase IIa trial, where the placebo
posterior is informed by historical data via a robust mixture prior, and
the treatment posterior uses a non-informative prior. A “Go” decision is
issued when the posterior probability of a meaningful treatment effect
(Δ \> δ) exceeds a predefined threshold γ.

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

`delta = 0.1` defines the minimum treatment-placebo difference of
clinical relevance — the effect must exceed this threshold for the drug
to be worth advancing. `gamma = 0.8` requires 80% posterior probability
of exceeding that threshold before a Go decision is issued, balancing
false-positive risk against the cost of a missed opportunity. The true
simulated effect `mu_treat = 0.30` lies above `delta`, so a well-powered
trial should frequently result in a Go.

``` r
# Trial design
sample_size <- 80
allocation  <- c(1, 1)
arms        <- c("placebo", "treatment")

# Data-generating truth for this simulation example
mu_placebo_true <- 0.00
delta_true      <- 0.30  # true treatment - placebo mean difference
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
  allocation  = list(allocation),
  delta_true  = delta_true,
  delta       = delta,
  gamma       = gamma
)
```

## Historical placebo dataset and priors

We include a small historical placebo dataset as study-level means and
sample sizes.

Studies H1, H2, and H3 contribute a combined 62 historical placebo
observations and are summarised into an informative prior via
[`postmix()`](https://opensource.nibr.com/RBesT/reference/postmix.html)
— a Bayesian update of the non-informative starting prior with the
pooled historical data.
[`robustify()`](https://opensource.nibr.com/RBesT/reference/robustify.html)
blends this informative component with a vague component (20% weight),
producing a mixture prior that protects against prior-data conflict: if
the current trial’s placebo behaves unexpectedly, the vague component
limits how strongly the historical data pulls the posterior. The
treatment arm uses `prior_noninf` (effectively a flat prior), so the
treatment posterior is driven entirely by the current data.

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
      y = rnorm(n, mean = mu_placebo_true + delta_true, sd = sigma_known),
      readout_time = 1
    )
  }
)
```

## Trial Parameters

`enrollment_fn` and `dropout_fn` supply the inter-arrival and
time-to-dropout distributions for individual subjects. These are passed
to
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md),
which uses them when building each replicate’s enrollment schedule.

``` r
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn <- function(n) rexp(n, rate = 0.01)
```

## Final analysis trigger (Bayesian Go/No-Go)

Only one trigger is added at final analysis.

`RBesT::pmixdiff(post_treat, post_placebo, delta, lower.tail = FALSE)`
evaluates P(μ_T − μ_P \> δ \| data) by numerical integration over the
mixture posterior distributions for treatment and placebo.
[`qmixdiff()`](https://opensource.nibr.com/RBesT/reference/mixdiff.html)
computes quantiles of the same posterior difference distribution,
yielding a 95% credible interval for Δ = μ_T − μ_P.

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

[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
row-binds analysis outputs across all replicates and prepends
`replicate` (integer index), `timepoint` (calendar time at which the
analysis fired), and `analysis` (the analysis name) to each row.
`post_prob_delta` is the key decision metric: values at or above
`gamma = 0.8` yield `decision_go = 1` (Go) and values below yield
`decision_go = 0` (No-Go). The credible interval columns `delta_q025`,
`delta_q500`, and `delta_q975` characterise the posterior uncertainty
about the treatment effect Δ; with n=80 and σ=1, the posterior is still
relatively wide. In practice you would run thousands of replicates and
report the Go probability under both the null (Δ = 0) and the
alternative (e.g., Δ = 0.30) to assess the design’s operating
characteristics. See [Example
8](https://boehringer-ingelheim.github.io/rxsim/articles/example-8.md)
for a seamless design that builds on this Bayesian decision rule.

``` r
replicate_results <- collect_results(trials)
replicate_results
#>     replicate  timepoint analysis sample_size allocation delta_true delta gamma
#> 1           1   78.01039    final          80       1, 1        0.3   0.1   0.8
#> 2           1  166.62671    final          80       1, 1        0.3   0.1   0.8
#> 3           1  274.74269    final          80       1, 1        0.3   0.1   0.8
#> 4           1  341.00793    final          80       1, 1        0.3   0.1   0.8
#> 5           1  397.64554    final          80       1, 1        0.3   0.1   0.8
#> 6           1  542.07960    final          80       1, 1        0.3   0.1   0.8
#> 7           1  993.52447    final          80       1, 1        0.3   0.1   0.8
#> 8           1 1030.14181    final          80       1, 1        0.3   0.1   0.8
#> 9           1 1046.09099    final          80       1, 1        0.3   0.1   0.8
#> 10          1 1266.02529    final          80       1, 1        0.3   0.1   0.8
#> 11          1 1600.54126    final          80       1, 1        0.3   0.1   0.8
#> 12          1 1700.98321    final          80       1, 1        0.3   0.1   0.8
#> 13          1 1765.90494    final          80       1, 1        0.3   0.1   0.8
#> 14          1 1777.88160    final          80       1, 1        0.3   0.1   0.8
#> 15          1 1853.14483    final          80       1, 1        0.3   0.1   0.8
#> 16          1 1890.49139    final          80       1, 1        0.3   0.1   0.8
#> 17          1 1988.66721    final          80       1, 1        0.3   0.1   0.8
#> 18          1 2048.00362    final          80       1, 1        0.3   0.1   0.8
#> 19          1 2204.84397    final          80       1, 1        0.3   0.1   0.8
#> 20          1 2205.45194    final          80       1, 1        0.3   0.1   0.8
#> 21          1 2274.05430    final          80       1, 1        0.3   0.1   0.8
#> 22          1 2304.29670    final          80       1, 1        0.3   0.1   0.8
#> 23          1 2494.47980    final          80       1, 1        0.3   0.1   0.8
#> 24          1 2589.90916    final          80       1, 1        0.3   0.1   0.8
#> 25          1 2609.87077    final          80       1, 1        0.3   0.1   0.8
#> 26          1 2802.56772    final          80       1, 1        0.3   0.1   0.8
#> 27          1 2814.75858    final          80       1, 1        0.3   0.1   0.8
#> 28          1 2875.69398    final          80       1, 1        0.3   0.1   0.8
#> 29          1 3386.66459    final          80       1, 1        0.3   0.1   0.8
#> 30          1 3619.37171    final          80       1, 1        0.3   0.1   0.8
#> 31          1 3691.53725    final          80       1, 1        0.3   0.1   0.8
#> 32          1 4082.85376    final          80       1, 1        0.3   0.1   0.8
#> 33          1 4268.87048    final          80       1, 1        0.3   0.1   0.8
#> 34          1 4320.64179    final          80       1, 1        0.3   0.1   0.8
#> 35          1 4526.06241    final          80       1, 1        0.3   0.1   0.8
#> 36          1 4567.62111    final          80       1, 1        0.3   0.1   0.8
#> 37          1 4668.73265    final          80       1, 1        0.3   0.1   0.8
#> 38          1 4755.41828    final          80       1, 1        0.3   0.1   0.8
#> 39          1 4811.44714    final          80       1, 1        0.3   0.1   0.8
#> 40          1 4848.88559    final          80       1, 1        0.3   0.1   0.8
#> 41          1 4948.88540    final          80       1, 1        0.3   0.1   0.8
#> 42          1 5227.15583    final          80       1, 1        0.3   0.1   0.8
#> 43          1 5257.56630    final          80       1, 1        0.3   0.1   0.8
#> 44          1 5355.23288    final          80       1, 1        0.3   0.1   0.8
#> 45          1 5416.18547    final          80       1, 1        0.3   0.1   0.8
#> 46          1 5584.38678    final          80       1, 1        0.3   0.1   0.8
#> 47          1 5620.74779    final          80       1, 1        0.3   0.1   0.8
#> 48          1 5835.61550    final          80       1, 1        0.3   0.1   0.8
#> 49          1 5888.72060    final          80       1, 1        0.3   0.1   0.8
#> 50          1 5898.19732    final          80       1, 1        0.3   0.1   0.8
#> 51          1 5915.49912    final          80       1, 1        0.3   0.1   0.8
#> 52          1 5966.22425    final          80       1, 1        0.3   0.1   0.8
#> 53          1 6001.27279    final          80       1, 1        0.3   0.1   0.8
#> 54          1 6216.09594    final          80       1, 1        0.3   0.1   0.8
#> 55          1 6461.56539    final          80       1, 1        0.3   0.1   0.8
#> 56          1 6604.41916    final          80       1, 1        0.3   0.1   0.8
#> 57          1 6625.90225    final          80       1, 1        0.3   0.1   0.8
#> 58          1 6626.04747    final          80       1, 1        0.3   0.1   0.8
#> 59          1 6912.07968    final          80       1, 1        0.3   0.1   0.8
#> 60          1 7008.96686    final          80       1, 1        0.3   0.1   0.8
#> 61          1 7076.90002    final          80       1, 1        0.3   0.1   0.8
#> 62          1 7187.95434    final          80       1, 1        0.3   0.1   0.8
#> 63          1 7275.23947    final          80       1, 1        0.3   0.1   0.8
#> 64          1 7276.16337    final          80       1, 1        0.3   0.1   0.8
#> 65          1 7312.20614    final          80       1, 1        0.3   0.1   0.8
#> 66          1 7338.43520    final          80       1, 1        0.3   0.1   0.8
#> 67          1 7339.54863    final          80       1, 1        0.3   0.1   0.8
#> 68          1 7504.97403    final          80       1, 1        0.3   0.1   0.8
#> 69          1 7569.27389    final          80       1, 1        0.3   0.1   0.8
#> 70          1 7708.18198    final          80       1, 1        0.3   0.1   0.8
#> 71          1 7776.55727    final          80       1, 1        0.3   0.1   0.8
#> 72          1 7862.20170    final          80       1, 1        0.3   0.1   0.8
#> 73          1 7909.30003    final          80       1, 1        0.3   0.1   0.8
#> 74          1 7959.38094    final          80       1, 1        0.3   0.1   0.8
#> 75          1 8061.10855    final          80       1, 1        0.3   0.1   0.8
#> 76          1 8100.11077    final          80       1, 1        0.3   0.1   0.8
#> 77          1 8198.15124    final          80       1, 1        0.3   0.1   0.8
#> 78          1 8257.56559    final          80       1, 1        0.3   0.1   0.8
#> 79          1 8379.88677    final          80       1, 1        0.3   0.1   0.8
#> 80          1 8420.19824    final          80       1, 1        0.3   0.1   0.8
#> 81          2   72.27821    final          80       1, 1        0.3   0.1   0.8
#> 82          2  232.93573    final          80       1, 1        0.3   0.1   0.8
#> 83          2  268.97090    final          80       1, 1        0.3   0.1   0.8
#> 84          2  303.69245    final          80       1, 1        0.3   0.1   0.8
#> 85          2  370.99856    final          80       1, 1        0.3   0.1   0.8
#> 86          2  425.48186    final          80       1, 1        0.3   0.1   0.8
#> 87          2  469.94651    final          80       1, 1        0.3   0.1   0.8
#> 88          2  523.53947    final          80       1, 1        0.3   0.1   0.8
#> 89          2  551.41397    final          80       1, 1        0.3   0.1   0.8
#> 90          2  760.20443    final          80       1, 1        0.3   0.1   0.8
#> 91          2  765.20555    final          80       1, 1        0.3   0.1   0.8
#> 92          2  876.66320    final          80       1, 1        0.3   0.1   0.8
#> 93          2  923.23818    final          80       1, 1        0.3   0.1   0.8
#> 94          2  957.10604    final          80       1, 1        0.3   0.1   0.8
#> 95          2 1110.50539    final          80       1, 1        0.3   0.1   0.8
#> 96          2 1412.34485    final          80       1, 1        0.3   0.1   0.8
#> 97          2 1611.75866    final          80       1, 1        0.3   0.1   0.8
#> 98          2 1642.23767    final          80       1, 1        0.3   0.1   0.8
#> 99          2 1756.04808    final          80       1, 1        0.3   0.1   0.8
#> 100         2 1779.66003    final          80       1, 1        0.3   0.1   0.8
#> 101         2 1820.23743    final          80       1, 1        0.3   0.1   0.8
#> 102         2 1841.89988    final          80       1, 1        0.3   0.1   0.8
#> 103         2 1886.16137    final          80       1, 1        0.3   0.1   0.8
#> 104         2 1887.57464    final          80       1, 1        0.3   0.1   0.8
#> 105         2 1961.61067    final          80       1, 1        0.3   0.1   0.8
#> 106         2 2164.22037    final          80       1, 1        0.3   0.1   0.8
#> 107         2 2175.48914    final          80       1, 1        0.3   0.1   0.8
#> 108         2 2457.28059    final          80       1, 1        0.3   0.1   0.8
#> 109         2 2458.24513    final          80       1, 1        0.3   0.1   0.8
#> 110         2 2461.51360    final          80       1, 1        0.3   0.1   0.8
#> 111         2 2485.92905    final          80       1, 1        0.3   0.1   0.8
#> 112         2 2525.24309    final          80       1, 1        0.3   0.1   0.8
#> 113         2 2622.41460    final          80       1, 1        0.3   0.1   0.8
#> 114         2 2830.12043    final          80       1, 1        0.3   0.1   0.8
#> 115         2 2833.16944    final          80       1, 1        0.3   0.1   0.8
#> 116         2 2906.41171    final          80       1, 1        0.3   0.1   0.8
#> 117         2 2926.71198    final          80       1, 1        0.3   0.1   0.8
#> 118         2 3053.92052    final          80       1, 1        0.3   0.1   0.8
#> 119         2 3153.91626    final          80       1, 1        0.3   0.1   0.8
#> 120         2 3174.25744    final          80       1, 1        0.3   0.1   0.8
#> 121         2 3264.25914    final          80       1, 1        0.3   0.1   0.8
#> 122         2 3423.66688    final          80       1, 1        0.3   0.1   0.8
#> 123         2 3685.21561    final          80       1, 1        0.3   0.1   0.8
#> 124         2 3826.62390    final          80       1, 1        0.3   0.1   0.8
#> 125         2 3835.03646    final          80       1, 1        0.3   0.1   0.8
#> 126         2 4064.78948    final          80       1, 1        0.3   0.1   0.8
#> 127         2 4076.08980    final          80       1, 1        0.3   0.1   0.8
#> 128         2 4092.68989    final          80       1, 1        0.3   0.1   0.8
#> 129         2 4132.61804    final          80       1, 1        0.3   0.1   0.8
#> 130         2 4198.34451    final          80       1, 1        0.3   0.1   0.8
#> 131         2 4246.89554    final          80       1, 1        0.3   0.1   0.8
#> 132         2 4295.48743    final          80       1, 1        0.3   0.1   0.8
#> 133         2 4368.81326    final          80       1, 1        0.3   0.1   0.8
#> 134         2 4462.32322    final          80       1, 1        0.3   0.1   0.8
#> 135         2 4652.33416    final          80       1, 1        0.3   0.1   0.8
#> 136         2 4655.58399    final          80       1, 1        0.3   0.1   0.8
#> 137         2 4962.54541    final          80       1, 1        0.3   0.1   0.8
#> 138         2 5019.79414    final          80       1, 1        0.3   0.1   0.8
#> 139         2 5116.77034    final          80       1, 1        0.3   0.1   0.8
#> 140         2 5209.60148    final          80       1, 1        0.3   0.1   0.8
#> 141         2 5275.25817    final          80       1, 1        0.3   0.1   0.8
#> 142         2 5452.44413    final          80       1, 1        0.3   0.1   0.8
#> 143         2 5538.56539    final          80       1, 1        0.3   0.1   0.8
#> 144         2 5636.65654    final          80       1, 1        0.3   0.1   0.8
#> 145         2 5717.70576    final          80       1, 1        0.3   0.1   0.8
#> 146         2 5769.81900    final          80       1, 1        0.3   0.1   0.8
#> 147         2 5775.07185    final          80       1, 1        0.3   0.1   0.8
#> 148         2 5825.52155    final          80       1, 1        0.3   0.1   0.8
#> 149         2 5832.12600    final          80       1, 1        0.3   0.1   0.8
#> 150         2 5835.52635    final          80       1, 1        0.3   0.1   0.8
#> 151         2 5865.30222    final          80       1, 1        0.3   0.1   0.8
#> 152         2 5895.14827    final          80       1, 1        0.3   0.1   0.8
#> 153         2 6320.35624    final          80       1, 1        0.3   0.1   0.8
#> 154         2 6334.37112    final          80       1, 1        0.3   0.1   0.8
#> 155         2 6478.15678    final          80       1, 1        0.3   0.1   0.8
#> 156         2 6509.34144    final          80       1, 1        0.3   0.1   0.8
#> 157         2 6511.41412    final          80       1, 1        0.3   0.1   0.8
#> 158         2 6553.27902    final          80       1, 1        0.3   0.1   0.8
#> 159         2 6584.02090    final          80       1, 1        0.3   0.1   0.8
#> 160         2 6833.71939    final          80       1, 1        0.3   0.1   0.8
#> 161         2 7185.92084    final          80       1, 1        0.3   0.1   0.8
#> 162         3   89.63501    final          80       1, 1        0.3   0.1   0.8
#> 163         3   92.44571    final          80       1, 1        0.3   0.1   0.8
#> 164         3  107.27336    final          80       1, 1        0.3   0.1   0.8
#> 165         3  155.03158    final          80       1, 1        0.3   0.1   0.8
#> 166         3  163.22004    final          80       1, 1        0.3   0.1   0.8
#> 167         3  211.50328    final          80       1, 1        0.3   0.1   0.8
#> 168         3  274.37196    final          80       1, 1        0.3   0.1   0.8
#> 169         3  350.40933    final          80       1, 1        0.3   0.1   0.8
#> 170         3  391.71272    final          80       1, 1        0.3   0.1   0.8
#> 171         3  439.60209    final          80       1, 1        0.3   0.1   0.8
#> 172         3  485.74316    final          80       1, 1        0.3   0.1   0.8
#> 173         3  526.20723    final          80       1, 1        0.3   0.1   0.8
#> 174         3  602.91463    final          80       1, 1        0.3   0.1   0.8
#> 175         3  653.13648    final          80       1, 1        0.3   0.1   0.8
#> 176         3  768.94793    final          80       1, 1        0.3   0.1   0.8
#> 177         3  786.14661    final          80       1, 1        0.3   0.1   0.8
#> 178         3  866.24212    final          80       1, 1        0.3   0.1   0.8
#> 179         3  872.28363    final          80       1, 1        0.3   0.1   0.8
#> 180         3  934.57354    final          80       1, 1        0.3   0.1   0.8
#> 181         3  987.19940    final          80       1, 1        0.3   0.1   0.8
#> 182         3 1017.71064    final          80       1, 1        0.3   0.1   0.8
#> 183         3 1200.11344    final          80       1, 1        0.3   0.1   0.8
#> 184         3 1306.67993    final          80       1, 1        0.3   0.1   0.8
#> 185         3 1360.44531    final          80       1, 1        0.3   0.1   0.8
#> 186         3 1362.00797    final          80       1, 1        0.3   0.1   0.8
#> 187         3 1388.21354    final          80       1, 1        0.3   0.1   0.8
#> 188         3 1578.95235    final          80       1, 1        0.3   0.1   0.8
#> 189         3 1599.65901    final          80       1, 1        0.3   0.1   0.8
#> 190         3 1605.60287    final          80       1, 1        0.3   0.1   0.8
#> 191         3 1609.82020    final          80       1, 1        0.3   0.1   0.8
#> 192         3 1648.13453    final          80       1, 1        0.3   0.1   0.8
#> 193         3 1748.21076    final          80       1, 1        0.3   0.1   0.8
#> 194         3 1844.85217    final          80       1, 1        0.3   0.1   0.8
#> 195         3 1899.01464    final          80       1, 1        0.3   0.1   0.8
#> 196         3 1918.59842    final          80       1, 1        0.3   0.1   0.8
#> 197         3 2025.58275    final          80       1, 1        0.3   0.1   0.8
#> 198         3 2088.78760    final          80       1, 1        0.3   0.1   0.8
#> 199         3 2161.93444    final          80       1, 1        0.3   0.1   0.8
#> 200         3 2299.87567    final          80       1, 1        0.3   0.1   0.8
#> 201         3 2508.37171    final          80       1, 1        0.3   0.1   0.8
#> 202         3 2552.18321    final          80       1, 1        0.3   0.1   0.8
#> 203         3 2623.58744    final          80       1, 1        0.3   0.1   0.8
#> 204         3 2807.08992    final          80       1, 1        0.3   0.1   0.8
#> 205         3 2895.90606    final          80       1, 1        0.3   0.1   0.8
#> 206         3 3291.70563    final          80       1, 1        0.3   0.1   0.8
#> 207         3 3312.89467    final          80       1, 1        0.3   0.1   0.8
#> 208         3 3405.24745    final          80       1, 1        0.3   0.1   0.8
#> 209         3 3422.07365    final          80       1, 1        0.3   0.1   0.8
#> 210         3 3532.45172    final          80       1, 1        0.3   0.1   0.8
#> 211         3 3560.52150    final          80       1, 1        0.3   0.1   0.8
#> 212         3 3620.08800    final          80       1, 1        0.3   0.1   0.8
#> 213         3 3938.12837    final          80       1, 1        0.3   0.1   0.8
#> 214         3 3963.60226    final          80       1, 1        0.3   0.1   0.8
#> 215         3 3987.83905    final          80       1, 1        0.3   0.1   0.8
#> 216         3 4221.72960    final          80       1, 1        0.3   0.1   0.8
#> 217         3 4349.81402    final          80       1, 1        0.3   0.1   0.8
#> 218         3 4719.78920    final          80       1, 1        0.3   0.1   0.8
#> 219         3 4750.43111    final          80       1, 1        0.3   0.1   0.8
#> 220         3 4880.62711    final          80       1, 1        0.3   0.1   0.8
#> 221         3 4934.73439    final          80       1, 1        0.3   0.1   0.8
#> 222         3 4960.69961    final          80       1, 1        0.3   0.1   0.8
#> 223         3 4965.92877    final          80       1, 1        0.3   0.1   0.8
#> 224         3 4974.81858    final          80       1, 1        0.3   0.1   0.8
#> 225         3 5021.99314    final          80       1, 1        0.3   0.1   0.8
#> 226         3 5072.50231    final          80       1, 1        0.3   0.1   0.8
#> 227         3 5072.95339    final          80       1, 1        0.3   0.1   0.8
#> 228         3 5121.03919    final          80       1, 1        0.3   0.1   0.8
#> 229         3 5190.72293    final          80       1, 1        0.3   0.1   0.8
#> 230         3 5257.31143    final          80       1, 1        0.3   0.1   0.8
#> 231         3 5338.31098    final          80       1, 1        0.3   0.1   0.8
#> 232         3 5359.66163    final          80       1, 1        0.3   0.1   0.8
#> 233         3 5413.62218    final          80       1, 1        0.3   0.1   0.8
#> 234         3 5528.46797    final          80       1, 1        0.3   0.1   0.8
#> 235         3 5584.10656    final          80       1, 1        0.3   0.1   0.8
#> 236         3 5608.75695    final          80       1, 1        0.3   0.1   0.8
#> 237         3 5717.07576    final          80       1, 1        0.3   0.1   0.8
#> 238         3 6044.72272    final          80       1, 1        0.3   0.1   0.8
#> 239         3 6105.39837    final          80       1, 1        0.3   0.1   0.8
#> 240         3 6201.62981    final          80       1, 1        0.3   0.1   0.8
#> 241         3 6625.67574    final          80       1, 1        0.3   0.1   0.8
#> 242         3 6701.99306    final          80       1, 1        0.3   0.1   0.8
#>     n_placebo n_treatment mean_placebo mean_treatment post_prob_delta
#> 1          40          40    0.1581546      0.3864720       0.8826182
#> 2          40          40    0.1581546      0.3864720       0.8826182
#> 3          40          40    0.1581546      0.3864720       0.8826182
#> 4          40          40    0.1581546      0.3864720       0.8826182
#> 5          40          40    0.1581546      0.3864720       0.8826182
#> 6          40          40    0.1581546      0.3864720       0.8826182
#> 7          40          40    0.1581546      0.3864720       0.8826182
#> 8          40          40    0.1581546      0.3864720       0.8826182
#> 9          40          40    0.1581546      0.3864720       0.8826182
#> 10         40          40    0.1581546      0.3864720       0.8826182
#> 11         40          40    0.1581546      0.3864720       0.8826182
#> 12         40          40    0.1581546      0.3864720       0.8826182
#> 13         40          40    0.1581546      0.3864720       0.8826182
#> 14         40          40    0.1581546      0.3864720       0.8826182
#> 15         40          40    0.1581546      0.3864720       0.8826182
#> 16         40          40    0.1581546      0.3864720       0.8826182
#> 17         40          40    0.1581546      0.3864720       0.8826182
#> 18         40          40    0.1581546      0.3864720       0.8826182
#> 19         40          40    0.1581546      0.3864720       0.8826182
#> 20         40          40    0.1581546      0.3864720       0.8826182
#> 21         40          40    0.1581546      0.3864720       0.8826182
#> 22         40          40    0.1581546      0.3864720       0.8826182
#> 23         40          40    0.1581546      0.3864720       0.8826182
#> 24         40          40    0.1581546      0.3864720       0.8826182
#> 25         40          40    0.1581546      0.3864720       0.8826182
#> 26         40          40    0.1581546      0.3864720       0.8826182
#> 27         40          40    0.1581546      0.3864720       0.8826182
#> 28         40          40    0.1581546      0.3864720       0.8826182
#> 29         40          40    0.1581546      0.3864720       0.8826182
#> 30         40          40    0.1581546      0.3864720       0.8826182
#> 31         40          40    0.1581546      0.3864720       0.8826182
#> 32         40          40    0.1581546      0.3864720       0.8826182
#> 33         40          40    0.1581546      0.3864720       0.8826182
#> 34         40          40    0.1581546      0.3864720       0.8826182
#> 35         40          40    0.1581546      0.3864720       0.8826182
#> 36         40          40    0.1581546      0.3864720       0.8826182
#> 37         40          40    0.1581546      0.3864720       0.8826182
#> 38         40          40    0.1581546      0.3864720       0.8826182
#> 39         40          40    0.1581546      0.3864720       0.8826182
#> 40         40          40    0.1581546      0.3864720       0.8826182
#> 41         40          40    0.1581546      0.3864720       0.8826182
#> 42         40          40    0.1581546      0.3864720       0.8826182
#> 43         40          40    0.1581546      0.3864720       0.8826182
#> 44         40          40    0.1581546      0.3864720       0.8826182
#> 45         40          40    0.1581546      0.3864720       0.8826182
#> 46         40          40    0.1581546      0.3864720       0.8826182
#> 47         40          40    0.1581546      0.3864720       0.8826182
#> 48         40          40    0.1581546      0.3864720       0.8826182
#> 49         40          40    0.1581546      0.3864720       0.8826182
#> 50         40          40    0.1581546      0.3864720       0.8826182
#> 51         40          40    0.1581546      0.3864720       0.8826182
#> 52         40          40    0.1581546      0.3864720       0.8826182
#> 53         40          40    0.1581546      0.3864720       0.8826182
#> 54         40          40    0.1581546      0.3864720       0.8826182
#> 55         40          40    0.1581546      0.3864720       0.8826182
#> 56         40          40    0.1581546      0.3864720       0.8826182
#> 57         40          40    0.1581546      0.3864720       0.8826182
#> 58         40          40    0.1581546      0.3864720       0.8826182
#> 59         40          40    0.1581546      0.3864720       0.8826182
#> 60         40          40    0.1581546      0.3864720       0.8826182
#> 61         40          40    0.1581546      0.3864720       0.8826182
#> 62         40          40    0.1581546      0.3864720       0.8826182
#> 63         40          40    0.1581546      0.3864720       0.8826182
#> 64         40          40    0.1581546      0.3864720       0.8826182
#> 65         40          40    0.1581546      0.3864720       0.8826182
#> 66         40          40    0.1581546      0.3864720       0.8826182
#> 67         40          40    0.1581546      0.3864720       0.8826182
#> 68         40          40    0.1581546      0.3864720       0.8826182
#> 69         40          40    0.1581546      0.3864720       0.8826182
#> 70         40          40    0.1581546      0.3864720       0.8826182
#> 71         40          40    0.1581546      0.3864720       0.8826182
#> 72         40          40    0.1581546      0.3864720       0.8826182
#> 73         40          40    0.1581546      0.3864720       0.8826182
#> 74         40          40    0.1581546      0.3864720       0.8826182
#> 75         40          40    0.1581546      0.3864720       0.8826182
#> 76         40          40    0.1581546      0.3864720       0.8826182
#> 77         40          40    0.1581546      0.3864720       0.8826182
#> 78         40          40    0.1581546      0.3864720       0.8826182
#> 79         40          40    0.1581546      0.3864720       0.8826182
#> 80         40          40    0.1581546      0.3864720       0.8826182
#> 81         40          40    0.1669418      0.1654918       0.5048670
#> 82         40          40    0.1669418      0.1654918       0.5048670
#> 83         40          40    0.1669418      0.1654918       0.5048670
#> 84         40          40    0.1669418      0.1654918       0.5048670
#> 85         40          40    0.1669418      0.1654918       0.5048670
#> 86         40          40    0.1669418      0.1654918       0.5048670
#> 87         40          40    0.1669418      0.1654918       0.5048670
#> 88         40          40    0.1669418      0.1654918       0.5048670
#> 89         40          40    0.1669418      0.1654918       0.5048670
#> 90         40          40    0.1669418      0.1654918       0.5048670
#> 91         40          40    0.1669418      0.1654918       0.5048670
#> 92         40          40    0.1669418      0.1654918       0.5048670
#> 93         40          40    0.1669418      0.1654918       0.5048670
#> 94         40          40    0.1669418      0.1654918       0.5048670
#> 95         40          40    0.1669418      0.1654918       0.5048670
#> 96         40          40    0.1669418      0.1654918       0.5048670
#> 97         40          40    0.1669418      0.1654918       0.5048670
#> 98         40          40    0.1669418      0.1654918       0.5048670
#> 99         40          40    0.1669418      0.1654918       0.5048670
#> 100        40          40    0.1669418      0.1654918       0.5048670
#> 101        40          40    0.1669418      0.1654918       0.5048670
#> 102        40          40    0.1669418      0.1654918       0.5048670
#> 103        40          40    0.1669418      0.1654918       0.5048670
#> 104        40          40    0.1669418      0.1654918       0.5048670
#> 105        40          40    0.1669418      0.1654918       0.5048670
#> 106        40          40    0.1669418      0.1654918       0.5048670
#> 107        40          40    0.1669418      0.1654918       0.5048670
#> 108        40          40    0.1669418      0.1654918       0.5048670
#> 109        40          40    0.1669418      0.1654918       0.5048670
#> 110        40          40    0.1669418      0.1654918       0.5048670
#> 111        40          40    0.1669418      0.1654918       0.5048670
#> 112        40          40    0.1669418      0.1654918       0.5048670
#> 113        40          40    0.1669418      0.1654918       0.5048670
#> 114        40          40    0.1669418      0.1654918       0.5048670
#> 115        40          40    0.1669418      0.1654918       0.5048670
#> 116        40          40    0.1669418      0.1654918       0.5048670
#> 117        40          40    0.1669418      0.1654918       0.5048670
#> 118        40          40    0.1669418      0.1654918       0.5048670
#> 119        40          40    0.1669418      0.1654918       0.5048670
#> 120        40          40    0.1669418      0.1654918       0.5048670
#> 121        40          40    0.1669418      0.1654918       0.5048670
#> 122        40          40    0.1669418      0.1654918       0.5048670
#> 123        40          40    0.1669418      0.1654918       0.5048670
#> 124        40          40    0.1669418      0.1654918       0.5048670
#> 125        40          40    0.1669418      0.1654918       0.5048670
#> 126        40          40    0.1669418      0.1654918       0.5048670
#> 127        40          40    0.1669418      0.1654918       0.5048670
#> 128        40          40    0.1669418      0.1654918       0.5048670
#> 129        40          40    0.1669418      0.1654918       0.5048670
#> 130        40          40    0.1669418      0.1654918       0.5048670
#> 131        40          40    0.1669418      0.1654918       0.5048670
#> 132        40          40    0.1669418      0.1654918       0.5048670
#> 133        40          40    0.1669418      0.1654918       0.5048670
#> 134        40          40    0.1669418      0.1654918       0.5048670
#> 135        40          40    0.1669418      0.1654918       0.5048670
#> 136        40          40    0.1669418      0.1654918       0.5048670
#> 137        40          40    0.1669418      0.1654918       0.5048670
#> 138        40          40    0.1669418      0.1654918       0.5048670
#> 139        40          40    0.1669418      0.1654918       0.5048670
#> 140        40          40    0.1669418      0.1654918       0.5048670
#> 141        40          40    0.1669418      0.1654918       0.5048670
#> 142        40          40    0.1669418      0.1654918       0.5048670
#> 143        40          40    0.1669418      0.1654918       0.5048670
#> 144        40          40    0.1669418      0.1654918       0.5048670
#> 145        40          40    0.1669418      0.1654918       0.5048670
#> 146        40          40    0.1669418      0.1654918       0.5048670
#> 147        40          40    0.1669418      0.1654918       0.5048670
#> 148        40          40    0.1669418      0.1654918       0.5048670
#> 149        40          40    0.1669418      0.1654918       0.5048670
#> 150        40          40    0.1669418      0.1654918       0.5048670
#> 151        40          40    0.1669418      0.1654918       0.5048670
#> 152        40          40    0.1669418      0.1654918       0.5048670
#> 153        40          40    0.1669418      0.1654918       0.5048670
#> 154        40          40    0.1669418      0.1654918       0.5048670
#> 155        40          40    0.1669418      0.1654918       0.5048670
#> 156        40          40    0.1669418      0.1654918       0.5048670
#> 157        40          40    0.1669418      0.1654918       0.5048670
#> 158        40          40    0.1669418      0.1654918       0.5048670
#> 159        40          40    0.1669418      0.1654918       0.5048670
#> 160        40          40    0.1669418      0.1654918       0.5048670
#> 161        40          40    0.1669418      0.1654918       0.5048670
#> 162        40          40   -0.4034483      0.3106430       0.9825414
#> 163        40          40   -0.4034483      0.3106430       0.9825414
#> 164        40          40   -0.4034483      0.3106430       0.9825414
#> 165        40          40   -0.4034483      0.3106430       0.9825414
#> 166        40          40   -0.4034483      0.3106430       0.9825414
#> 167        40          40   -0.4034483      0.3106430       0.9825414
#> 168        40          40   -0.4034483      0.3106430       0.9825414
#> 169        40          40   -0.4034483      0.3106430       0.9825414
#> 170        40          40   -0.4034483      0.3106430       0.9825414
#> 171        40          40   -0.4034483      0.3106430       0.9825414
#> 172        40          40   -0.4034483      0.3106430       0.9825414
#> 173        40          40   -0.4034483      0.3106430       0.9825414
#> 174        40          40   -0.4034483      0.3106430       0.9825414
#> 175        40          40   -0.4034483      0.3106430       0.9825414
#> 176        40          40   -0.4034483      0.3106430       0.9825414
#> 177        40          40   -0.4034483      0.3106430       0.9825414
#> 178        40          40   -0.4034483      0.3106430       0.9825414
#> 179        40          40   -0.4034483      0.3106430       0.9825414
#> 180        40          40   -0.4034483      0.3106430       0.9825414
#> 181        40          40   -0.4034483      0.3106430       0.9825414
#> 182        40          40   -0.4034483      0.3106430       0.9825414
#> 183        40          40   -0.4034483      0.3106430       0.9825414
#> 184        40          40   -0.4034483      0.3106430       0.9825414
#> 185        40          40   -0.4034483      0.3106430       0.9825414
#> 186        40          40   -0.4034483      0.3106430       0.9825414
#> 187        40          40   -0.4034483      0.3106430       0.9825414
#> 188        40          40   -0.4034483      0.3106430       0.9825414
#> 189        40          40   -0.4034483      0.3106430       0.9825414
#> 190        40          40   -0.4034483      0.3106430       0.9825414
#> 191        40          40   -0.4034483      0.3106430       0.9825414
#> 192        40          40   -0.4034483      0.3106430       0.9825414
#> 193        40          40   -0.4034483      0.3106430       0.9825414
#> 194        40          40   -0.4034483      0.3106430       0.9825414
#> 195        40          40   -0.4034483      0.3106430       0.9825414
#> 196        40          40   -0.4034483      0.3106430       0.9825414
#> 197        40          40   -0.4034483      0.3106430       0.9825414
#> 198        40          40   -0.4034483      0.3106430       0.9825414
#> 199        40          40   -0.4034483      0.3106430       0.9825414
#> 200        40          40   -0.4034483      0.3106430       0.9825414
#> 201        40          40   -0.4034483      0.3106430       0.9825414
#> 202        40          40   -0.4034483      0.3106430       0.9825414
#> 203        40          40   -0.4034483      0.3106430       0.9825414
#> 204        40          40   -0.4034483      0.3106430       0.9825414
#> 205        40          40   -0.4034483      0.3106430       0.9825414
#> 206        40          40   -0.4034483      0.3106430       0.9825414
#> 207        40          40   -0.4034483      0.3106430       0.9825414
#> 208        40          40   -0.4034483      0.3106430       0.9825414
#> 209        40          40   -0.4034483      0.3106430       0.9825414
#> 210        40          40   -0.4034483      0.3106430       0.9825414
#> 211        40          40   -0.4034483      0.3106430       0.9825414
#> 212        40          40   -0.4034483      0.3106430       0.9825414
#> 213        40          40   -0.4034483      0.3106430       0.9825414
#> 214        40          40   -0.4034483      0.3106430       0.9825414
#> 215        40          40   -0.4034483      0.3106430       0.9825414
#> 216        40          40   -0.4034483      0.3106430       0.9825414
#> 217        40          40   -0.4034483      0.3106430       0.9825414
#> 218        40          40   -0.4034483      0.3106430       0.9825414
#> 219        40          40   -0.4034483      0.3106430       0.9825414
#> 220        40          40   -0.4034483      0.3106430       0.9825414
#> 221        40          40   -0.4034483      0.3106430       0.9825414
#> 222        40          40   -0.4034483      0.3106430       0.9825414
#> 223        40          40   -0.4034483      0.3106430       0.9825414
#> 224        40          40   -0.4034483      0.3106430       0.9825414
#> 225        40          40   -0.4034483      0.3106430       0.9825414
#> 226        40          40   -0.4034483      0.3106430       0.9825414
#> 227        40          40   -0.4034483      0.3106430       0.9825414
#> 228        40          40   -0.4034483      0.3106430       0.9825414
#> 229        40          40   -0.4034483      0.3106430       0.9825414
#> 230        40          40   -0.4034483      0.3106430       0.9825414
#> 231        40          40   -0.4034483      0.3106430       0.9825414
#> 232        40          40   -0.4034483      0.3106430       0.9825414
#> 233        40          40   -0.4034483      0.3106430       0.9825414
#> 234        40          40   -0.4034483      0.3106430       0.9825414
#> 235        40          40   -0.4034483      0.3106430       0.9825414
#> 236        40          40   -0.4034483      0.3106430       0.9825414
#> 237        40          40   -0.4034483      0.3106430       0.9825414
#> 238        40          40   -0.4034483      0.3106430       0.9825414
#> 239        40          40   -0.4034483      0.3106430       0.9825414
#> 240        40          40   -0.4034483      0.3106430       0.9825414
#> 241        40          40   -0.4034483      0.3106430       0.9825414
#> 242        40          40   -0.4034483      0.3106430       0.9825414
#>      delta_q025 delta_q500 delta_q975 decision_go
#> 1   -0.05183437  0.3272105  0.6968002           1
#> 2   -0.05183437  0.3272105  0.6968002           1
#> 3   -0.05183437  0.3272105  0.6968002           1
#> 4   -0.05183437  0.3272105  0.6968002           1
#> 5   -0.05183437  0.3272105  0.6968002           1
#> 6   -0.05183437  0.3272105  0.6968002           1
#> 7   -0.05183437  0.3272105  0.6968002           1
#> 8   -0.05183437  0.3272105  0.6968002           1
#> 9   -0.05183437  0.3272105  0.6968002           1
#> 10  -0.05183437  0.3272105  0.6968002           1
#> 11  -0.05183437  0.3272105  0.6968002           1
#> 12  -0.05183437  0.3272105  0.6968002           1
#> 13  -0.05183437  0.3272105  0.6968002           1
#> 14  -0.05183437  0.3272105  0.6968002           1
#> 15  -0.05183437  0.3272105  0.6968002           1
#> 16  -0.05183437  0.3272105  0.6968002           1
#> 17  -0.05183437  0.3272105  0.6968002           1
#> 18  -0.05183437  0.3272105  0.6968002           1
#> 19  -0.05183437  0.3272105  0.6968002           1
#> 20  -0.05183437  0.3272105  0.6968002           1
#> 21  -0.05183437  0.3272105  0.6968002           1
#> 22  -0.05183437  0.3272105  0.6968002           1
#> 23  -0.05183437  0.3272105  0.6968002           1
#> 24  -0.05183437  0.3272105  0.6968002           1
#> 25  -0.05183437  0.3272105  0.6968002           1
#> 26  -0.05183437  0.3272105  0.6968002           1
#> 27  -0.05183437  0.3272105  0.6968002           1
#> 28  -0.05183437  0.3272105  0.6968002           1
#> 29  -0.05183437  0.3272105  0.6968002           1
#> 30  -0.05183437  0.3272105  0.6968002           1
#> 31  -0.05183437  0.3272105  0.6968002           1
#> 32  -0.05183437  0.3272105  0.6968002           1
#> 33  -0.05183437  0.3272105  0.6968002           1
#> 34  -0.05183437  0.3272105  0.6968002           1
#> 35  -0.05183437  0.3272105  0.6968002           1
#> 36  -0.05183437  0.3272105  0.6968002           1
#> 37  -0.05183437  0.3272105  0.6968002           1
#> 38  -0.05183437  0.3272105  0.6968002           1
#> 39  -0.05183437  0.3272105  0.6968002           1
#> 40  -0.05183437  0.3272105  0.6968002           1
#> 41  -0.05183437  0.3272105  0.6968002           1
#> 42  -0.05183437  0.3272105  0.6968002           1
#> 43  -0.05183437  0.3272105  0.6968002           1
#> 44  -0.05183437  0.3272105  0.6968002           1
#> 45  -0.05183437  0.3272105  0.6968002           1
#> 46  -0.05183437  0.3272105  0.6968002           1
#> 47  -0.05183437  0.3272105  0.6968002           1
#> 48  -0.05183437  0.3272105  0.6968002           1
#> 49  -0.05183437  0.3272105  0.6968002           1
#> 50  -0.05183437  0.3272105  0.6968002           1
#> 51  -0.05183437  0.3272105  0.6968002           1
#> 52  -0.05183437  0.3272105  0.6968002           1
#> 53  -0.05183437  0.3272105  0.6968002           1
#> 54  -0.05183437  0.3272105  0.6968002           1
#> 55  -0.05183437  0.3272105  0.6968002           1
#> 56  -0.05183437  0.3272105  0.6968002           1
#> 57  -0.05183437  0.3272105  0.6968002           1
#> 58  -0.05183437  0.3272105  0.6968002           1
#> 59  -0.05183437  0.3272105  0.6968002           1
#> 60  -0.05183437  0.3272105  0.6968002           1
#> 61  -0.05183437  0.3272105  0.6968002           1
#> 62  -0.05183437  0.3272105  0.6968002           1
#> 63  -0.05183437  0.3272105  0.6968002           1
#> 64  -0.05183437  0.3272105  0.6968002           1
#> 65  -0.05183437  0.3272105  0.6968002           1
#> 66  -0.05183437  0.3272105  0.6968002           1
#> 67  -0.05183437  0.3272105  0.6968002           1
#> 68  -0.05183437  0.3272105  0.6968002           1
#> 69  -0.05183437  0.3272105  0.6968002           1
#> 70  -0.05183437  0.3272105  0.6968002           1
#> 71  -0.05183437  0.3272105  0.6968002           1
#> 72  -0.05183437  0.3272105  0.6968002           1
#> 73  -0.05183437  0.3272105  0.6968002           1
#> 74  -0.05183437  0.3272105  0.6968002           1
#> 75  -0.05183437  0.3272105  0.6968002           1
#> 76  -0.05183437  0.3272105  0.6968002           1
#> 77  -0.05183437  0.3272105  0.6968002           1
#> 78  -0.05183437  0.3272105  0.6968002           1
#> 79  -0.05183437  0.3272105  0.6968002           1
#> 80  -0.05183437  0.3272105  0.6968002           1
#> 81  -0.27784403  0.1023216  0.4720979           0
#> 82  -0.27784403  0.1023216  0.4720979           0
#> 83  -0.27784403  0.1023216  0.4720979           0
#> 84  -0.27784403  0.1023216  0.4720979           0
#> 85  -0.27784403  0.1023216  0.4720979           0
#> 86  -0.27784403  0.1023216  0.4720979           0
#> 87  -0.27784403  0.1023216  0.4720979           0
#> 88  -0.27784403  0.1023216  0.4720979           0
#> 89  -0.27784403  0.1023216  0.4720979           0
#> 90  -0.27784403  0.1023216  0.4720979           0
#> 91  -0.27784403  0.1023216  0.4720979           0
#> 92  -0.27784403  0.1023216  0.4720979           0
#> 93  -0.27784403  0.1023216  0.4720979           0
#> 94  -0.27784403  0.1023216  0.4720979           0
#> 95  -0.27784403  0.1023216  0.4720979           0
#> 96  -0.27784403  0.1023216  0.4720979           0
#> 97  -0.27784403  0.1023216  0.4720979           0
#> 98  -0.27784403  0.1023216  0.4720979           0
#> 99  -0.27784403  0.1023216  0.4720979           0
#> 100 -0.27784403  0.1023216  0.4720979           0
#> 101 -0.27784403  0.1023216  0.4720979           0
#> 102 -0.27784403  0.1023216  0.4720979           0
#> 103 -0.27784403  0.1023216  0.4720979           0
#> 104 -0.27784403  0.1023216  0.4720979           0
#> 105 -0.27784403  0.1023216  0.4720979           0
#> 106 -0.27784403  0.1023216  0.4720979           0
#> 107 -0.27784403  0.1023216  0.4720979           0
#> 108 -0.27784403  0.1023216  0.4720979           0
#> 109 -0.27784403  0.1023216  0.4720979           0
#> 110 -0.27784403  0.1023216  0.4720979           0
#> 111 -0.27784403  0.1023216  0.4720979           0
#> 112 -0.27784403  0.1023216  0.4720979           0
#> 113 -0.27784403  0.1023216  0.4720979           0
#> 114 -0.27784403  0.1023216  0.4720979           0
#> 115 -0.27784403  0.1023216  0.4720979           0
#> 116 -0.27784403  0.1023216  0.4720979           0
#> 117 -0.27784403  0.1023216  0.4720979           0
#> 118 -0.27784403  0.1023216  0.4720979           0
#> 119 -0.27784403  0.1023216  0.4720979           0
#> 120 -0.27784403  0.1023216  0.4720979           0
#> 121 -0.27784403  0.1023216  0.4720979           0
#> 122 -0.27784403  0.1023216  0.4720979           0
#> 123 -0.27784403  0.1023216  0.4720979           0
#> 124 -0.27784403  0.1023216  0.4720979           0
#> 125 -0.27784403  0.1023216  0.4720979           0
#> 126 -0.27784403  0.1023216  0.4720979           0
#> 127 -0.27784403  0.1023216  0.4720979           0
#> 128 -0.27784403  0.1023216  0.4720979           0
#> 129 -0.27784403  0.1023216  0.4720979           0
#> 130 -0.27784403  0.1023216  0.4720979           0
#> 131 -0.27784403  0.1023216  0.4720979           0
#> 132 -0.27784403  0.1023216  0.4720979           0
#> 133 -0.27784403  0.1023216  0.4720979           0
#> 134 -0.27784403  0.1023216  0.4720979           0
#> 135 -0.27784403  0.1023216  0.4720979           0
#> 136 -0.27784403  0.1023216  0.4720979           0
#> 137 -0.27784403  0.1023216  0.4720979           0
#> 138 -0.27784403  0.1023216  0.4720979           0
#> 139 -0.27784403  0.1023216  0.4720979           0
#> 140 -0.27784403  0.1023216  0.4720979           0
#> 141 -0.27784403  0.1023216  0.4720979           0
#> 142 -0.27784403  0.1023216  0.4720979           0
#> 143 -0.27784403  0.1023216  0.4720979           0
#> 144 -0.27784403  0.1023216  0.4720979           0
#> 145 -0.27784403  0.1023216  0.4720979           0
#> 146 -0.27784403  0.1023216  0.4720979           0
#> 147 -0.27784403  0.1023216  0.4720979           0
#> 148 -0.27784403  0.1023216  0.4720979           0
#> 149 -0.27784403  0.1023216  0.4720979           0
#> 150 -0.27784403  0.1023216  0.4720979           0
#> 151 -0.27784403  0.1023216  0.4720979           0
#> 152 -0.27784403  0.1023216  0.4720979           0
#> 153 -0.27784403  0.1023216  0.4720979           0
#> 154 -0.27784403  0.1023216  0.4720979           0
#> 155 -0.27784403  0.1023216  0.4720979           0
#> 156 -0.27784403  0.1023216  0.4720979           0
#> 157 -0.27784403  0.1023216  0.4720979           0
#> 158 -0.27784403  0.1023216  0.4720979           0
#> 159 -0.27784403  0.1023216  0.4720979           0
#> 160 -0.27784403  0.1023216  0.4720979           0
#> 161 -0.27784403  0.1023216  0.4720979           0
#> 162  0.12871903  0.5184149  0.9892527           1
#> 163  0.12871903  0.5184149  0.9892527           1
#> 164  0.12871903  0.5184149  0.9892527           1
#> 165  0.12871903  0.5184149  0.9892527           1
#> 166  0.12871903  0.5184149  0.9892527           1
#> 167  0.12871903  0.5184149  0.9892527           1
#> 168  0.12871903  0.5184149  0.9892527           1
#> 169  0.12871903  0.5184149  0.9892527           1
#> 170  0.12871903  0.5184149  0.9892527           1
#> 171  0.12871903  0.5184149  0.9892527           1
#> 172  0.12871903  0.5184149  0.9892527           1
#> 173  0.12871903  0.5184149  0.9892527           1
#> 174  0.12871903  0.5184149  0.9892527           1
#> 175  0.12871903  0.5184149  0.9892527           1
#> 176  0.12871903  0.5184149  0.9892527           1
#> 177  0.12871903  0.5184149  0.9892527           1
#> 178  0.12871903  0.5184149  0.9892527           1
#> 179  0.12871903  0.5184149  0.9892527           1
#> 180  0.12871903  0.5184149  0.9892527           1
#> 181  0.12871903  0.5184149  0.9892527           1
#> 182  0.12871903  0.5184149  0.9892527           1
#> 183  0.12871903  0.5184149  0.9892527           1
#> 184  0.12871903  0.5184149  0.9892527           1
#> 185  0.12871903  0.5184149  0.9892527           1
#> 186  0.12871903  0.5184149  0.9892527           1
#> 187  0.12871903  0.5184149  0.9892527           1
#> 188  0.12871903  0.5184149  0.9892527           1
#> 189  0.12871903  0.5184149  0.9892527           1
#> 190  0.12871903  0.5184149  0.9892527           1
#> 191  0.12871903  0.5184149  0.9892527           1
#> 192  0.12871903  0.5184149  0.9892527           1
#> 193  0.12871903  0.5184149  0.9892527           1
#> 194  0.12871903  0.5184149  0.9892527           1
#> 195  0.12871903  0.5184149  0.9892527           1
#> 196  0.12871903  0.5184149  0.9892527           1
#> 197  0.12871903  0.5184149  0.9892527           1
#> 198  0.12871903  0.5184149  0.9892527           1
#> 199  0.12871903  0.5184149  0.9892527           1
#> 200  0.12871903  0.5184149  0.9892527           1
#> 201  0.12871903  0.5184149  0.9892527           1
#> 202  0.12871903  0.5184149  0.9892527           1
#> 203  0.12871903  0.5184149  0.9892527           1
#> 204  0.12871903  0.5184149  0.9892527           1
#> 205  0.12871903  0.5184149  0.9892527           1
#> 206  0.12871903  0.5184149  0.9892527           1
#> 207  0.12871903  0.5184149  0.9892527           1
#> 208  0.12871903  0.5184149  0.9892527           1
#> 209  0.12871903  0.5184149  0.9892527           1
#> 210  0.12871903  0.5184149  0.9892527           1
#> 211  0.12871903  0.5184149  0.9892527           1
#> 212  0.12871903  0.5184149  0.9892527           1
#> 213  0.12871903  0.5184149  0.9892527           1
#> 214  0.12871903  0.5184149  0.9892527           1
#> 215  0.12871903  0.5184149  0.9892527           1
#> 216  0.12871903  0.5184149  0.9892527           1
#> 217  0.12871903  0.5184149  0.9892527           1
#> 218  0.12871903  0.5184149  0.9892527           1
#> 219  0.12871903  0.5184149  0.9892527           1
#> 220  0.12871903  0.5184149  0.9892527           1
#> 221  0.12871903  0.5184149  0.9892527           1
#> 222  0.12871903  0.5184149  0.9892527           1
#> 223  0.12871903  0.5184149  0.9892527           1
#> 224  0.12871903  0.5184149  0.9892527           1
#> 225  0.12871903  0.5184149  0.9892527           1
#> 226  0.12871903  0.5184149  0.9892527           1
#> 227  0.12871903  0.5184149  0.9892527           1
#> 228  0.12871903  0.5184149  0.9892527           1
#> 229  0.12871903  0.5184149  0.9892527           1
#> 230  0.12871903  0.5184149  0.9892527           1
#> 231  0.12871903  0.5184149  0.9892527           1
#> 232  0.12871903  0.5184149  0.9892527           1
#> 233  0.12871903  0.5184149  0.9892527           1
#> 234  0.12871903  0.5184149  0.9892527           1
#> 235  0.12871903  0.5184149  0.9892527           1
#> 236  0.12871903  0.5184149  0.9892527           1
#> 237  0.12871903  0.5184149  0.9892527           1
#> 238  0.12871903  0.5184149  0.9892527           1
#> 239  0.12871903  0.5184149  0.9892527           1
#> 240  0.12871903  0.5184149  0.9892527           1
#> 241  0.12871903  0.5184149  0.9892527           1
#> 242  0.12871903  0.5184149  0.9892527           1
```

`decision_go = 1` indicates **Go**, and `decision_go = 0` indicates
**No-Go** under the criterion
$P\left( \Delta > 0.1 \mid \text{data} \right) \geq 0.8$.

# Example 4: Two-arm \| Fixed design \| Single continuous endpoint \| Bayesian Go/No-Go with placebo borrowing

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
(delta \> threshold) exceeds a predefined threshold.

**Unique focus:** Bayesian Go/No-Go criterion, historical data borrowing
via
[`RBesT::robustify()`](https://opensource.nibr.com/RBesT/reference/robustify.html),
Go probability as a function of true effect size.

The setup pattern stays aligned with [Example
1](https://boehringer-ingelheim.github.io/rxsim/articles/example-1.md),
while the analysis is replaced by a Bayesian decision rule with
borrowing.

## Scenario

A **two-arm, fixed design** trial (placebo vs treatment) with a **single
continuous endpoint** and a **single final analysis**.

- Placebo arm uses **historical borrowing** through an informative prior
  derived from historical placebo data and robustified with a
  non-informative component.
- Treatment arm uses a **non-informative prior**.
- Go/No-Go rule: Go if P(mu_T - mu_P \> delta \| data) \>= gamma.

`delta = 0.1` defines the minimum treatment-placebo difference of
clinical relevance. `gamma = 0.8` requires 80% posterior probability of
exceeding that threshold before a Go decision is issued. The true
simulated effect `delta_true = 0.30` lies above `delta`, so a
well-powered trial should frequently result in a Go.

``` r
# Trial design
sample_size <- 80
allocation  <- c(1, 1)
arms        <- c("placebo", "treatment")

# Data-generating truth for this simulation example
mu_placebo_true <- 0.00
delta_true      <- 0.30  # true treatment - placebo mean difference
sigma_known     <- 1.0

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

enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn    <- function(n) rexp(n, rate = 0.01)
```

## Historical placebo dataset and priors

We include a small historical placebo dataset as study-level means and
sample sizes.

Studies H1, H2, and H3 contribute a combined 62 historical placebo
observations and are summarised into an informative prior via
[`postmix()`](https://opensource.nibr.com/RBesT/reference/postmix.html) -
a Bayesian update of the non-informative starting prior with the pooled
historical data.
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

## Conditions

`RBesT::pmixdiff(post_treat, post_placebo, delta, lower.tail = FALSE)`
evaluates P(mu_T - mu_P \> delta \| data) by numerical integration over
the mixture posterior distributions for treatment and placebo.
[`qmixdiff()`](https://opensource.nibr.com/RBesT/reference/mixdiff.html)
computes quantiles of the same posterior difference distribution,
yielding a 95% credible interval for Delta = mu_T - mu_P.

``` r
analysis_generators <- list(
  final = list(
    trigger = enroll_trigger(1.0, sample_size),
    analysis = function(df, current_time) {
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

## Simulate

``` r
set.seed(7)
trials <- replicate_trial(
  trial_name            = "example_4_bayes_two_arm_fixed",
  sample_size           = sample_size,
  arms                  = arms,
  allocation            = allocation,
  enrollment            = enrollment_fn,
  dropout               = dropout_fn,
  analysis_generators   = analysis_generators,
  population_generators = population_generators,
  n                     = 5
)
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
5](https://boehringer-ingelheim.github.io/rxsim/articles/example-5.md)
for a seamless design that builds on this Bayesian decision rule.

``` r
replicate_results <- collect_results(trials)
replicate_results
#>   replicate timepoint analysis sample_size allocation delta_true delta gamma
#> 1         1  81.22545    final          80       1, 1        0.3   0.1   0.8
#> 2         2  85.38045    final          80       1, 1        0.3   0.1   0.8
#> 3         3  74.54773    final          80       1, 1        0.3   0.1   0.8
#> 4         4  80.51910    final          80       1, 1        0.3   0.1   0.8
#> 5         5  82.06812    final          80       1, 1        0.3   0.1   0.8
#>   n_placebo n_treatment mean_placebo mean_treatment post_prob_delta
#> 1        40          40  -0.02043065      0.1359496       0.6096597
#> 2        40          40  -0.14954828      0.5987274       0.9987468
#> 3        40          40  -0.16555987      0.2343530       0.8695478
#> 4        40          40   0.21295567      0.1868780       0.5054026
#> 5        40          40   0.01465888      0.3709480       0.9260954
#>     delta_q025 delta_q500 delta_q975 decision_go
#> 1 -0.216985237  0.1523456  0.5218935           0
#> 2  0.300417230  0.6694469  1.0447389           1
#> 3 -0.057222457  0.3120097  0.6887534           1
#> 4 -0.285650120  0.1026001  0.4737831           0
#> 5  0.002710835  0.3727511  0.7417912           1
```

`decision_go = 1` indicates Go; `decision_go = 0` indicates No-Go under
the criterion P(Delta \> 0.1 \| data) \>= 0.8.

## Go probability curve

The operating characteristics of a Bayesian Go/No-Go design are
summarised by the Go probability under a range of true treatment
effects: near zero for small effects (false Go = type I error analog)
and approaching 1 for large effects.

``` r
set.seed(77)
n_reps_pw   <- 200
delta_sweep <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)

go_df <- do.call(rbind, lapply(delta_sweep, function(dt) {
  pop_pw <- list(
    placebo = function(n) {
      data.frame(
        id = 1:n, y = rnorm(n, 0, sigma_known), readout_time = 1
      )
    },
    treatment = local({
      d <- dt
      function(n) {
        data.frame(
          id = 1:n,
          y  = rnorm(n, d, sigma_known),
          readout_time = 1
        )
      }
    })
  )
  an_pw <- list(final = list(
    trigger  = enroll_trigger(1.0, sample_size),
    analysis = function(df, ct) {
      dat <- dplyr::filter(df, !is.na(enroll_time)) |>
        dplyr::mutate(arm = factor(arm, levels = c("placebo", "treatment")))
      y_p <- dat$y[dat$arm == "placebo"]
      y_t <- dat$y[dat$arm == "treatment"]
      pp  <- RBesT::postmix(prior_placebo,   n = length(y_p), m = mean(y_p))
      pt  <- RBesT::postmix(prior_treatment, n = length(y_t), m = mean(y_t))
      data.frame(
        go = as.integer(
          RBesT::pmixdiff(pt, pp, delta, lower.tail = FALSE) >= gamma
        )
      )
    }
  ))
  tr <- replicate_trial(
    "pw7", sample_size, arms, allocation,
    enrollment_fn, dropout_fn, an_pw, pop_pw, n_reps_pw
  )
  invisible(run_trials(tr))
  data.frame(
    delta_true   = dt,
    go_prob      = mean(collect_results(tr)$go)
  )
}))
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
#> Using default prior reference scale 1
```

``` r
plot(go_df$delta_true, go_df$go_prob, type = "b", pch = 19,
     col = "steelblue",
     xlab = "True treatment - placebo difference",
     ylab = "Go probability",
     main = sprintf("Go/No-Go operating characteristics (n = %d)", sample_size),
     ylim = c(0, 1))
abline(v = delta, lty = 3, col = "grey50")  # decision threshold
abline(h = gamma, lty = 2, col = "grey50")  # gamma threshold
legend("topleft",
       legend = c(
         sprintf("delta threshold = %.1f", delta),
         sprintf("gamma threshold = %.1f", gamma)
       ),
       lty = c(3, 2), col = "grey50")
```

![](example-4_files/figure-html/plot_go-1.png)

## Next steps

- [Example
  5](https://boehringer-ingelheim.github.io/rxsim/articles/example-5.md) -
  seamless Phase IIa/IIb design using BayesianMCPMod
- [Conditions and
  Triggers](https://boehringer-ingelheim.github.io/rxsim/articles/conditions.md) -
  adding interim looks to a Bayesian design

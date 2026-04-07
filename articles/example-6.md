# Example 6: Two-arm \| Fixed design \| Single continuous endpoint \| Overall + exploratory subgroup analysis

``` r
# Core simulation framework
library(rxsim)

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

set.seed(6606)
```

Subgroup analyses are a standard exploratory component of clinical trial
reporting. This example simulates a trial where the treatment effect
differs between two patient subgroups — a form of treatment effect
heterogeneity. We evaluate the overall treatment effect,
subgroup-specific effects, and formally test whether the arm-by-subgroup
interaction is statistically significant.

## Scenario

A **two-arm, fixed design** trial (placebo vs treatment) with a **single
continuous endpoint** and one **final analysis**.

At final analysis, we report:

- **Overall analysis**: treatment effect across all subjects.
- **Exploratory subgroup analysis**: treatment effect within each
  subgroup + arm-by-subgroup interaction test.

Treatment effect heterogeneity arises when patients with different
baseline characteristics — such as disease subtype or biomarker status —
respond differently to treatment. Here, `delta_A = 0.40` and
`delta_B = 0.15` encode a larger benefit in subgroup A; treatment is
effective in both subgroups, but considerably more so in A.
`subgroup_prob = c(0.6, 0.4)` means that on average 60% of enrolled
patients belong to subgroup A and 40% to subgroup B, though actual
proportions vary across replicates due to random assignment.

``` r
# Trial design
sample_size <- 120
allocation  <- c(1, 1)
arms        <- c("placebo", "treatment")

# Subgroup settings
subgroup_levels <- c("A", "B")
subgroup_prob   <- c(0.6, 0.4)

# Data-generating truth (example values)
mu_placebo_A <- 0.00
mu_placebo_B <- 0.00
delta_A      <- 0.40  # true treatment effect in subgroup A
delta_B      <- 0.15  # true treatment effect in subgroup B
sigma        <- 1.0

# Enrollment/dropout profile
enrollment <- list(
  end_time = c(4, 8),
  rate     = c(8, 7)
)
dropout <- list(
  end_time = c(4, 8),
  rate     = c(0, 0)
)

scenario <- tidyr::expand_grid(
  sample_size  = sample_size,
  allocation   = list(allocation),
  delta_A      = delta_A,
  delta_B      = delta_B,
  p_subgroup_A = subgroup_prob[1],
  p_subgroup_B = subgroup_prob[2]
)
```

## Time points

[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
uses a two-interval piecewise-constant schedule (rates 8 and 7 per time
unit) to reach the target sample size of 120. Both dropout rates are
zero, so all enrolled subjects complete the trial and contribute to the
final analysis.

``` r
timepoints <- gen_timepoints(
  sample_size = sample_size,
  arms        = arms,
  allocation  = allocation,
  enrollment  = enrollment,
  dropout     = dropout
)

tr_timer <- Timer$new(name = "timer_example_6")
add_timepoints(tr_timer, timepoints)

final_time <- tr_timer$get_end_timepoint()
final_time
#> [1] 9
```

## Populations

Create two populations and include a `subgroup` column in each
population dataset.

In each generator, subgroup membership is randomly drawn for each
subject with probabilities `c(0.6, 0.4)`, and the arm-specific mean is
selected conditional on that assignment. Because subgroup is sampled
independently in each replicate, realized subgroup counts — and hence
per-subgroup sample sizes — vary across replicates, reflecting a
realistic pre-stratified (but not block-randomised) design.

``` r
# Create generator for arm-specific data with subgroups
population_generators <- list(
  placebo = function(n) {
    subgroup <- sample(subgroup_levels, size = n, replace = TRUE, prob = subgroup_prob)
    mu <- ifelse(subgroup == "A", mu_placebo_A, mu_placebo_B)
    data.frame(
      id = seq_len(n),
      subgroup = subgroup,
      y = rnorm(n, mean = mu, sd = sigma),
      readout_time = 1
    )
  },
  treatment = function(n) {
    subgroup <- sample(subgroup_levels, size = n, replace = TRUE, prob = subgroup_prob)
    mu <- ifelse(subgroup == "A", mu_placebo_A + delta_A, mu_placebo_B + delta_B)
    data.frame(
      id = seq_len(n),
      subgroup = subgroup,
      y = rnorm(n, mean = mu, sd = sigma),
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

## Final analysis trigger

`lm(y ~ arm)` estimates the overall treatment effect across all
subjects. Per-subgroup effects are obtained by fitting `lm(y ~ arm)`
separately within each subgroup subset. The interaction model
`lm(y ~ arm * subgroup)` includes an `arm:subgroupB` coefficient that
tests directly whether the treatment effect in subgroup B differs from
that in the reference subgroup A — a statistically significant result
indicates genuine treatment effect heterogeneity.

``` r
analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer) {
      dat <- df |>
        dplyr::filter(!is.na(enroll_time)) |>
        dplyr::mutate(
          arm = factor(arm, levels = c("placebo", "treatment")),
          subgroup = factor(subgroup, levels = subgroup_levels)
        )

      fit_overall <- lm(y ~ arm, data = dat)
      coef_overall <- summary(fit_overall)$coefficients

      overall_est <- unname(coef_overall["armtreatment", "Estimate"])
      overall_p   <- unname(coef_overall["armtreatment", "Pr(>|t|)"])

      subgroup_stats <- lapply(levels(dat$subgroup), function(sg) {
        dsg <- dat[dat$subgroup == sg, , drop = FALSE]
        fit_sg <- lm(y ~ arm, data = dsg)
        coef_sg <- summary(fit_sg)$coefficients

        c(
          n_total = nrow(dsg),
          estimate = unname(coef_sg["armtreatment", "Estimate"]),
          p_value = unname(coef_sg["armtreatment", "Pr(>|t|)"])
        )
      })
      names(subgroup_stats) <- levels(dat$subgroup)

      fit_interaction <- lm(y ~ arm * subgroup, data = dat)
      coef_int <- summary(fit_interaction)$coefficients
      int_term <- grep("^arm.*:subgroup", rownames(coef_int), value = TRUE)
      p_interaction <- if (length(int_term) > 0) {
        unname(coef_int[int_term[1], "Pr(>|t|)"])
      } else {
        NA_real_
      }

      data.frame(
        scenario,
        n_total = nrow(dat),
        overall_estimate = overall_est,
        overall_p_value = overall_p,
        subgroup_A_n = unname(subgroup_stats[["A"]]["n_total"]),
        subgroup_A_estimate = unname(subgroup_stats[["A"]]["estimate"]),
        subgroup_A_p_value = unname(subgroup_stats[["A"]]["p_value"]),
        subgroup_B_n = unname(subgroup_stats[["B"]]["n_total"]),
        subgroup_B_estimate = unname(subgroup_stats[["B"]]["estimate"]),
        subgroup_B_p_value = unname(subgroup_stats[["B"]]["p_value"]),
        interaction_p_value = p_interaction,
        stringsAsFactors = FALSE
      )
    }
  )
)
```

## Trial

``` r
trials <- replicate_trial(
  trial_name = "example_6_two_arm_fixed_subgroup",
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
#>     name: example_6_two_arm_fixed_subgroup_1
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
#>     name: example_6_two_arm_fixed_subgroup_2
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
#>     name: example_6_two_arm_fixed_subgroup_3
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
`overall_estimate` and `overall_p_value` summarise the pooled treatment
effect; the per-subgroup columns reveal within-group signals. Note that
with approximately 24 patients per arm in subgroup B (40% of 60
treated), subgroup-level power is limited and p-values will be highly
variable — this is intentional, since exploratory subgroup analyses are
hypothesis-generating and should not support confirmatory conclusions. A
small `interaction_p_value` suggests genuine treatment effect
heterogeneity, but interpretation should be made cautiously given the
sample sizes. See [Getting
Started](https://boehringer-ingelheim.github.io/rxsim/articles/getting-started.md)
for simulation setup and [Core
Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.md)
for background.

``` r
replicate_results <- collect_results(trials)
replicate_results
#>   replicate timepoint analysis sample_size allocation delta_A delta_B
#> 1         1  126.9942    final         120       1, 1     0.4    0.15
#> 2         2  116.8608    final         120       1, 1     0.4    0.15
#> 3         3  117.7798    final         120       1, 1     0.4    0.15
#>   p_subgroup_A p_subgroup_B n_total overall_estimate overall_p_value
#> 1          0.6          0.4     120        0.3414139      0.06906957
#> 2          0.6          0.4     120        0.1649366      0.37184435
#> 3          0.6          0.4     120        0.3223337      0.08392587
#>   subgroup_A_n subgroup_A_estimate subgroup_A_p_value subgroup_B_n
#> 1           76           0.4707117         0.05403434           44
#> 2           64           0.3320609         0.18188722           56
#> 3           71           0.1774698         0.50062690           49
#>   subgroup_B_estimate subgroup_B_p_value interaction_p_value
#> 1          0.11808122         0.68769806           0.3632263
#> 2         -0.01071538         0.96958098           0.3573194
#> 3          0.45256210         0.06346625           0.4623181
```

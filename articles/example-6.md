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

## Scenario

A **two-arm, fixed design** trial (placebo vs treatment) with a **single
continuous endpoint** and one **final analysis**.

At final analysis, we report:

- **Overall analysis**: treatment effect across all subjects.
- **Exploratory subgroup analysis**: treatment effect within each
  subgroup + arm-by-subgroup interaction test.

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
mu_treat_A   <- 0.40
mu_treat_B   <- 0.15
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
  sample_size = sample_size,
  allocation = list(allocation),
  p_subgroup_A = subgroup_prob[1],
  p_subgroup_B = subgroup_prob[2]
)
```

## Time points

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
    mu <- ifelse(subgroup == "A", mu_treat_A, mu_treat_B)
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

``` r
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn <- function(n) rexp(n, rate = 0.01)
```

## Final analysis trigger

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

``` r
replicate_results <- do.call(rbind, lapply(seq_along(trials), function(i) {
  data.frame(
    replicate = i,
    trials[[i]]$results[[1]]$final,
    stringsAsFactors = FALSE
  )
}))

replicate_results
#>   replicate sample_size allocation p_subgroup_A p_subgroup_B n_total
#> 1         1         120       1, 1          0.6          0.4     120
#> 2         2         120       1, 1          0.6          0.4     120
#> 3         3         120       1, 1          0.6          0.4     120
#>   overall_estimate overall_p_value subgroup_A_n subgroup_A_estimate
#> 1        0.3414139      0.06906957           76           0.4707117
#> 2        0.1649366      0.37184435           64           0.3320609
#> 3        0.3223337      0.08392587           71           0.1774698
#>   subgroup_A_p_value subgroup_B_n subgroup_B_estimate subgroup_B_p_value
#> 1         0.05403434           44          0.11808122         0.68769806
#> 2         0.18188722           56         -0.01071538         0.96958098
#> 3         0.50062690           49          0.45256210         0.06346625
#>   interaction_p_value
#> 1           0.3632263
#> 2           0.3573194
#> 3           0.4623181
```

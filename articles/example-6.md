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
#>     replicate  timepoint analysis sample_size allocation delta_A delta_B
#> 1           1   126.9942    final         120       1, 1     0.4    0.15
#> 2           1   159.9026    final         120       1, 1     0.4    0.15
#> 3           1   262.9726    final         120       1, 1     0.4    0.15
#> 4           1   360.8773    final         120       1, 1     0.4    0.15
#> 5           1   418.0143    final         120       1, 1     0.4    0.15
#> 6           1   705.6521    final         120       1, 1     0.4    0.15
#> 7           1   714.8538    final         120       1, 1     0.4    0.15
#> 8           1   859.4674    final         120       1, 1     0.4    0.15
#> 9           1   879.8699    final         120       1, 1     0.4    0.15
#> 10          1   957.4728    final         120       1, 1     0.4    0.15
#> 11          1  1205.5393    final         120       1, 1     0.4    0.15
#> 12          1  1241.9691    final         120       1, 1     0.4    0.15
#> 13          1  1344.9184    final         120       1, 1     0.4    0.15
#> 14          1  1388.0477    final         120       1, 1     0.4    0.15
#> 15          1  1411.8025    final         120       1, 1     0.4    0.15
#> 16          1  1568.4174    final         120       1, 1     0.4    0.15
#> 17          1  1796.6362    final         120       1, 1     0.4    0.15
#> 18          1  1798.9369    final         120       1, 1     0.4    0.15
#> 19          1  1914.8075    final         120       1, 1     0.4    0.15
#> 20          1  1946.1549    final         120       1, 1     0.4    0.15
#> 21          1  1955.1874    final         120       1, 1     0.4    0.15
#> 22          1  2080.9565    final         120       1, 1     0.4    0.15
#> 23          1  2104.2159    final         120       1, 1     0.4    0.15
#> 24          1  2105.4026    final         120       1, 1     0.4    0.15
#> 25          1  2336.8627    final         120       1, 1     0.4    0.15
#> 26          1  2948.8485    final         120       1, 1     0.4    0.15
#> 27          1  3095.3305    final         120       1, 1     0.4    0.15
#> 28          1  3162.5072    final         120       1, 1     0.4    0.15
#> 29          1  3164.0664    final         120       1, 1     0.4    0.15
#> 30          1  3213.6982    final         120       1, 1     0.4    0.15
#> 31          1  3484.4899    final         120       1, 1     0.4    0.15
#> 32          1  3517.2528    final         120       1, 1     0.4    0.15
#> 33          1  3926.6322    final         120       1, 1     0.4    0.15
#> 34          1  4033.4085    final         120       1, 1     0.4    0.15
#> 35          1  4116.4265    final         120       1, 1     0.4    0.15
#> 36          1  4147.9365    final         120       1, 1     0.4    0.15
#> 37          1  4278.2886    final         120       1, 1     0.4    0.15
#> 38          1  4442.8165    final         120       1, 1     0.4    0.15
#> 39          1  4521.5402    final         120       1, 1     0.4    0.15
#> 40          1  4629.5562    final         120       1, 1     0.4    0.15
#> 41          1  4677.8727    final         120       1, 1     0.4    0.15
#> 42          1  4703.9073    final         120       1, 1     0.4    0.15
#> 43          1  4910.0092    final         120       1, 1     0.4    0.15
#> 44          1  5097.2426    final         120       1, 1     0.4    0.15
#> 45          1  5166.9992    final         120       1, 1     0.4    0.15
#> 46          1  5194.1642    final         120       1, 1     0.4    0.15
#> 47          1  5249.6239    final         120       1, 1     0.4    0.15
#> 48          1  5292.6248    final         120       1, 1     0.4    0.15
#> 49          1  5320.8122    final         120       1, 1     0.4    0.15
#> 50          1  5361.6571    final         120       1, 1     0.4    0.15
#> 51          1  5506.0544    final         120       1, 1     0.4    0.15
#> 52          1  5522.3163    final         120       1, 1     0.4    0.15
#> 53          1  5549.5819    final         120       1, 1     0.4    0.15
#> 54          1  5913.5980    final         120       1, 1     0.4    0.15
#> 55          1  5957.2678    final         120       1, 1     0.4    0.15
#> 56          1  6003.2059    final         120       1, 1     0.4    0.15
#> 57          1  6096.8808    final         120       1, 1     0.4    0.15
#> 58          1  6181.6053    final         120       1, 1     0.4    0.15
#> 59          1  6278.8775    final         120       1, 1     0.4    0.15
#> 60          1  6322.3404    final         120       1, 1     0.4    0.15
#> 61          1  6378.0450    final         120       1, 1     0.4    0.15
#> 62          1  6419.6931    final         120       1, 1     0.4    0.15
#> 63          1  6470.2407    final         120       1, 1     0.4    0.15
#> 64          1  6477.9458    final         120       1, 1     0.4    0.15
#> 65          1  6523.5139    final         120       1, 1     0.4    0.15
#> 66          1  6737.0848    final         120       1, 1     0.4    0.15
#> 67          1  6747.2130    final         120       1, 1     0.4    0.15
#> 68          1  6755.3474    final         120       1, 1     0.4    0.15
#> 69          1  6874.2393    final         120       1, 1     0.4    0.15
#> 70          1  6895.6774    final         120       1, 1     0.4    0.15
#> 71          1  6914.2226    final         120       1, 1     0.4    0.15
#> 72          1  7018.7020    final         120       1, 1     0.4    0.15
#> 73          1  7056.0540    final         120       1, 1     0.4    0.15
#> 74          1  7079.4071    final         120       1, 1     0.4    0.15
#> 75          1  7125.7133    final         120       1, 1     0.4    0.15
#> 76          1  7196.7340    final         120       1, 1     0.4    0.15
#> 77          1  7219.2676    final         120       1, 1     0.4    0.15
#> 78          1  7235.6936    final         120       1, 1     0.4    0.15
#> 79          1  7276.6637    final         120       1, 1     0.4    0.15
#> 80          1  7327.4466    final         120       1, 1     0.4    0.15
#> 81          1  7341.7815    final         120       1, 1     0.4    0.15
#> 82          1  7589.8428    final         120       1, 1     0.4    0.15
#> 83          1  7688.1042    final         120       1, 1     0.4    0.15
#> 84          1  7689.3229    final         120       1, 1     0.4    0.15
#> 85          1  7753.7958    final         120       1, 1     0.4    0.15
#> 86          1  7960.9844    final         120       1, 1     0.4    0.15
#> 87          1  8070.0121    final         120       1, 1     0.4    0.15
#> 88          1  8252.8266    final         120       1, 1     0.4    0.15
#> 89          1  8292.5254    final         120       1, 1     0.4    0.15
#> 90          1  8307.8861    final         120       1, 1     0.4    0.15
#> 91          1  8405.1250    final         120       1, 1     0.4    0.15
#> 92          1  8453.0335    final         120       1, 1     0.4    0.15
#> 93          1  8574.1773    final         120       1, 1     0.4    0.15
#> 94          1  8579.9216    final         120       1, 1     0.4    0.15
#> 95          1  8599.0932    final         120       1, 1     0.4    0.15
#> 96          1  8612.4215    final         120       1, 1     0.4    0.15
#> 97          1  8638.1610    final         120       1, 1     0.4    0.15
#> 98          1  8651.0329    final         120       1, 1     0.4    0.15
#> 99          1  8663.7031    final         120       1, 1     0.4    0.15
#> 100         1  8738.0257    final         120       1, 1     0.4    0.15
#> 101         1  8744.8145    final         120       1, 1     0.4    0.15
#> 102         1  8787.4906    final         120       1, 1     0.4    0.15
#> 103         1  8790.0346    final         120       1, 1     0.4    0.15
#> 104         1  8958.1746    final         120       1, 1     0.4    0.15
#> 105         1  8984.1668    final         120       1, 1     0.4    0.15
#> 106         1  9051.0283    final         120       1, 1     0.4    0.15
#> 107         1  9205.8030    final         120       1, 1     0.4    0.15
#> 108         1  9227.6409    final         120       1, 1     0.4    0.15
#> 109         1  9285.0251    final         120       1, 1     0.4    0.15
#> 110         1  9391.7396    final         120       1, 1     0.4    0.15
#> 111         1  9430.9920    final         120       1, 1     0.4    0.15
#> 112         1  9626.5986    final         120       1, 1     0.4    0.15
#> 113         1  9632.8382    final         120       1, 1     0.4    0.15
#> 114         1  9882.6013    final         120       1, 1     0.4    0.15
#> 115         1  9943.1782    final         120       1, 1     0.4    0.15
#> 116         1  9972.5609    final         120       1, 1     0.4    0.15
#> 117         1  9994.8197    final         120       1, 1     0.4    0.15
#> 118         1 10127.0338    final         120       1, 1     0.4    0.15
#> 119         1 10140.1540    final         120       1, 1     0.4    0.15
#> 120         2   116.8608    final         120       1, 1     0.4    0.15
#> 121         2   118.3699    final         120       1, 1     0.4    0.15
#> 122         2   152.9387    final         120       1, 1     0.4    0.15
#> 123         2   266.3100    final         120       1, 1     0.4    0.15
#> 124         2   315.5000    final         120       1, 1     0.4    0.15
#> 125         2   316.3716    final         120       1, 1     0.4    0.15
#> 126         2   383.3559    final         120       1, 1     0.4    0.15
#> 127         2   545.9930    final         120       1, 1     0.4    0.15
#> 128         2   574.4208    final         120       1, 1     0.4    0.15
#> 129         2   919.6843    final         120       1, 1     0.4    0.15
#> 130         2   937.3382    final         120       1, 1     0.4    0.15
#> 131         2  1235.2473    final         120       1, 1     0.4    0.15
#> 132         2  1400.4183    final         120       1, 1     0.4    0.15
#> 133         2  1471.3563    final         120       1, 1     0.4    0.15
#> 134         2  1720.3584    final         120       1, 1     0.4    0.15
#> 135         2  1888.0007    final         120       1, 1     0.4    0.15
#> 136         2  2082.5300    final         120       1, 1     0.4    0.15
#> 137         2  2315.2835    final         120       1, 1     0.4    0.15
#> 138         2  2372.9135    final         120       1, 1     0.4    0.15
#> 139         2  2411.5581    final         120       1, 1     0.4    0.15
#> 140         2  2457.7624    final         120       1, 1     0.4    0.15
#> 141         2  2549.3934    final         120       1, 1     0.4    0.15
#> 142         2  2641.3126    final         120       1, 1     0.4    0.15
#> 143         2  2648.6126    final         120       1, 1     0.4    0.15
#> 144         2  2654.4026    final         120       1, 1     0.4    0.15
#> 145         2  2682.7957    final         120       1, 1     0.4    0.15
#> 146         2  2687.8458    final         120       1, 1     0.4    0.15
#> 147         2  2822.3858    final         120       1, 1     0.4    0.15
#> 148         2  2864.9707    final         120       1, 1     0.4    0.15
#> 149         2  2988.4824    final         120       1, 1     0.4    0.15
#> 150         2  3013.1278    final         120       1, 1     0.4    0.15
#> 151         2  3073.4366    final         120       1, 1     0.4    0.15
#> 152         2  3493.2524    final         120       1, 1     0.4    0.15
#> 153         2  3819.6447    final         120       1, 1     0.4    0.15
#> 154         2  3840.1027    final         120       1, 1     0.4    0.15
#> 155         2  4303.5245    final         120       1, 1     0.4    0.15
#> 156         2  4557.6462    final         120       1, 1     0.4    0.15
#> 157         2  4804.9065    final         120       1, 1     0.4    0.15
#> 158         2  4981.6575    final         120       1, 1     0.4    0.15
#> 159         2  5012.7323    final         120       1, 1     0.4    0.15
#> 160         2  5145.9628    final         120       1, 1     0.4    0.15
#> 161         2  5206.0924    final         120       1, 1     0.4    0.15
#> 162         2  5457.0055    final         120       1, 1     0.4    0.15
#> 163         2  5638.0848    final         120       1, 1     0.4    0.15
#> 164         2  5796.7006    final         120       1, 1     0.4    0.15
#> 165         2  6067.1212    final         120       1, 1     0.4    0.15
#> 166         2  6115.2388    final         120       1, 1     0.4    0.15
#> 167         2  6353.3218    final         120       1, 1     0.4    0.15
#> 168         2  6448.7033    final         120       1, 1     0.4    0.15
#> 169         2  6547.5079    final         120       1, 1     0.4    0.15
#> 170         2  6689.2264    final         120       1, 1     0.4    0.15
#> 171         2  6742.7127    final         120       1, 1     0.4    0.15
#> 172         2  7059.9069    final         120       1, 1     0.4    0.15
#> 173         2  7119.8497    final         120       1, 1     0.4    0.15
#> 174         2  7201.3364    final         120       1, 1     0.4    0.15
#> 175         2  7218.6617    final         120       1, 1     0.4    0.15
#> 176         2  7312.9388    final         120       1, 1     0.4    0.15
#> 177         2  7340.3264    final         120       1, 1     0.4    0.15
#> 178         2  7389.3001    final         120       1, 1     0.4    0.15
#> 179         2  7441.8725    final         120       1, 1     0.4    0.15
#> 180         2  7589.3324    final         120       1, 1     0.4    0.15
#> 181         2  7801.4383    final         120       1, 1     0.4    0.15
#> 182         2  7884.9057    final         120       1, 1     0.4    0.15
#> 183         2  8318.1828    final         120       1, 1     0.4    0.15
#> 184         2  8417.5600    final         120       1, 1     0.4    0.15
#> 185         2  8418.6991    final         120       1, 1     0.4    0.15
#> 186         2  8514.5156    final         120       1, 1     0.4    0.15
#> 187         2  8574.4622    final         120       1, 1     0.4    0.15
#> 188         2  8716.9881    final         120       1, 1     0.4    0.15
#> 189         2  8719.5517    final         120       1, 1     0.4    0.15
#> 190         2  8806.2886    final         120       1, 1     0.4    0.15
#> 191         2  8844.4981    final         120       1, 1     0.4    0.15
#> 192         2  8908.4505    final         120       1, 1     0.4    0.15
#> 193         2  9096.9945    final         120       1, 1     0.4    0.15
#> 194         2  9119.8056    final         120       1, 1     0.4    0.15
#> 195         2  9189.0470    final         120       1, 1     0.4    0.15
#> 196         2  9398.8525    final         120       1, 1     0.4    0.15
#> 197         2  9548.7794    final         120       1, 1     0.4    0.15
#> 198         2  9645.6699    final         120       1, 1     0.4    0.15
#> 199         2  9701.7895    final         120       1, 1     0.4    0.15
#> 200         2  9702.4022    final         120       1, 1     0.4    0.15
#> 201         2  9704.5875    final         120       1, 1     0.4    0.15
#> 202         2  9721.1057    final         120       1, 1     0.4    0.15
#> 203         2  9791.4548    final         120       1, 1     0.4    0.15
#> 204         2  9865.2968    final         120       1, 1     0.4    0.15
#> 205         2  9974.9392    final         120       1, 1     0.4    0.15
#> 206         2 10161.4662    final         120       1, 1     0.4    0.15
#> 207         2 10407.8681    final         120       1, 1     0.4    0.15
#> 208         2 10452.5790    final         120       1, 1     0.4    0.15
#> 209         2 10475.8074    final         120       1, 1     0.4    0.15
#> 210         2 10626.0228    final         120       1, 1     0.4    0.15
#> 211         2 10664.2511    final         120       1, 1     0.4    0.15
#> 212         2 10696.2755    final         120       1, 1     0.4    0.15
#> 213         2 10798.0761    final         120       1, 1     0.4    0.15
#> 214         2 10802.6704    final         120       1, 1     0.4    0.15
#> 215         2 10855.2878    final         120       1, 1     0.4    0.15
#> 216         2 10877.7262    final         120       1, 1     0.4    0.15
#> 217         2 10896.5621    final         120       1, 1     0.4    0.15
#> 218         2 11008.2027    final         120       1, 1     0.4    0.15
#> 219         2 11034.3357    final         120       1, 1     0.4    0.15
#> 220         2 11062.8649    final         120       1, 1     0.4    0.15
#> 221         2 11141.9174    final         120       1, 1     0.4    0.15
#> 222         2 11242.4553    final         120       1, 1     0.4    0.15
#> 223         2 11481.2195    final         120       1, 1     0.4    0.15
#> 224         2 11547.8767    final         120       1, 1     0.4    0.15
#> 225         2 11555.8539    final         120       1, 1     0.4    0.15
#> 226         2 11605.5359    final         120       1, 1     0.4    0.15
#> 227         2 11777.7078    final         120       1, 1     0.4    0.15
#> 228         2 11917.8879    final         120       1, 1     0.4    0.15
#> 229         2 12042.2040    final         120       1, 1     0.4    0.15
#> 230         2 12269.5897    final         120       1, 1     0.4    0.15
#> 231         2 12506.0034    final         120       1, 1     0.4    0.15
#> 232         2 12967.8777    final         120       1, 1     0.4    0.15
#> 233         2 13496.4057    final         120       1, 1     0.4    0.15
#> 234         2 13547.5968    final         120       1, 1     0.4    0.15
#> 235         2 13662.9178    final         120       1, 1     0.4    0.15
#> 236         2 13851.9485    final         120       1, 1     0.4    0.15
#> 237         2 14399.7511    final         120       1, 1     0.4    0.15
#> 238         2 14422.1896    final         120       1, 1     0.4    0.15
#> 239         2 14497.2541    final         120       1, 1     0.4    0.15
#> 240         3   117.7798    final         120       1, 1     0.4    0.15
#> 241         3   220.2221    final         120       1, 1     0.4    0.15
#> 242         3   431.2903    final         120       1, 1     0.4    0.15
#> 243         3   475.0459    final         120       1, 1     0.4    0.15
#> 244         3   551.2730    final         120       1, 1     0.4    0.15
#> 245         3   774.7013    final         120       1, 1     0.4    0.15
#> 246         3  1147.7307    final         120       1, 1     0.4    0.15
#> 247         3  1216.7620    final         120       1, 1     0.4    0.15
#> 248         3  1360.7081    final         120       1, 1     0.4    0.15
#> 249         3  1380.7569    final         120       1, 1     0.4    0.15
#> 250         3  1493.9849    final         120       1, 1     0.4    0.15
#> 251         3  1552.5583    final         120       1, 1     0.4    0.15
#> 252         3  1593.3474    final         120       1, 1     0.4    0.15
#> 253         3  1598.6306    final         120       1, 1     0.4    0.15
#> 254         3  1841.7779    final         120       1, 1     0.4    0.15
#> 255         3  2229.7360    final         120       1, 1     0.4    0.15
#> 256         3  2339.6198    final         120       1, 1     0.4    0.15
#> 257         3  2344.2087    final         120       1, 1     0.4    0.15
#> 258         3  2556.9493    final         120       1, 1     0.4    0.15
#> 259         3  2570.7601    final         120       1, 1     0.4    0.15
#> 260         3  2621.2186    final         120       1, 1     0.4    0.15
#> 261         3  2659.3525    final         120       1, 1     0.4    0.15
#> 262         3  2771.4882    final         120       1, 1     0.4    0.15
#> 263         3  2842.2807    final         120       1, 1     0.4    0.15
#> 264         3  2954.9869    final         120       1, 1     0.4    0.15
#> 265         3  3143.7629    final         120       1, 1     0.4    0.15
#> 266         3  3267.4622    final         120       1, 1     0.4    0.15
#> 267         3  3296.2668    final         120       1, 1     0.4    0.15
#> 268         3  3414.3779    final         120       1, 1     0.4    0.15
#> 269         3  3587.2647    final         120       1, 1     0.4    0.15
#> 270         3  3588.9781    final         120       1, 1     0.4    0.15
#> 271         3  3669.6895    final         120       1, 1     0.4    0.15
#> 272         3  3696.0434    final         120       1, 1     0.4    0.15
#> 273         3  3799.7831    final         120       1, 1     0.4    0.15
#> 274         3  3825.5802    final         120       1, 1     0.4    0.15
#> 275         3  3867.1478    final         120       1, 1     0.4    0.15
#> 276         3  4113.2838    final         120       1, 1     0.4    0.15
#> 277         3  4249.1165    final         120       1, 1     0.4    0.15
#> 278         3  4300.4661    final         120       1, 1     0.4    0.15
#> 279         3  4302.3988    final         120       1, 1     0.4    0.15
#> 280         3  4331.6992    final         120       1, 1     0.4    0.15
#> 281         3  4338.6772    final         120       1, 1     0.4    0.15
#> 282         3  4398.5320    final         120       1, 1     0.4    0.15
#> 283         3  4438.8326    final         120       1, 1     0.4    0.15
#> 284         3  4469.1345    final         120       1, 1     0.4    0.15
#> 285         3  4525.0403    final         120       1, 1     0.4    0.15
#> 286         3  4543.7775    final         120       1, 1     0.4    0.15
#> 287         3  4592.2930    final         120       1, 1     0.4    0.15
#> 288         3  4641.7373    final         120       1, 1     0.4    0.15
#> 289         3  4709.8040    final         120       1, 1     0.4    0.15
#> 290         3  4876.4647    final         120       1, 1     0.4    0.15
#> 291         3  4918.2614    final         120       1, 1     0.4    0.15
#> 292         3  5163.3083    final         120       1, 1     0.4    0.15
#> 293         3  5232.2349    final         120       1, 1     0.4    0.15
#> 294         3  5338.9924    final         120       1, 1     0.4    0.15
#> 295         3  5381.1980    final         120       1, 1     0.4    0.15
#> 296         3  5432.2330    final         120       1, 1     0.4    0.15
#> 297         3  5725.6024    final         120       1, 1     0.4    0.15
#> 298         3  5790.0258    final         120       1, 1     0.4    0.15
#> 299         3  5924.0980    final         120       1, 1     0.4    0.15
#> 300         3  6122.4618    final         120       1, 1     0.4    0.15
#> 301         3  6174.8339    final         120       1, 1     0.4    0.15
#> 302         3  6413.3257    final         120       1, 1     0.4    0.15
#> 303         3  6659.2488    final         120       1, 1     0.4    0.15
#> 304         3  6673.3416    final         120       1, 1     0.4    0.15
#> 305         3  6696.6901    final         120       1, 1     0.4    0.15
#> 306         3  7023.7990    final         120       1, 1     0.4    0.15
#> 307         3  7235.0918    final         120       1, 1     0.4    0.15
#> 308         3  7314.5256    final         120       1, 1     0.4    0.15
#> 309         3  7426.7892    final         120       1, 1     0.4    0.15
#> 310         3  7442.8183    final         120       1, 1     0.4    0.15
#> 311         3  7541.5965    final         120       1, 1     0.4    0.15
#> 312         3  7704.4783    final         120       1, 1     0.4    0.15
#> 313         3  7816.1523    final         120       1, 1     0.4    0.15
#> 314         3  7895.4647    final         120       1, 1     0.4    0.15
#> 315         3  7940.8404    final         120       1, 1     0.4    0.15
#> 316         3  8034.1345    final         120       1, 1     0.4    0.15
#> 317         3  8070.9385    final         120       1, 1     0.4    0.15
#> 318         3  8083.8012    final         120       1, 1     0.4    0.15
#> 319         3  8275.9989    final         120       1, 1     0.4    0.15
#> 320         3  8399.6192    final         120       1, 1     0.4    0.15
#> 321         3  8429.7178    final         120       1, 1     0.4    0.15
#> 322         3  8447.3282    final         120       1, 1     0.4    0.15
#> 323         3  8578.3312    final         120       1, 1     0.4    0.15
#> 324         3  8598.7597    final         120       1, 1     0.4    0.15
#> 325         3  8809.2365    final         120       1, 1     0.4    0.15
#> 326         3  8822.8053    final         120       1, 1     0.4    0.15
#> 327         3  9068.5278    final         120       1, 1     0.4    0.15
#> 328         3  9188.8585    final         120       1, 1     0.4    0.15
#> 329         3  9484.1493    final         120       1, 1     0.4    0.15
#> 330         3  9621.6985    final         120       1, 1     0.4    0.15
#> 331         3  9667.3301    final         120       1, 1     0.4    0.15
#> 332         3  9674.9429    final         120       1, 1     0.4    0.15
#> 333         3  9774.6793    final         120       1, 1     0.4    0.15
#> 334         3  9799.3427    final         120       1, 1     0.4    0.15
#> 335         3  9809.9191    final         120       1, 1     0.4    0.15
#> 336         3  9891.6097    final         120       1, 1     0.4    0.15
#> 337         3 10264.6379    final         120       1, 1     0.4    0.15
#> 338         3 10276.3116    final         120       1, 1     0.4    0.15
#> 339         3 10308.4945    final         120       1, 1     0.4    0.15
#> 340         3 10627.3193    final         120       1, 1     0.4    0.15
#> 341         3 10771.4948    final         120       1, 1     0.4    0.15
#> 342         3 10846.0150    final         120       1, 1     0.4    0.15
#> 343         3 10863.7261    final         120       1, 1     0.4    0.15
#> 344         3 11006.2680    final         120       1, 1     0.4    0.15
#> 345         3 11105.1478    final         120       1, 1     0.4    0.15
#> 346         3 11112.0992    final         120       1, 1     0.4    0.15
#> 347         3 11167.9745    final         120       1, 1     0.4    0.15
#> 348         3 11304.5211    final         120       1, 1     0.4    0.15
#> 349         3 11335.2773    final         120       1, 1     0.4    0.15
#> 350         3 11399.6479    final         120       1, 1     0.4    0.15
#> 351         3 11530.6670    final         120       1, 1     0.4    0.15
#> 352         3 11576.6874    final         120       1, 1     0.4    0.15
#> 353         3 11602.1352    final         120       1, 1     0.4    0.15
#> 354         3 11738.0804    final         120       1, 1     0.4    0.15
#> 355         3 11821.5127    final         120       1, 1     0.4    0.15
#> 356         3 11900.9172    final         120       1, 1     0.4    0.15
#> 357         3 12145.8338    final         120       1, 1     0.4    0.15
#> 358         3 12203.8313    final         120       1, 1     0.4    0.15
#> 359         3 12378.8956    final         120       1, 1     0.4    0.15
#> 360         3 12393.7843    final         120       1, 1     0.4    0.15
#>     p_subgroup_A p_subgroup_B n_total overall_estimate overall_p_value
#> 1            0.6          0.4     120        0.3414139      0.06906957
#> 2            0.6          0.4     120        0.3414139      0.06906957
#> 3            0.6          0.4     120        0.3414139      0.06906957
#> 4            0.6          0.4     120        0.3414139      0.06906957
#> 5            0.6          0.4     120        0.3414139      0.06906957
#> 6            0.6          0.4     120        0.3414139      0.06906957
#> 7            0.6          0.4     120        0.3414139      0.06906957
#> 8            0.6          0.4     120        0.3414139      0.06906957
#> 9            0.6          0.4     120        0.3414139      0.06906957
#> 10           0.6          0.4     120        0.3414139      0.06906957
#> 11           0.6          0.4     120        0.3414139      0.06906957
#> 12           0.6          0.4     120        0.3414139      0.06906957
#> 13           0.6          0.4     120        0.3414139      0.06906957
#> 14           0.6          0.4     120        0.3414139      0.06906957
#> 15           0.6          0.4     120        0.3414139      0.06906957
#> 16           0.6          0.4     120        0.3414139      0.06906957
#> 17           0.6          0.4     120        0.3414139      0.06906957
#> 18           0.6          0.4     120        0.3414139      0.06906957
#> 19           0.6          0.4     120        0.3414139      0.06906957
#> 20           0.6          0.4     120        0.3414139      0.06906957
#> 21           0.6          0.4     120        0.3414139      0.06906957
#> 22           0.6          0.4     120        0.3414139      0.06906957
#> 23           0.6          0.4     120        0.3414139      0.06906957
#> 24           0.6          0.4     120        0.3414139      0.06906957
#> 25           0.6          0.4     120        0.3414139      0.06906957
#> 26           0.6          0.4     120        0.3414139      0.06906957
#> 27           0.6          0.4     120        0.3414139      0.06906957
#> 28           0.6          0.4     120        0.3414139      0.06906957
#> 29           0.6          0.4     120        0.3414139      0.06906957
#> 30           0.6          0.4     120        0.3414139      0.06906957
#> 31           0.6          0.4     120        0.3414139      0.06906957
#> 32           0.6          0.4     120        0.3414139      0.06906957
#> 33           0.6          0.4     120        0.3414139      0.06906957
#> 34           0.6          0.4     120        0.3414139      0.06906957
#> 35           0.6          0.4     120        0.3414139      0.06906957
#> 36           0.6          0.4     120        0.3414139      0.06906957
#> 37           0.6          0.4     120        0.3414139      0.06906957
#> 38           0.6          0.4     120        0.3414139      0.06906957
#> 39           0.6          0.4     120        0.3414139      0.06906957
#> 40           0.6          0.4     120        0.3414139      0.06906957
#> 41           0.6          0.4     120        0.3414139      0.06906957
#> 42           0.6          0.4     120        0.3414139      0.06906957
#> 43           0.6          0.4     120        0.3414139      0.06906957
#> 44           0.6          0.4     120        0.3414139      0.06906957
#> 45           0.6          0.4     120        0.3414139      0.06906957
#> 46           0.6          0.4     120        0.3414139      0.06906957
#> 47           0.6          0.4     120        0.3414139      0.06906957
#> 48           0.6          0.4     120        0.3414139      0.06906957
#> 49           0.6          0.4     120        0.3414139      0.06906957
#> 50           0.6          0.4     120        0.3414139      0.06906957
#> 51           0.6          0.4     120        0.3414139      0.06906957
#> 52           0.6          0.4     120        0.3414139      0.06906957
#> 53           0.6          0.4     120        0.3414139      0.06906957
#> 54           0.6          0.4     120        0.3414139      0.06906957
#> 55           0.6          0.4     120        0.3414139      0.06906957
#> 56           0.6          0.4     120        0.3414139      0.06906957
#> 57           0.6          0.4     120        0.3414139      0.06906957
#> 58           0.6          0.4     120        0.3414139      0.06906957
#> 59           0.6          0.4     120        0.3414139      0.06906957
#> 60           0.6          0.4     120        0.3414139      0.06906957
#> 61           0.6          0.4     120        0.3414139      0.06906957
#> 62           0.6          0.4     120        0.3414139      0.06906957
#> 63           0.6          0.4     120        0.3414139      0.06906957
#> 64           0.6          0.4     120        0.3414139      0.06906957
#> 65           0.6          0.4     120        0.3414139      0.06906957
#> 66           0.6          0.4     120        0.3414139      0.06906957
#> 67           0.6          0.4     120        0.3414139      0.06906957
#> 68           0.6          0.4     120        0.3414139      0.06906957
#> 69           0.6          0.4     120        0.3414139      0.06906957
#> 70           0.6          0.4     120        0.3414139      0.06906957
#> 71           0.6          0.4     120        0.3414139      0.06906957
#> 72           0.6          0.4     120        0.3414139      0.06906957
#> 73           0.6          0.4     120        0.3414139      0.06906957
#> 74           0.6          0.4     120        0.3414139      0.06906957
#> 75           0.6          0.4     120        0.3414139      0.06906957
#> 76           0.6          0.4     120        0.3414139      0.06906957
#> 77           0.6          0.4     120        0.3414139      0.06906957
#> 78           0.6          0.4     120        0.3414139      0.06906957
#> 79           0.6          0.4     120        0.3414139      0.06906957
#> 80           0.6          0.4     120        0.3414139      0.06906957
#> 81           0.6          0.4     120        0.3414139      0.06906957
#> 82           0.6          0.4     120        0.3414139      0.06906957
#> 83           0.6          0.4     120        0.3414139      0.06906957
#> 84           0.6          0.4     120        0.3414139      0.06906957
#> 85           0.6          0.4     120        0.3414139      0.06906957
#> 86           0.6          0.4     120        0.3414139      0.06906957
#> 87           0.6          0.4     120        0.3414139      0.06906957
#> 88           0.6          0.4     120        0.3414139      0.06906957
#> 89           0.6          0.4     120        0.3414139      0.06906957
#> 90           0.6          0.4     120        0.3414139      0.06906957
#> 91           0.6          0.4     120        0.3414139      0.06906957
#> 92           0.6          0.4     120        0.3414139      0.06906957
#> 93           0.6          0.4     120        0.3414139      0.06906957
#> 94           0.6          0.4     120        0.3414139      0.06906957
#> 95           0.6          0.4     120        0.3414139      0.06906957
#> 96           0.6          0.4     120        0.3414139      0.06906957
#> 97           0.6          0.4     120        0.3414139      0.06906957
#> 98           0.6          0.4     120        0.3414139      0.06906957
#> 99           0.6          0.4     120        0.3414139      0.06906957
#> 100          0.6          0.4     120        0.3414139      0.06906957
#> 101          0.6          0.4     120        0.3414139      0.06906957
#> 102          0.6          0.4     120        0.3414139      0.06906957
#> 103          0.6          0.4     120        0.3414139      0.06906957
#> 104          0.6          0.4     120        0.3414139      0.06906957
#> 105          0.6          0.4     120        0.3414139      0.06906957
#> 106          0.6          0.4     120        0.3414139      0.06906957
#> 107          0.6          0.4     120        0.3414139      0.06906957
#> 108          0.6          0.4     120        0.3414139      0.06906957
#> 109          0.6          0.4     120        0.3414139      0.06906957
#> 110          0.6          0.4     120        0.3414139      0.06906957
#> 111          0.6          0.4     120        0.3414139      0.06906957
#> 112          0.6          0.4     120        0.3414139      0.06906957
#> 113          0.6          0.4     120        0.3414139      0.06906957
#> 114          0.6          0.4     120        0.3414139      0.06906957
#> 115          0.6          0.4     120        0.3414139      0.06906957
#> 116          0.6          0.4     120        0.3414139      0.06906957
#> 117          0.6          0.4     120        0.3414139      0.06906957
#> 118          0.6          0.4     120        0.3414139      0.06906957
#> 119          0.6          0.4     120        0.3414139      0.06906957
#> 120          0.6          0.4     120        0.1649366      0.37184435
#> 121          0.6          0.4     120        0.1649366      0.37184435
#> 122          0.6          0.4     120        0.1649366      0.37184435
#> 123          0.6          0.4     120        0.1649366      0.37184435
#> 124          0.6          0.4     120        0.1649366      0.37184435
#> 125          0.6          0.4     120        0.1649366      0.37184435
#> 126          0.6          0.4     120        0.1649366      0.37184435
#> 127          0.6          0.4     120        0.1649366      0.37184435
#> 128          0.6          0.4     120        0.1649366      0.37184435
#> 129          0.6          0.4     120        0.1649366      0.37184435
#> 130          0.6          0.4     120        0.1649366      0.37184435
#> 131          0.6          0.4     120        0.1649366      0.37184435
#> 132          0.6          0.4     120        0.1649366      0.37184435
#> 133          0.6          0.4     120        0.1649366      0.37184435
#> 134          0.6          0.4     120        0.1649366      0.37184435
#> 135          0.6          0.4     120        0.1649366      0.37184435
#> 136          0.6          0.4     120        0.1649366      0.37184435
#> 137          0.6          0.4     120        0.1649366      0.37184435
#> 138          0.6          0.4     120        0.1649366      0.37184435
#> 139          0.6          0.4     120        0.1649366      0.37184435
#> 140          0.6          0.4     120        0.1649366      0.37184435
#> 141          0.6          0.4     120        0.1649366      0.37184435
#> 142          0.6          0.4     120        0.1649366      0.37184435
#> 143          0.6          0.4     120        0.1649366      0.37184435
#> 144          0.6          0.4     120        0.1649366      0.37184435
#> 145          0.6          0.4     120        0.1649366      0.37184435
#> 146          0.6          0.4     120        0.1649366      0.37184435
#> 147          0.6          0.4     120        0.1649366      0.37184435
#> 148          0.6          0.4     120        0.1649366      0.37184435
#> 149          0.6          0.4     120        0.1649366      0.37184435
#> 150          0.6          0.4     120        0.1649366      0.37184435
#> 151          0.6          0.4     120        0.1649366      0.37184435
#> 152          0.6          0.4     120        0.1649366      0.37184435
#> 153          0.6          0.4     120        0.1649366      0.37184435
#> 154          0.6          0.4     120        0.1649366      0.37184435
#> 155          0.6          0.4     120        0.1649366      0.37184435
#> 156          0.6          0.4     120        0.1649366      0.37184435
#> 157          0.6          0.4     120        0.1649366      0.37184435
#> 158          0.6          0.4     120        0.1649366      0.37184435
#> 159          0.6          0.4     120        0.1649366      0.37184435
#> 160          0.6          0.4     120        0.1649366      0.37184435
#> 161          0.6          0.4     120        0.1649366      0.37184435
#> 162          0.6          0.4     120        0.1649366      0.37184435
#> 163          0.6          0.4     120        0.1649366      0.37184435
#> 164          0.6          0.4     120        0.1649366      0.37184435
#> 165          0.6          0.4     120        0.1649366      0.37184435
#> 166          0.6          0.4     120        0.1649366      0.37184435
#> 167          0.6          0.4     120        0.1649366      0.37184435
#> 168          0.6          0.4     120        0.1649366      0.37184435
#> 169          0.6          0.4     120        0.1649366      0.37184435
#> 170          0.6          0.4     120        0.1649366      0.37184435
#> 171          0.6          0.4     120        0.1649366      0.37184435
#> 172          0.6          0.4     120        0.1649366      0.37184435
#> 173          0.6          0.4     120        0.1649366      0.37184435
#> 174          0.6          0.4     120        0.1649366      0.37184435
#> 175          0.6          0.4     120        0.1649366      0.37184435
#> 176          0.6          0.4     120        0.1649366      0.37184435
#> 177          0.6          0.4     120        0.1649366      0.37184435
#> 178          0.6          0.4     120        0.1649366      0.37184435
#> 179          0.6          0.4     120        0.1649366      0.37184435
#> 180          0.6          0.4     120        0.1649366      0.37184435
#> 181          0.6          0.4     120        0.1649366      0.37184435
#> 182          0.6          0.4     120        0.1649366      0.37184435
#> 183          0.6          0.4     120        0.1649366      0.37184435
#> 184          0.6          0.4     120        0.1649366      0.37184435
#> 185          0.6          0.4     120        0.1649366      0.37184435
#> 186          0.6          0.4     120        0.1649366      0.37184435
#> 187          0.6          0.4     120        0.1649366      0.37184435
#> 188          0.6          0.4     120        0.1649366      0.37184435
#> 189          0.6          0.4     120        0.1649366      0.37184435
#> 190          0.6          0.4     120        0.1649366      0.37184435
#> 191          0.6          0.4     120        0.1649366      0.37184435
#> 192          0.6          0.4     120        0.1649366      0.37184435
#> 193          0.6          0.4     120        0.1649366      0.37184435
#> 194          0.6          0.4     120        0.1649366      0.37184435
#> 195          0.6          0.4     120        0.1649366      0.37184435
#> 196          0.6          0.4     120        0.1649366      0.37184435
#> 197          0.6          0.4     120        0.1649366      0.37184435
#> 198          0.6          0.4     120        0.1649366      0.37184435
#> 199          0.6          0.4     120        0.1649366      0.37184435
#> 200          0.6          0.4     120        0.1649366      0.37184435
#> 201          0.6          0.4     120        0.1649366      0.37184435
#> 202          0.6          0.4     120        0.1649366      0.37184435
#> 203          0.6          0.4     120        0.1649366      0.37184435
#> 204          0.6          0.4     120        0.1649366      0.37184435
#> 205          0.6          0.4     120        0.1649366      0.37184435
#> 206          0.6          0.4     120        0.1649366      0.37184435
#> 207          0.6          0.4     120        0.1649366      0.37184435
#> 208          0.6          0.4     120        0.1649366      0.37184435
#> 209          0.6          0.4     120        0.1649366      0.37184435
#> 210          0.6          0.4     120        0.1649366      0.37184435
#> 211          0.6          0.4     120        0.1649366      0.37184435
#> 212          0.6          0.4     120        0.1649366      0.37184435
#> 213          0.6          0.4     120        0.1649366      0.37184435
#> 214          0.6          0.4     120        0.1649366      0.37184435
#> 215          0.6          0.4     120        0.1649366      0.37184435
#> 216          0.6          0.4     120        0.1649366      0.37184435
#> 217          0.6          0.4     120        0.1649366      0.37184435
#> 218          0.6          0.4     120        0.1649366      0.37184435
#> 219          0.6          0.4     120        0.1649366      0.37184435
#> 220          0.6          0.4     120        0.1649366      0.37184435
#> 221          0.6          0.4     120        0.1649366      0.37184435
#> 222          0.6          0.4     120        0.1649366      0.37184435
#> 223          0.6          0.4     120        0.1649366      0.37184435
#> 224          0.6          0.4     120        0.1649366      0.37184435
#> 225          0.6          0.4     120        0.1649366      0.37184435
#> 226          0.6          0.4     120        0.1649366      0.37184435
#> 227          0.6          0.4     120        0.1649366      0.37184435
#> 228          0.6          0.4     120        0.1649366      0.37184435
#> 229          0.6          0.4     120        0.1649366      0.37184435
#> 230          0.6          0.4     120        0.1649366      0.37184435
#> 231          0.6          0.4     120        0.1649366      0.37184435
#> 232          0.6          0.4     120        0.1649366      0.37184435
#> 233          0.6          0.4     120        0.1649366      0.37184435
#> 234          0.6          0.4     120        0.1649366      0.37184435
#> 235          0.6          0.4     120        0.1649366      0.37184435
#> 236          0.6          0.4     120        0.1649366      0.37184435
#> 237          0.6          0.4     120        0.1649366      0.37184435
#> 238          0.6          0.4     120        0.1649366      0.37184435
#> 239          0.6          0.4     120        0.1649366      0.37184435
#> 240          0.6          0.4     120        0.3223337      0.08392587
#> 241          0.6          0.4     120        0.3223337      0.08392587
#> 242          0.6          0.4     120        0.3223337      0.08392587
#> 243          0.6          0.4     120        0.3223337      0.08392587
#> 244          0.6          0.4     120        0.3223337      0.08392587
#> 245          0.6          0.4     120        0.3223337      0.08392587
#> 246          0.6          0.4     120        0.3223337      0.08392587
#> 247          0.6          0.4     120        0.3223337      0.08392587
#> 248          0.6          0.4     120        0.3223337      0.08392587
#> 249          0.6          0.4     120        0.3223337      0.08392587
#> 250          0.6          0.4     120        0.3223337      0.08392587
#> 251          0.6          0.4     120        0.3223337      0.08392587
#> 252          0.6          0.4     120        0.3223337      0.08392587
#> 253          0.6          0.4     120        0.3223337      0.08392587
#> 254          0.6          0.4     120        0.3223337      0.08392587
#> 255          0.6          0.4     120        0.3223337      0.08392587
#> 256          0.6          0.4     120        0.3223337      0.08392587
#> 257          0.6          0.4     120        0.3223337      0.08392587
#> 258          0.6          0.4     120        0.3223337      0.08392587
#> 259          0.6          0.4     120        0.3223337      0.08392587
#> 260          0.6          0.4     120        0.3223337      0.08392587
#> 261          0.6          0.4     120        0.3223337      0.08392587
#> 262          0.6          0.4     120        0.3223337      0.08392587
#> 263          0.6          0.4     120        0.3223337      0.08392587
#> 264          0.6          0.4     120        0.3223337      0.08392587
#> 265          0.6          0.4     120        0.3223337      0.08392587
#> 266          0.6          0.4     120        0.3223337      0.08392587
#> 267          0.6          0.4     120        0.3223337      0.08392587
#> 268          0.6          0.4     120        0.3223337      0.08392587
#> 269          0.6          0.4     120        0.3223337      0.08392587
#> 270          0.6          0.4     120        0.3223337      0.08392587
#> 271          0.6          0.4     120        0.3223337      0.08392587
#> 272          0.6          0.4     120        0.3223337      0.08392587
#> 273          0.6          0.4     120        0.3223337      0.08392587
#> 274          0.6          0.4     120        0.3223337      0.08392587
#> 275          0.6          0.4     120        0.3223337      0.08392587
#> 276          0.6          0.4     120        0.3223337      0.08392587
#> 277          0.6          0.4     120        0.3223337      0.08392587
#> 278          0.6          0.4     120        0.3223337      0.08392587
#> 279          0.6          0.4     120        0.3223337      0.08392587
#> 280          0.6          0.4     120        0.3223337      0.08392587
#> 281          0.6          0.4     120        0.3223337      0.08392587
#> 282          0.6          0.4     120        0.3223337      0.08392587
#> 283          0.6          0.4     120        0.3223337      0.08392587
#> 284          0.6          0.4     120        0.3223337      0.08392587
#> 285          0.6          0.4     120        0.3223337      0.08392587
#> 286          0.6          0.4     120        0.3223337      0.08392587
#> 287          0.6          0.4     120        0.3223337      0.08392587
#> 288          0.6          0.4     120        0.3223337      0.08392587
#> 289          0.6          0.4     120        0.3223337      0.08392587
#> 290          0.6          0.4     120        0.3223337      0.08392587
#> 291          0.6          0.4     120        0.3223337      0.08392587
#> 292          0.6          0.4     120        0.3223337      0.08392587
#> 293          0.6          0.4     120        0.3223337      0.08392587
#> 294          0.6          0.4     120        0.3223337      0.08392587
#> 295          0.6          0.4     120        0.3223337      0.08392587
#> 296          0.6          0.4     120        0.3223337      0.08392587
#> 297          0.6          0.4     120        0.3223337      0.08392587
#> 298          0.6          0.4     120        0.3223337      0.08392587
#> 299          0.6          0.4     120        0.3223337      0.08392587
#> 300          0.6          0.4     120        0.3223337      0.08392587
#> 301          0.6          0.4     120        0.3223337      0.08392587
#> 302          0.6          0.4     120        0.3223337      0.08392587
#> 303          0.6          0.4     120        0.3223337      0.08392587
#> 304          0.6          0.4     120        0.3223337      0.08392587
#> 305          0.6          0.4     120        0.3223337      0.08392587
#> 306          0.6          0.4     120        0.3223337      0.08392587
#> 307          0.6          0.4     120        0.3223337      0.08392587
#> 308          0.6          0.4     120        0.3223337      0.08392587
#> 309          0.6          0.4     120        0.3223337      0.08392587
#> 310          0.6          0.4     120        0.3223337      0.08392587
#> 311          0.6          0.4     120        0.3223337      0.08392587
#> 312          0.6          0.4     120        0.3223337      0.08392587
#> 313          0.6          0.4     120        0.3223337      0.08392587
#> 314          0.6          0.4     120        0.3223337      0.08392587
#> 315          0.6          0.4     120        0.3223337      0.08392587
#> 316          0.6          0.4     120        0.3223337      0.08392587
#> 317          0.6          0.4     120        0.3223337      0.08392587
#> 318          0.6          0.4     120        0.3223337      0.08392587
#> 319          0.6          0.4     120        0.3223337      0.08392587
#> 320          0.6          0.4     120        0.3223337      0.08392587
#> 321          0.6          0.4     120        0.3223337      0.08392587
#> 322          0.6          0.4     120        0.3223337      0.08392587
#> 323          0.6          0.4     120        0.3223337      0.08392587
#> 324          0.6          0.4     120        0.3223337      0.08392587
#> 325          0.6          0.4     120        0.3223337      0.08392587
#> 326          0.6          0.4     120        0.3223337      0.08392587
#> 327          0.6          0.4     120        0.3223337      0.08392587
#> 328          0.6          0.4     120        0.3223337      0.08392587
#> 329          0.6          0.4     120        0.3223337      0.08392587
#> 330          0.6          0.4     120        0.3223337      0.08392587
#> 331          0.6          0.4     120        0.3223337      0.08392587
#> 332          0.6          0.4     120        0.3223337      0.08392587
#> 333          0.6          0.4     120        0.3223337      0.08392587
#> 334          0.6          0.4     120        0.3223337      0.08392587
#> 335          0.6          0.4     120        0.3223337      0.08392587
#> 336          0.6          0.4     120        0.3223337      0.08392587
#> 337          0.6          0.4     120        0.3223337      0.08392587
#> 338          0.6          0.4     120        0.3223337      0.08392587
#> 339          0.6          0.4     120        0.3223337      0.08392587
#> 340          0.6          0.4     120        0.3223337      0.08392587
#> 341          0.6          0.4     120        0.3223337      0.08392587
#> 342          0.6          0.4     120        0.3223337      0.08392587
#> 343          0.6          0.4     120        0.3223337      0.08392587
#> 344          0.6          0.4     120        0.3223337      0.08392587
#> 345          0.6          0.4     120        0.3223337      0.08392587
#> 346          0.6          0.4     120        0.3223337      0.08392587
#> 347          0.6          0.4     120        0.3223337      0.08392587
#> 348          0.6          0.4     120        0.3223337      0.08392587
#> 349          0.6          0.4     120        0.3223337      0.08392587
#> 350          0.6          0.4     120        0.3223337      0.08392587
#> 351          0.6          0.4     120        0.3223337      0.08392587
#> 352          0.6          0.4     120        0.3223337      0.08392587
#> 353          0.6          0.4     120        0.3223337      0.08392587
#> 354          0.6          0.4     120        0.3223337      0.08392587
#> 355          0.6          0.4     120        0.3223337      0.08392587
#> 356          0.6          0.4     120        0.3223337      0.08392587
#> 357          0.6          0.4     120        0.3223337      0.08392587
#> 358          0.6          0.4     120        0.3223337      0.08392587
#> 359          0.6          0.4     120        0.3223337      0.08392587
#> 360          0.6          0.4     120        0.3223337      0.08392587
#>     subgroup_A_n subgroup_A_estimate subgroup_A_p_value subgroup_B_n
#> 1             76           0.4707117         0.05403434           44
#> 2             76           0.4707117         0.05403434           44
#> 3             76           0.4707117         0.05403434           44
#> 4             76           0.4707117         0.05403434           44
#> 5             76           0.4707117         0.05403434           44
#> 6             76           0.4707117         0.05403434           44
#> 7             76           0.4707117         0.05403434           44
#> 8             76           0.4707117         0.05403434           44
#> 9             76           0.4707117         0.05403434           44
#> 10            76           0.4707117         0.05403434           44
#> 11            76           0.4707117         0.05403434           44
#> 12            76           0.4707117         0.05403434           44
#> 13            76           0.4707117         0.05403434           44
#> 14            76           0.4707117         0.05403434           44
#> 15            76           0.4707117         0.05403434           44
#> 16            76           0.4707117         0.05403434           44
#> 17            76           0.4707117         0.05403434           44
#> 18            76           0.4707117         0.05403434           44
#> 19            76           0.4707117         0.05403434           44
#> 20            76           0.4707117         0.05403434           44
#> 21            76           0.4707117         0.05403434           44
#> 22            76           0.4707117         0.05403434           44
#> 23            76           0.4707117         0.05403434           44
#> 24            76           0.4707117         0.05403434           44
#> 25            76           0.4707117         0.05403434           44
#> 26            76           0.4707117         0.05403434           44
#> 27            76           0.4707117         0.05403434           44
#> 28            76           0.4707117         0.05403434           44
#> 29            76           0.4707117         0.05403434           44
#> 30            76           0.4707117         0.05403434           44
#> 31            76           0.4707117         0.05403434           44
#> 32            76           0.4707117         0.05403434           44
#> 33            76           0.4707117         0.05403434           44
#> 34            76           0.4707117         0.05403434           44
#> 35            76           0.4707117         0.05403434           44
#> 36            76           0.4707117         0.05403434           44
#> 37            76           0.4707117         0.05403434           44
#> 38            76           0.4707117         0.05403434           44
#> 39            76           0.4707117         0.05403434           44
#> 40            76           0.4707117         0.05403434           44
#> 41            76           0.4707117         0.05403434           44
#> 42            76           0.4707117         0.05403434           44
#> 43            76           0.4707117         0.05403434           44
#> 44            76           0.4707117         0.05403434           44
#> 45            76           0.4707117         0.05403434           44
#> 46            76           0.4707117         0.05403434           44
#> 47            76           0.4707117         0.05403434           44
#> 48            76           0.4707117         0.05403434           44
#> 49            76           0.4707117         0.05403434           44
#> 50            76           0.4707117         0.05403434           44
#> 51            76           0.4707117         0.05403434           44
#> 52            76           0.4707117         0.05403434           44
#> 53            76           0.4707117         0.05403434           44
#> 54            76           0.4707117         0.05403434           44
#> 55            76           0.4707117         0.05403434           44
#> 56            76           0.4707117         0.05403434           44
#> 57            76           0.4707117         0.05403434           44
#> 58            76           0.4707117         0.05403434           44
#> 59            76           0.4707117         0.05403434           44
#> 60            76           0.4707117         0.05403434           44
#> 61            76           0.4707117         0.05403434           44
#> 62            76           0.4707117         0.05403434           44
#> 63            76           0.4707117         0.05403434           44
#> 64            76           0.4707117         0.05403434           44
#> 65            76           0.4707117         0.05403434           44
#> 66            76           0.4707117         0.05403434           44
#> 67            76           0.4707117         0.05403434           44
#> 68            76           0.4707117         0.05403434           44
#> 69            76           0.4707117         0.05403434           44
#> 70            76           0.4707117         0.05403434           44
#> 71            76           0.4707117         0.05403434           44
#> 72            76           0.4707117         0.05403434           44
#> 73            76           0.4707117         0.05403434           44
#> 74            76           0.4707117         0.05403434           44
#> 75            76           0.4707117         0.05403434           44
#> 76            76           0.4707117         0.05403434           44
#> 77            76           0.4707117         0.05403434           44
#> 78            76           0.4707117         0.05403434           44
#> 79            76           0.4707117         0.05403434           44
#> 80            76           0.4707117         0.05403434           44
#> 81            76           0.4707117         0.05403434           44
#> 82            76           0.4707117         0.05403434           44
#> 83            76           0.4707117         0.05403434           44
#> 84            76           0.4707117         0.05403434           44
#> 85            76           0.4707117         0.05403434           44
#> 86            76           0.4707117         0.05403434           44
#> 87            76           0.4707117         0.05403434           44
#> 88            76           0.4707117         0.05403434           44
#> 89            76           0.4707117         0.05403434           44
#> 90            76           0.4707117         0.05403434           44
#> 91            76           0.4707117         0.05403434           44
#> 92            76           0.4707117         0.05403434           44
#> 93            76           0.4707117         0.05403434           44
#> 94            76           0.4707117         0.05403434           44
#> 95            76           0.4707117         0.05403434           44
#> 96            76           0.4707117         0.05403434           44
#> 97            76           0.4707117         0.05403434           44
#> 98            76           0.4707117         0.05403434           44
#> 99            76           0.4707117         0.05403434           44
#> 100           76           0.4707117         0.05403434           44
#> 101           76           0.4707117         0.05403434           44
#> 102           76           0.4707117         0.05403434           44
#> 103           76           0.4707117         0.05403434           44
#> 104           76           0.4707117         0.05403434           44
#> 105           76           0.4707117         0.05403434           44
#> 106           76           0.4707117         0.05403434           44
#> 107           76           0.4707117         0.05403434           44
#> 108           76           0.4707117         0.05403434           44
#> 109           76           0.4707117         0.05403434           44
#> 110           76           0.4707117         0.05403434           44
#> 111           76           0.4707117         0.05403434           44
#> 112           76           0.4707117         0.05403434           44
#> 113           76           0.4707117         0.05403434           44
#> 114           76           0.4707117         0.05403434           44
#> 115           76           0.4707117         0.05403434           44
#> 116           76           0.4707117         0.05403434           44
#> 117           76           0.4707117         0.05403434           44
#> 118           76           0.4707117         0.05403434           44
#> 119           76           0.4707117         0.05403434           44
#> 120           64           0.3320609         0.18188722           56
#> 121           64           0.3320609         0.18188722           56
#> 122           64           0.3320609         0.18188722           56
#> 123           64           0.3320609         0.18188722           56
#> 124           64           0.3320609         0.18188722           56
#> 125           64           0.3320609         0.18188722           56
#> 126           64           0.3320609         0.18188722           56
#> 127           64           0.3320609         0.18188722           56
#> 128           64           0.3320609         0.18188722           56
#> 129           64           0.3320609         0.18188722           56
#> 130           64           0.3320609         0.18188722           56
#> 131           64           0.3320609         0.18188722           56
#> 132           64           0.3320609         0.18188722           56
#> 133           64           0.3320609         0.18188722           56
#> 134           64           0.3320609         0.18188722           56
#> 135           64           0.3320609         0.18188722           56
#> 136           64           0.3320609         0.18188722           56
#> 137           64           0.3320609         0.18188722           56
#> 138           64           0.3320609         0.18188722           56
#> 139           64           0.3320609         0.18188722           56
#> 140           64           0.3320609         0.18188722           56
#> 141           64           0.3320609         0.18188722           56
#> 142           64           0.3320609         0.18188722           56
#> 143           64           0.3320609         0.18188722           56
#> 144           64           0.3320609         0.18188722           56
#> 145           64           0.3320609         0.18188722           56
#> 146           64           0.3320609         0.18188722           56
#> 147           64           0.3320609         0.18188722           56
#> 148           64           0.3320609         0.18188722           56
#> 149           64           0.3320609         0.18188722           56
#> 150           64           0.3320609         0.18188722           56
#> 151           64           0.3320609         0.18188722           56
#> 152           64           0.3320609         0.18188722           56
#> 153           64           0.3320609         0.18188722           56
#> 154           64           0.3320609         0.18188722           56
#> 155           64           0.3320609         0.18188722           56
#> 156           64           0.3320609         0.18188722           56
#> 157           64           0.3320609         0.18188722           56
#> 158           64           0.3320609         0.18188722           56
#> 159           64           0.3320609         0.18188722           56
#> 160           64           0.3320609         0.18188722           56
#> 161           64           0.3320609         0.18188722           56
#> 162           64           0.3320609         0.18188722           56
#> 163           64           0.3320609         0.18188722           56
#> 164           64           0.3320609         0.18188722           56
#> 165           64           0.3320609         0.18188722           56
#> 166           64           0.3320609         0.18188722           56
#> 167           64           0.3320609         0.18188722           56
#> 168           64           0.3320609         0.18188722           56
#> 169           64           0.3320609         0.18188722           56
#> 170           64           0.3320609         0.18188722           56
#> 171           64           0.3320609         0.18188722           56
#> 172           64           0.3320609         0.18188722           56
#> 173           64           0.3320609         0.18188722           56
#> 174           64           0.3320609         0.18188722           56
#> 175           64           0.3320609         0.18188722           56
#> 176           64           0.3320609         0.18188722           56
#> 177           64           0.3320609         0.18188722           56
#> 178           64           0.3320609         0.18188722           56
#> 179           64           0.3320609         0.18188722           56
#> 180           64           0.3320609         0.18188722           56
#> 181           64           0.3320609         0.18188722           56
#> 182           64           0.3320609         0.18188722           56
#> 183           64           0.3320609         0.18188722           56
#> 184           64           0.3320609         0.18188722           56
#> 185           64           0.3320609         0.18188722           56
#> 186           64           0.3320609         0.18188722           56
#> 187           64           0.3320609         0.18188722           56
#> 188           64           0.3320609         0.18188722           56
#> 189           64           0.3320609         0.18188722           56
#> 190           64           0.3320609         0.18188722           56
#> 191           64           0.3320609         0.18188722           56
#> 192           64           0.3320609         0.18188722           56
#> 193           64           0.3320609         0.18188722           56
#> 194           64           0.3320609         0.18188722           56
#> 195           64           0.3320609         0.18188722           56
#> 196           64           0.3320609         0.18188722           56
#> 197           64           0.3320609         0.18188722           56
#> 198           64           0.3320609         0.18188722           56
#> 199           64           0.3320609         0.18188722           56
#> 200           64           0.3320609         0.18188722           56
#> 201           64           0.3320609         0.18188722           56
#> 202           64           0.3320609         0.18188722           56
#> 203           64           0.3320609         0.18188722           56
#> 204           64           0.3320609         0.18188722           56
#> 205           64           0.3320609         0.18188722           56
#> 206           64           0.3320609         0.18188722           56
#> 207           64           0.3320609         0.18188722           56
#> 208           64           0.3320609         0.18188722           56
#> 209           64           0.3320609         0.18188722           56
#> 210           64           0.3320609         0.18188722           56
#> 211           64           0.3320609         0.18188722           56
#> 212           64           0.3320609         0.18188722           56
#> 213           64           0.3320609         0.18188722           56
#> 214           64           0.3320609         0.18188722           56
#> 215           64           0.3320609         0.18188722           56
#> 216           64           0.3320609         0.18188722           56
#> 217           64           0.3320609         0.18188722           56
#> 218           64           0.3320609         0.18188722           56
#> 219           64           0.3320609         0.18188722           56
#> 220           64           0.3320609         0.18188722           56
#> 221           64           0.3320609         0.18188722           56
#> 222           64           0.3320609         0.18188722           56
#> 223           64           0.3320609         0.18188722           56
#> 224           64           0.3320609         0.18188722           56
#> 225           64           0.3320609         0.18188722           56
#> 226           64           0.3320609         0.18188722           56
#> 227           64           0.3320609         0.18188722           56
#> 228           64           0.3320609         0.18188722           56
#> 229           64           0.3320609         0.18188722           56
#> 230           64           0.3320609         0.18188722           56
#> 231           64           0.3320609         0.18188722           56
#> 232           64           0.3320609         0.18188722           56
#> 233           64           0.3320609         0.18188722           56
#> 234           64           0.3320609         0.18188722           56
#> 235           64           0.3320609         0.18188722           56
#> 236           64           0.3320609         0.18188722           56
#> 237           64           0.3320609         0.18188722           56
#> 238           64           0.3320609         0.18188722           56
#> 239           64           0.3320609         0.18188722           56
#> 240           71           0.1774698         0.50062690           49
#> 241           71           0.1774698         0.50062690           49
#> 242           71           0.1774698         0.50062690           49
#> 243           71           0.1774698         0.50062690           49
#> 244           71           0.1774698         0.50062690           49
#> 245           71           0.1774698         0.50062690           49
#> 246           71           0.1774698         0.50062690           49
#> 247           71           0.1774698         0.50062690           49
#> 248           71           0.1774698         0.50062690           49
#> 249           71           0.1774698         0.50062690           49
#> 250           71           0.1774698         0.50062690           49
#> 251           71           0.1774698         0.50062690           49
#> 252           71           0.1774698         0.50062690           49
#> 253           71           0.1774698         0.50062690           49
#> 254           71           0.1774698         0.50062690           49
#> 255           71           0.1774698         0.50062690           49
#> 256           71           0.1774698         0.50062690           49
#> 257           71           0.1774698         0.50062690           49
#> 258           71           0.1774698         0.50062690           49
#> 259           71           0.1774698         0.50062690           49
#> 260           71           0.1774698         0.50062690           49
#> 261           71           0.1774698         0.50062690           49
#> 262           71           0.1774698         0.50062690           49
#> 263           71           0.1774698         0.50062690           49
#> 264           71           0.1774698         0.50062690           49
#> 265           71           0.1774698         0.50062690           49
#> 266           71           0.1774698         0.50062690           49
#> 267           71           0.1774698         0.50062690           49
#> 268           71           0.1774698         0.50062690           49
#> 269           71           0.1774698         0.50062690           49
#> 270           71           0.1774698         0.50062690           49
#> 271           71           0.1774698         0.50062690           49
#> 272           71           0.1774698         0.50062690           49
#> 273           71           0.1774698         0.50062690           49
#> 274           71           0.1774698         0.50062690           49
#> 275           71           0.1774698         0.50062690           49
#> 276           71           0.1774698         0.50062690           49
#> 277           71           0.1774698         0.50062690           49
#> 278           71           0.1774698         0.50062690           49
#> 279           71           0.1774698         0.50062690           49
#> 280           71           0.1774698         0.50062690           49
#> 281           71           0.1774698         0.50062690           49
#> 282           71           0.1774698         0.50062690           49
#> 283           71           0.1774698         0.50062690           49
#> 284           71           0.1774698         0.50062690           49
#> 285           71           0.1774698         0.50062690           49
#> 286           71           0.1774698         0.50062690           49
#> 287           71           0.1774698         0.50062690           49
#> 288           71           0.1774698         0.50062690           49
#> 289           71           0.1774698         0.50062690           49
#> 290           71           0.1774698         0.50062690           49
#> 291           71           0.1774698         0.50062690           49
#> 292           71           0.1774698         0.50062690           49
#> 293           71           0.1774698         0.50062690           49
#> 294           71           0.1774698         0.50062690           49
#> 295           71           0.1774698         0.50062690           49
#> 296           71           0.1774698         0.50062690           49
#> 297           71           0.1774698         0.50062690           49
#> 298           71           0.1774698         0.50062690           49
#> 299           71           0.1774698         0.50062690           49
#> 300           71           0.1774698         0.50062690           49
#> 301           71           0.1774698         0.50062690           49
#> 302           71           0.1774698         0.50062690           49
#> 303           71           0.1774698         0.50062690           49
#> 304           71           0.1774698         0.50062690           49
#> 305           71           0.1774698         0.50062690           49
#> 306           71           0.1774698         0.50062690           49
#> 307           71           0.1774698         0.50062690           49
#> 308           71           0.1774698         0.50062690           49
#> 309           71           0.1774698         0.50062690           49
#> 310           71           0.1774698         0.50062690           49
#> 311           71           0.1774698         0.50062690           49
#> 312           71           0.1774698         0.50062690           49
#> 313           71           0.1774698         0.50062690           49
#> 314           71           0.1774698         0.50062690           49
#> 315           71           0.1774698         0.50062690           49
#> 316           71           0.1774698         0.50062690           49
#> 317           71           0.1774698         0.50062690           49
#> 318           71           0.1774698         0.50062690           49
#> 319           71           0.1774698         0.50062690           49
#> 320           71           0.1774698         0.50062690           49
#> 321           71           0.1774698         0.50062690           49
#> 322           71           0.1774698         0.50062690           49
#> 323           71           0.1774698         0.50062690           49
#> 324           71           0.1774698         0.50062690           49
#> 325           71           0.1774698         0.50062690           49
#> 326           71           0.1774698         0.50062690           49
#> 327           71           0.1774698         0.50062690           49
#> 328           71           0.1774698         0.50062690           49
#> 329           71           0.1774698         0.50062690           49
#> 330           71           0.1774698         0.50062690           49
#> 331           71           0.1774698         0.50062690           49
#> 332           71           0.1774698         0.50062690           49
#> 333           71           0.1774698         0.50062690           49
#> 334           71           0.1774698         0.50062690           49
#> 335           71           0.1774698         0.50062690           49
#> 336           71           0.1774698         0.50062690           49
#> 337           71           0.1774698         0.50062690           49
#> 338           71           0.1774698         0.50062690           49
#> 339           71           0.1774698         0.50062690           49
#> 340           71           0.1774698         0.50062690           49
#> 341           71           0.1774698         0.50062690           49
#> 342           71           0.1774698         0.50062690           49
#> 343           71           0.1774698         0.50062690           49
#> 344           71           0.1774698         0.50062690           49
#> 345           71           0.1774698         0.50062690           49
#> 346           71           0.1774698         0.50062690           49
#> 347           71           0.1774698         0.50062690           49
#> 348           71           0.1774698         0.50062690           49
#> 349           71           0.1774698         0.50062690           49
#> 350           71           0.1774698         0.50062690           49
#> 351           71           0.1774698         0.50062690           49
#> 352           71           0.1774698         0.50062690           49
#> 353           71           0.1774698         0.50062690           49
#> 354           71           0.1774698         0.50062690           49
#> 355           71           0.1774698         0.50062690           49
#> 356           71           0.1774698         0.50062690           49
#> 357           71           0.1774698         0.50062690           49
#> 358           71           0.1774698         0.50062690           49
#> 359           71           0.1774698         0.50062690           49
#> 360           71           0.1774698         0.50062690           49
#>     subgroup_B_estimate subgroup_B_p_value interaction_p_value
#> 1            0.11808122         0.68769806           0.3632263
#> 2            0.11808122         0.68769806           0.3632263
#> 3            0.11808122         0.68769806           0.3632263
#> 4            0.11808122         0.68769806           0.3632263
#> 5            0.11808122         0.68769806           0.3632263
#> 6            0.11808122         0.68769806           0.3632263
#> 7            0.11808122         0.68769806           0.3632263
#> 8            0.11808122         0.68769806           0.3632263
#> 9            0.11808122         0.68769806           0.3632263
#> 10           0.11808122         0.68769806           0.3632263
#> 11           0.11808122         0.68769806           0.3632263
#> 12           0.11808122         0.68769806           0.3632263
#> 13           0.11808122         0.68769806           0.3632263
#> 14           0.11808122         0.68769806           0.3632263
#> 15           0.11808122         0.68769806           0.3632263
#> 16           0.11808122         0.68769806           0.3632263
#> 17           0.11808122         0.68769806           0.3632263
#> 18           0.11808122         0.68769806           0.3632263
#> 19           0.11808122         0.68769806           0.3632263
#> 20           0.11808122         0.68769806           0.3632263
#> 21           0.11808122         0.68769806           0.3632263
#> 22           0.11808122         0.68769806           0.3632263
#> 23           0.11808122         0.68769806           0.3632263
#> 24           0.11808122         0.68769806           0.3632263
#> 25           0.11808122         0.68769806           0.3632263
#> 26           0.11808122         0.68769806           0.3632263
#> 27           0.11808122         0.68769806           0.3632263
#> 28           0.11808122         0.68769806           0.3632263
#> 29           0.11808122         0.68769806           0.3632263
#> 30           0.11808122         0.68769806           0.3632263
#> 31           0.11808122         0.68769806           0.3632263
#> 32           0.11808122         0.68769806           0.3632263
#> 33           0.11808122         0.68769806           0.3632263
#> 34           0.11808122         0.68769806           0.3632263
#> 35           0.11808122         0.68769806           0.3632263
#> 36           0.11808122         0.68769806           0.3632263
#> 37           0.11808122         0.68769806           0.3632263
#> 38           0.11808122         0.68769806           0.3632263
#> 39           0.11808122         0.68769806           0.3632263
#> 40           0.11808122         0.68769806           0.3632263
#> 41           0.11808122         0.68769806           0.3632263
#> 42           0.11808122         0.68769806           0.3632263
#> 43           0.11808122         0.68769806           0.3632263
#> 44           0.11808122         0.68769806           0.3632263
#> 45           0.11808122         0.68769806           0.3632263
#> 46           0.11808122         0.68769806           0.3632263
#> 47           0.11808122         0.68769806           0.3632263
#> 48           0.11808122         0.68769806           0.3632263
#> 49           0.11808122         0.68769806           0.3632263
#> 50           0.11808122         0.68769806           0.3632263
#> 51           0.11808122         0.68769806           0.3632263
#> 52           0.11808122         0.68769806           0.3632263
#> 53           0.11808122         0.68769806           0.3632263
#> 54           0.11808122         0.68769806           0.3632263
#> 55           0.11808122         0.68769806           0.3632263
#> 56           0.11808122         0.68769806           0.3632263
#> 57           0.11808122         0.68769806           0.3632263
#> 58           0.11808122         0.68769806           0.3632263
#> 59           0.11808122         0.68769806           0.3632263
#> 60           0.11808122         0.68769806           0.3632263
#> 61           0.11808122         0.68769806           0.3632263
#> 62           0.11808122         0.68769806           0.3632263
#> 63           0.11808122         0.68769806           0.3632263
#> 64           0.11808122         0.68769806           0.3632263
#> 65           0.11808122         0.68769806           0.3632263
#> 66           0.11808122         0.68769806           0.3632263
#> 67           0.11808122         0.68769806           0.3632263
#> 68           0.11808122         0.68769806           0.3632263
#> 69           0.11808122         0.68769806           0.3632263
#> 70           0.11808122         0.68769806           0.3632263
#> 71           0.11808122         0.68769806           0.3632263
#> 72           0.11808122         0.68769806           0.3632263
#> 73           0.11808122         0.68769806           0.3632263
#> 74           0.11808122         0.68769806           0.3632263
#> 75           0.11808122         0.68769806           0.3632263
#> 76           0.11808122         0.68769806           0.3632263
#> 77           0.11808122         0.68769806           0.3632263
#> 78           0.11808122         0.68769806           0.3632263
#> 79           0.11808122         0.68769806           0.3632263
#> 80           0.11808122         0.68769806           0.3632263
#> 81           0.11808122         0.68769806           0.3632263
#> 82           0.11808122         0.68769806           0.3632263
#> 83           0.11808122         0.68769806           0.3632263
#> 84           0.11808122         0.68769806           0.3632263
#> 85           0.11808122         0.68769806           0.3632263
#> 86           0.11808122         0.68769806           0.3632263
#> 87           0.11808122         0.68769806           0.3632263
#> 88           0.11808122         0.68769806           0.3632263
#> 89           0.11808122         0.68769806           0.3632263
#> 90           0.11808122         0.68769806           0.3632263
#> 91           0.11808122         0.68769806           0.3632263
#> 92           0.11808122         0.68769806           0.3632263
#> 93           0.11808122         0.68769806           0.3632263
#> 94           0.11808122         0.68769806           0.3632263
#> 95           0.11808122         0.68769806           0.3632263
#> 96           0.11808122         0.68769806           0.3632263
#> 97           0.11808122         0.68769806           0.3632263
#> 98           0.11808122         0.68769806           0.3632263
#> 99           0.11808122         0.68769806           0.3632263
#> 100          0.11808122         0.68769806           0.3632263
#> 101          0.11808122         0.68769806           0.3632263
#> 102          0.11808122         0.68769806           0.3632263
#> 103          0.11808122         0.68769806           0.3632263
#> 104          0.11808122         0.68769806           0.3632263
#> 105          0.11808122         0.68769806           0.3632263
#> 106          0.11808122         0.68769806           0.3632263
#> 107          0.11808122         0.68769806           0.3632263
#> 108          0.11808122         0.68769806           0.3632263
#> 109          0.11808122         0.68769806           0.3632263
#> 110          0.11808122         0.68769806           0.3632263
#> 111          0.11808122         0.68769806           0.3632263
#> 112          0.11808122         0.68769806           0.3632263
#> 113          0.11808122         0.68769806           0.3632263
#> 114          0.11808122         0.68769806           0.3632263
#> 115          0.11808122         0.68769806           0.3632263
#> 116          0.11808122         0.68769806           0.3632263
#> 117          0.11808122         0.68769806           0.3632263
#> 118          0.11808122         0.68769806           0.3632263
#> 119          0.11808122         0.68769806           0.3632263
#> 120         -0.01071538         0.96958098           0.3573194
#> 121         -0.01071538         0.96958098           0.3573194
#> 122         -0.01071538         0.96958098           0.3573194
#> 123         -0.01071538         0.96958098           0.3573194
#> 124         -0.01071538         0.96958098           0.3573194
#> 125         -0.01071538         0.96958098           0.3573194
#> 126         -0.01071538         0.96958098           0.3573194
#> 127         -0.01071538         0.96958098           0.3573194
#> 128         -0.01071538         0.96958098           0.3573194
#> 129         -0.01071538         0.96958098           0.3573194
#> 130         -0.01071538         0.96958098           0.3573194
#> 131         -0.01071538         0.96958098           0.3573194
#> 132         -0.01071538         0.96958098           0.3573194
#> 133         -0.01071538         0.96958098           0.3573194
#> 134         -0.01071538         0.96958098           0.3573194
#> 135         -0.01071538         0.96958098           0.3573194
#> 136         -0.01071538         0.96958098           0.3573194
#> 137         -0.01071538         0.96958098           0.3573194
#> 138         -0.01071538         0.96958098           0.3573194
#> 139         -0.01071538         0.96958098           0.3573194
#> 140         -0.01071538         0.96958098           0.3573194
#> 141         -0.01071538         0.96958098           0.3573194
#> 142         -0.01071538         0.96958098           0.3573194
#> 143         -0.01071538         0.96958098           0.3573194
#> 144         -0.01071538         0.96958098           0.3573194
#> 145         -0.01071538         0.96958098           0.3573194
#> 146         -0.01071538         0.96958098           0.3573194
#> 147         -0.01071538         0.96958098           0.3573194
#> 148         -0.01071538         0.96958098           0.3573194
#> 149         -0.01071538         0.96958098           0.3573194
#> 150         -0.01071538         0.96958098           0.3573194
#> 151         -0.01071538         0.96958098           0.3573194
#> 152         -0.01071538         0.96958098           0.3573194
#> 153         -0.01071538         0.96958098           0.3573194
#> 154         -0.01071538         0.96958098           0.3573194
#> 155         -0.01071538         0.96958098           0.3573194
#> 156         -0.01071538         0.96958098           0.3573194
#> 157         -0.01071538         0.96958098           0.3573194
#> 158         -0.01071538         0.96958098           0.3573194
#> 159         -0.01071538         0.96958098           0.3573194
#> 160         -0.01071538         0.96958098           0.3573194
#> 161         -0.01071538         0.96958098           0.3573194
#> 162         -0.01071538         0.96958098           0.3573194
#> 163         -0.01071538         0.96958098           0.3573194
#> 164         -0.01071538         0.96958098           0.3573194
#> 165         -0.01071538         0.96958098           0.3573194
#> 166         -0.01071538         0.96958098           0.3573194
#> 167         -0.01071538         0.96958098           0.3573194
#> 168         -0.01071538         0.96958098           0.3573194
#> 169         -0.01071538         0.96958098           0.3573194
#> 170         -0.01071538         0.96958098           0.3573194
#> 171         -0.01071538         0.96958098           0.3573194
#> 172         -0.01071538         0.96958098           0.3573194
#> 173         -0.01071538         0.96958098           0.3573194
#> 174         -0.01071538         0.96958098           0.3573194
#> 175         -0.01071538         0.96958098           0.3573194
#> 176         -0.01071538         0.96958098           0.3573194
#> 177         -0.01071538         0.96958098           0.3573194
#> 178         -0.01071538         0.96958098           0.3573194
#> 179         -0.01071538         0.96958098           0.3573194
#> 180         -0.01071538         0.96958098           0.3573194
#> 181         -0.01071538         0.96958098           0.3573194
#> 182         -0.01071538         0.96958098           0.3573194
#> 183         -0.01071538         0.96958098           0.3573194
#> 184         -0.01071538         0.96958098           0.3573194
#> 185         -0.01071538         0.96958098           0.3573194
#> 186         -0.01071538         0.96958098           0.3573194
#> 187         -0.01071538         0.96958098           0.3573194
#> 188         -0.01071538         0.96958098           0.3573194
#> 189         -0.01071538         0.96958098           0.3573194
#> 190         -0.01071538         0.96958098           0.3573194
#> 191         -0.01071538         0.96958098           0.3573194
#> 192         -0.01071538         0.96958098           0.3573194
#> 193         -0.01071538         0.96958098           0.3573194
#> 194         -0.01071538         0.96958098           0.3573194
#> 195         -0.01071538         0.96958098           0.3573194
#> 196         -0.01071538         0.96958098           0.3573194
#> 197         -0.01071538         0.96958098           0.3573194
#> 198         -0.01071538         0.96958098           0.3573194
#> 199         -0.01071538         0.96958098           0.3573194
#> 200         -0.01071538         0.96958098           0.3573194
#> 201         -0.01071538         0.96958098           0.3573194
#> 202         -0.01071538         0.96958098           0.3573194
#> 203         -0.01071538         0.96958098           0.3573194
#> 204         -0.01071538         0.96958098           0.3573194
#> 205         -0.01071538         0.96958098           0.3573194
#> 206         -0.01071538         0.96958098           0.3573194
#> 207         -0.01071538         0.96958098           0.3573194
#> 208         -0.01071538         0.96958098           0.3573194
#> 209         -0.01071538         0.96958098           0.3573194
#> 210         -0.01071538         0.96958098           0.3573194
#> 211         -0.01071538         0.96958098           0.3573194
#> 212         -0.01071538         0.96958098           0.3573194
#> 213         -0.01071538         0.96958098           0.3573194
#> 214         -0.01071538         0.96958098           0.3573194
#> 215         -0.01071538         0.96958098           0.3573194
#> 216         -0.01071538         0.96958098           0.3573194
#> 217         -0.01071538         0.96958098           0.3573194
#> 218         -0.01071538         0.96958098           0.3573194
#> 219         -0.01071538         0.96958098           0.3573194
#> 220         -0.01071538         0.96958098           0.3573194
#> 221         -0.01071538         0.96958098           0.3573194
#> 222         -0.01071538         0.96958098           0.3573194
#> 223         -0.01071538         0.96958098           0.3573194
#> 224         -0.01071538         0.96958098           0.3573194
#> 225         -0.01071538         0.96958098           0.3573194
#> 226         -0.01071538         0.96958098           0.3573194
#> 227         -0.01071538         0.96958098           0.3573194
#> 228         -0.01071538         0.96958098           0.3573194
#> 229         -0.01071538         0.96958098           0.3573194
#> 230         -0.01071538         0.96958098           0.3573194
#> 231         -0.01071538         0.96958098           0.3573194
#> 232         -0.01071538         0.96958098           0.3573194
#> 233         -0.01071538         0.96958098           0.3573194
#> 234         -0.01071538         0.96958098           0.3573194
#> 235         -0.01071538         0.96958098           0.3573194
#> 236         -0.01071538         0.96958098           0.3573194
#> 237         -0.01071538         0.96958098           0.3573194
#> 238         -0.01071538         0.96958098           0.3573194
#> 239         -0.01071538         0.96958098           0.3573194
#> 240          0.45256210         0.06346625           0.4623181
#> 241          0.45256210         0.06346625           0.4623181
#> 242          0.45256210         0.06346625           0.4623181
#> 243          0.45256210         0.06346625           0.4623181
#> 244          0.45256210         0.06346625           0.4623181
#> 245          0.45256210         0.06346625           0.4623181
#> 246          0.45256210         0.06346625           0.4623181
#> 247          0.45256210         0.06346625           0.4623181
#> 248          0.45256210         0.06346625           0.4623181
#> 249          0.45256210         0.06346625           0.4623181
#> 250          0.45256210         0.06346625           0.4623181
#> 251          0.45256210         0.06346625           0.4623181
#> 252          0.45256210         0.06346625           0.4623181
#> 253          0.45256210         0.06346625           0.4623181
#> 254          0.45256210         0.06346625           0.4623181
#> 255          0.45256210         0.06346625           0.4623181
#> 256          0.45256210         0.06346625           0.4623181
#> 257          0.45256210         0.06346625           0.4623181
#> 258          0.45256210         0.06346625           0.4623181
#> 259          0.45256210         0.06346625           0.4623181
#> 260          0.45256210         0.06346625           0.4623181
#> 261          0.45256210         0.06346625           0.4623181
#> 262          0.45256210         0.06346625           0.4623181
#> 263          0.45256210         0.06346625           0.4623181
#> 264          0.45256210         0.06346625           0.4623181
#> 265          0.45256210         0.06346625           0.4623181
#> 266          0.45256210         0.06346625           0.4623181
#> 267          0.45256210         0.06346625           0.4623181
#> 268          0.45256210         0.06346625           0.4623181
#> 269          0.45256210         0.06346625           0.4623181
#> 270          0.45256210         0.06346625           0.4623181
#> 271          0.45256210         0.06346625           0.4623181
#> 272          0.45256210         0.06346625           0.4623181
#> 273          0.45256210         0.06346625           0.4623181
#> 274          0.45256210         0.06346625           0.4623181
#> 275          0.45256210         0.06346625           0.4623181
#> 276          0.45256210         0.06346625           0.4623181
#> 277          0.45256210         0.06346625           0.4623181
#> 278          0.45256210         0.06346625           0.4623181
#> 279          0.45256210         0.06346625           0.4623181
#> 280          0.45256210         0.06346625           0.4623181
#> 281          0.45256210         0.06346625           0.4623181
#> 282          0.45256210         0.06346625           0.4623181
#> 283          0.45256210         0.06346625           0.4623181
#> 284          0.45256210         0.06346625           0.4623181
#> 285          0.45256210         0.06346625           0.4623181
#> 286          0.45256210         0.06346625           0.4623181
#> 287          0.45256210         0.06346625           0.4623181
#> 288          0.45256210         0.06346625           0.4623181
#> 289          0.45256210         0.06346625           0.4623181
#> 290          0.45256210         0.06346625           0.4623181
#> 291          0.45256210         0.06346625           0.4623181
#> 292          0.45256210         0.06346625           0.4623181
#> 293          0.45256210         0.06346625           0.4623181
#> 294          0.45256210         0.06346625           0.4623181
#> 295          0.45256210         0.06346625           0.4623181
#> 296          0.45256210         0.06346625           0.4623181
#> 297          0.45256210         0.06346625           0.4623181
#> 298          0.45256210         0.06346625           0.4623181
#> 299          0.45256210         0.06346625           0.4623181
#> 300          0.45256210         0.06346625           0.4623181
#> 301          0.45256210         0.06346625           0.4623181
#> 302          0.45256210         0.06346625           0.4623181
#> 303          0.45256210         0.06346625           0.4623181
#> 304          0.45256210         0.06346625           0.4623181
#> 305          0.45256210         0.06346625           0.4623181
#> 306          0.45256210         0.06346625           0.4623181
#> 307          0.45256210         0.06346625           0.4623181
#> 308          0.45256210         0.06346625           0.4623181
#> 309          0.45256210         0.06346625           0.4623181
#> 310          0.45256210         0.06346625           0.4623181
#> 311          0.45256210         0.06346625           0.4623181
#> 312          0.45256210         0.06346625           0.4623181
#> 313          0.45256210         0.06346625           0.4623181
#> 314          0.45256210         0.06346625           0.4623181
#> 315          0.45256210         0.06346625           0.4623181
#> 316          0.45256210         0.06346625           0.4623181
#> 317          0.45256210         0.06346625           0.4623181
#> 318          0.45256210         0.06346625           0.4623181
#> 319          0.45256210         0.06346625           0.4623181
#> 320          0.45256210         0.06346625           0.4623181
#> 321          0.45256210         0.06346625           0.4623181
#> 322          0.45256210         0.06346625           0.4623181
#> 323          0.45256210         0.06346625           0.4623181
#> 324          0.45256210         0.06346625           0.4623181
#> 325          0.45256210         0.06346625           0.4623181
#> 326          0.45256210         0.06346625           0.4623181
#> 327          0.45256210         0.06346625           0.4623181
#> 328          0.45256210         0.06346625           0.4623181
#> 329          0.45256210         0.06346625           0.4623181
#> 330          0.45256210         0.06346625           0.4623181
#> 331          0.45256210         0.06346625           0.4623181
#> 332          0.45256210         0.06346625           0.4623181
#> 333          0.45256210         0.06346625           0.4623181
#> 334          0.45256210         0.06346625           0.4623181
#> 335          0.45256210         0.06346625           0.4623181
#> 336          0.45256210         0.06346625           0.4623181
#> 337          0.45256210         0.06346625           0.4623181
#> 338          0.45256210         0.06346625           0.4623181
#> 339          0.45256210         0.06346625           0.4623181
#> 340          0.45256210         0.06346625           0.4623181
#> 341          0.45256210         0.06346625           0.4623181
#> 342          0.45256210         0.06346625           0.4623181
#> 343          0.45256210         0.06346625           0.4623181
#> 344          0.45256210         0.06346625           0.4623181
#> 345          0.45256210         0.06346625           0.4623181
#> 346          0.45256210         0.06346625           0.4623181
#> 347          0.45256210         0.06346625           0.4623181
#> 348          0.45256210         0.06346625           0.4623181
#> 349          0.45256210         0.06346625           0.4623181
#> 350          0.45256210         0.06346625           0.4623181
#> 351          0.45256210         0.06346625           0.4623181
#> 352          0.45256210         0.06346625           0.4623181
#> 353          0.45256210         0.06346625           0.4623181
#> 354          0.45256210         0.06346625           0.4623181
#> 355          0.45256210         0.06346625           0.4623181
#> 356          0.45256210         0.06346625           0.4623181
#> 357          0.45256210         0.06346625           0.4623181
#> 358          0.45256210         0.06346625           0.4623181
#> 359          0.45256210         0.06346625           0.4623181
#> 360          0.45256210         0.06346625           0.4623181
```

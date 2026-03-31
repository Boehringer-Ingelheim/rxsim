# Example 5: Multi-arm \| Fixed design \| Single continuous endpoint \| Dunnett test + MCP-Mod

``` r
# Core simulation framework (Timer, Population, Trial, gen_timepoints, add_timepoints, ...)
library(rxsim)

# Analyses
library(multcomp)     # Dunnett
#> Loading required package: mvtnorm
#> Loading required package: survival
#> Loading required package: TH.data
#> Loading required package: MASS
#> 
#> Attaching package: 'TH.data'
#> The following object is masked from 'package:MASS':
#> 
#>     geyser
library(DoseFinding)  # MCP-Mod

# Helpers
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following object is masked from 'package:MASS':
#> 
#>     select
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
set.seed(4566)
```

Dose-finding trials aim to identify the dose-response relationship and
establish an effective dose level. This example simulates a multi-arm
fixed design across a placebo and four active doses, with the endpoint
generated from an Emax dose-response model. At the final analysis, two
methods are applied: a Dunnett test (comparing each active dose against
placebo while controlling family-wise error rate) and MCP-Mod (a
model-based contrast approach that selects the best-fitting
dose-response shape from a candidate set).

## Scenario

A **multi-arm, fixed design** (placebo + 4 active doses) with a **single
continuous endpoint**. At the **final analysis**, we perform:

- **Dunnett test**: all active doses vs placebo using a normal-theory
  linear model.
- **MCP-Mod**: model-based multiple contrast test + model fitting on a
  candidate set of dose–response shapes.

The Emax model generates mean responses via
`E(d) = e0 + emax × d / (ed50 + d)`, where `e0 = 0` (placebo baseline),
`emax = 1` (maximum achievable effect), and `ed50 = 20` (dose at
half-maximum effect). `arm_names = paste0("d", doses)` creates
human-readable labels (`d0`, `d5`, `d10`, `d20`, `d50`) that propagate
through trial data and analysis results. `delta = 0.1` is the minimum
effect size of clinical relevance used by MCP-Mod, and `alpha = 0.05`
controls the family-wise Type I error rate for both testing procedures.

``` r
# Dose levels (placebo + 4 actives)
doses <- c(0, 5, 10, 20, 50)

# Total N and allocation (balanced across arms)
sample_size <- 150
allocation  <- rep(1, length(doses))
arm_names   <- paste0("d", doses)  # e.g., d0, d5, ... used as arm labels

# Enrollment & dropout profile (piece-wise constant)
enrollment <- list(
  end_time = c(4, 8, 12, 16),
  rate     = c(6, 12, 18, 24)
)
dropout <- list(
  end_time = c(4, 8, 12, 16),
  rate     = c(0,  1,  2,  3)
)

# Data-generating model: Emax with homoscedastic noise
e0    <- 0.0
emax  <- 1.0
ed50  <- 20.0
sigma <- 1.0

mean_fun <- function(d) e0 + emax * d / (ed50 + d)

# operating characteristics
alpha <- 0.05
delta <- 0.1

scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation = list(allocation),
  n_arms = length(doses),
  alpha = alpha,
  delta = delta
)
```

## Time points

[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)
generates a piecewise-constant enrollment schedule with rates that
increase in steps (6, 12, 18, 24 subjects per unit time), reflecting a
realistic ramp-up in site activation. `tr_timer$get_end_timepoint()`
returns the earliest calendar time at which all `sample_size` subjects
are expected to be enrolled — this value is used as the trigger time for
the final analysis.

``` r
timepoints <- gen_timepoints(
  sample_size = sample_size,
  arms        = arm_names,
  allocation  = allocation,
  enrollment  = enrollment,
  dropout     = dropout
)

# Timer
tr_timer <- Timer$new(name = "timer_multiarm")
add_timepoints(tr_timer, timepoints)
final_time <- tr_timer$get_end_timepoint()
final_time
#> [1] 14
```

## Populations

`mk_pop_gen` is a closure factory: it captures the dose value `d` and
returns a generator function that, when called with `n`, draws `n`
responses from N(E(d), σ²). Each arm’s generator stores `dose` as a
numeric column because MCP-Mod requires actual dose values — not just
arm labels — to fit dose-response models and evaluate contrasts at
analysis time.

``` r
mk_pop_gen <- function(d) {
  function(n) {
    mu <- mean_fun(d)
    y  <- rnorm(n, mean = mu, sd = sigma)
    data.frame(
      id = seq_len(n),
      dose = d,
      y = y,
      readout_time = 1
    )
  }
}

population_generators <- lapply(seq_along(doses), function(i) {
  mk_pop_gen(doses[i])
})
names(population_generators) <- arm_names
```

## Trial Parameters

`enrollment_fn` and `dropout_fn` supply the inter-arrival and
time-to-dropout distributions for individual subjects. These are passed
to
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md),
which uses them internally when building each replicate’s enrollment
schedule.

``` r
enrollment_fn <- function(n) rexp(n, rate = 1)
dropout_fn <- function(n) rexp(n, rate = 0.01)
```

## Triggers & Analyses

At the **final time**, run:

1.  **Dunnett test** (active vs placebo) using
    [`multcomp::glht`](https://rdrr.io/pkg/multcomp/man/glht.html) on
    `lm(y ~ arm)`.

2.  **MCP-Mod** using
    [`DoseFinding::MCPMod`](https://openpharma.github.io/DoseFinding/reference/MCPMod.html)
    with a candidate model set (`Mods`).

The Dunnett test
([`multcomp::glht`](https://rdrr.io/pkg/multcomp/man/glht.html) with
`mcp(arm = "Dunnett")`) simultaneously compares each active dose against
placebo while controlling the family-wise error rate at `alpha`. MCP-Mod
is run with five candidate dose-response shapes; `selModel = "aveAIC"`
selects the best model by AIC-weighted averaging, and `Delta = delta`
sets the minimum effect size of clinical relevance for the contrast
step. `mct_min_p` extracts the minimum p-value across all candidate
model contrast tests — a small value indicates that at least one
candidate model detects a dose-response signal.

``` r
# Candidate model set for MCP-Mod
mods <- Mods(
  linear     = NULL,
  emax       = 20,
  exponential= 50,
  sigEmax    = c(20, 3),
  quadratic  = -0.2,
  doses      = doses
)

analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer){
      df_e <- df |>
        dplyr::filter(!is.na(enroll_time)) |>
        dplyr::mutate(
          arm  = factor(arm, levels = paste0("d", doses)),
          dose = as.numeric(dose)
        )

      # 1) Dunnett (active vs placebo)
      fit <- lm(y ~ arm, data = df_e)
      dun  <- multcomp::glht(fit, linfct = multcomp::mcp(arm = "Dunnett"))
      summ <- summary(dun)

      # 2) MCP-Mod (one-step)
      mm <- DoseFinding::MCPMod(
          dose   = df_e$dose,
          resp   = df_e$y,
          models = mods,
          type   = "normal",
          Delta  = delta,
          alpha  = alpha,
          selModel = "aveAIC"
      )

      mct_min_p <- min(attr(mm$MCTtest$tStat, "pVal"), na.rm = TRUE)

      data.frame(
        scenario,
        n_total      = nrow(df_e),
        dunn_min_p   = summ$test$pvalues |> as.numeric() |> min(),
        mcpmod_min_p = mct_min_p,
        stringsAsFactors = FALSE
      )
    }
  )
)
```

## Trial

``` r
trials <- replicate_trial(
  trial_name = "multiarm_dunnett_mcpmod",
  sample_size = sample_size,
  arms = arm_names,
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
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> Warning in MCTpval(contMat, corMat, df, tStat, alternative, mvtcontrol):
#> Warning from mvtnorm::pmvt: Completion with error > abseps.
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: multiarm_dunnett_mcpmod_1
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
#>     name: multiarm_dunnett_mcpmod_2
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
#>     name: multiarm_dunnett_mcpmod_3
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
```

## Results

[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
row-binds analysis outputs across all replicates and prepends three
bookkeeping columns: `replicate` (integer index), `timepoint` (calendar
time at which the analysis fired), and `analysis` (the analysis name,
here `"final"`). `dunn_min_p` is the smallest Dunnett-adjusted p-value
across the four active-dose comparisons — a value below `alpha`
indicates that at least one dose is significantly different from placebo
after multiplicity correction. `mcpmod_min_p` is the minimum MCP-Mod
contrast test p-value across candidate models; both metrics test for the
presence of a dose-response signal, but MCP-Mod is generally more
powerful when the true dose-response shape is well-captured by one of
the candidate models. See [Core
Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.md)
for background on trial parameters.

One row per replicate, row-bound across all replicates:

``` r
replicate_results <- collect_results(trials)
replicate_results
#>     replicate  timepoint analysis sample_size    allocation n_arms alpha delta
#> 1           1   143.0045    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 2           1   339.9395    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 3           1   360.7831    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 4           1   588.7844    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 5           1   670.0388    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 6           1   772.1356    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 7           1   775.1476    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 8           1   787.1425    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 9           1   798.6437    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 10          1  1028.4725    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 11          1  1080.1645    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 12          1  1253.6509    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 13          1  1317.8008    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 14          1  1397.0989    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 15          1  1492.6022    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 16          1  1775.4050    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 17          1  1887.0859    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 18          1  1915.6536    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 19          1  1930.4359    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 20          1  2106.0079    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 21          1  2153.1301    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 22          1  2154.5230    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 23          1  2159.6204    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 24          1  2246.0420    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 25          1  2284.4002    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 26          1  2394.9954    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 27          1  2437.6287    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 28          1  2988.2379    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 29          1  3058.7220    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 30          1  3073.4531    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 31          1  3108.8572    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 32          1  3318.4561    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 33          1  3331.7474    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 34          1  3449.3426    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 35          1  3541.0817    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 36          1  3710.8414    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 37          1  3892.1890    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 38          1  3965.8668    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 39          1  4212.2229    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 40          1  4374.7451    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 41          1  4443.1475    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 42          1  4654.6850    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 43          1  4750.0769    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 44          1  4778.3686    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 45          1  4786.5390    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 46          1  4786.6340    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 47          1  5027.6661    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 48          1  5052.1873    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 49          1  5371.1915    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 50          1  5435.0975    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 51          1  5540.2019    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 52          1  5567.7204    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 53          1  5663.2639    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 54          1  5722.4476    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 55          1  5767.2693    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 56          1  5900.9525    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 57          1  6228.5340    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 58          1  6328.1474    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 59          1  6361.2112    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 60          1  6450.6320    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 61          1  6499.9763    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 62          1  6519.2181    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 63          1  6572.0482    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 64          1  6847.9736    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 65          1  6854.6817    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 66          1  6999.4844    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 67          1  7048.1833    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 68          1  7143.2179    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 69          1  7156.0925    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 70          1  7323.9579    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 71          1  7542.9430    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 72          1  7647.3063    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 73          1  7710.8601    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 74          1  7857.1631    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 75          1  7892.4612    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 76          1  7989.6922    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 77          1  8021.6614    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 78          1  8066.3794    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 79          1  8093.6099    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 80          1  8315.4158    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 81          1  8359.9900    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 82          1  8456.6798    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 83          1  8691.9279    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 84          1  8700.0742    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 85          1  8715.5817    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 86          1  8739.9677    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 87          1  9030.9311    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 88          1  9070.9733    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 89          1  9229.1156    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 90          1  9236.6338    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 91          1  9287.9738    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 92          1  9296.3791    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 93          1  9319.3655    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 94          1  9327.6958    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 95          1  9404.7676    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 96          1  9458.4023    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 97          1  9498.5916    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 98          1  9523.2242    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 99          1  9543.2741    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 100         1  9627.7292    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 101         1  9777.8174    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 102         1  9798.5206    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 103         1 10006.4467    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 104         1 10007.1916    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 105         1 10020.3462    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 106         1 10033.1938    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 107         1 10129.7089    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 108         1 10132.8448    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 109         1 10183.4848    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 110         1 10262.3121    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 111         1 10298.2363    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 112         1 10381.1510    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 113         1 10449.3668    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 114         1 10595.6510    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 115         1 10599.5382    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 116         1 10994.7581    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 117         1 11156.1662    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 118         1 11245.4584    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 119         1 11370.0855    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 120         1 11464.1854    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 121         1 11489.9480    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 122         1 11615.6723    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 123         1 11818.7085    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 124         1 11896.2410    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 125         1 12035.7341    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 126         1 12155.2739    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 127         1 12171.0963    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 128         1 12435.5286    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 129         1 12721.4746    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 130         1 12785.5565    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 131         1 12991.3102    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 132         1 13073.8836    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 133         1 13238.9815    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 134         1 13382.9275    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 135         1 13433.1607    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 136         1 13491.5979    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 137         1 13613.6308    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 138         1 13620.4157    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 139         1 13738.3349    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 140         1 13763.0853    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 141         1 13779.3080    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 142         1 13990.4016    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 143         1 14018.8871    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 144         1 14037.5864    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 145         1 14040.5193    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 146         1 14113.0028    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 147         1 14211.0777    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 148         1 14232.8339    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 149         1 14242.4501    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 150         1 14256.6441    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 151         1 14260.4322    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 152         2   148.2639    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 153         2   183.0158    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 154         2   228.1675    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 155         2   272.9651    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 156         2   540.9145    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 157         2   591.6585    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 158         2   653.8085    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 159         2   693.7968    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 160         2   696.8845    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 161         2   796.6384    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 162         2   815.3134    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 163         2   943.4248    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 164         2   988.3653    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 165         2  1041.9727    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 166         2  1110.7420    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 167         2  1596.0697    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 168         2  1680.2417    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 169         2  1917.6488    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 170         2  1923.6892    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 171         2  2185.3293    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 172         2  2309.3059    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 173         2  2516.3584    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 174         2  2579.2572    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 175         2  2658.6346    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 176         2  2672.1026    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 177         2  2852.5090    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 178         2  3227.3399    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 179         2  3267.4390    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 180         2  3346.3101    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 181         2  3354.5941    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 182         2  3377.7399    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 183         2  3411.9433    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 184         2  3467.1993    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 185         2  3471.0882    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 186         2  3501.8996    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 187         2  3658.3615    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 188         2  3767.1029    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 189         2  4328.9796    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 190         2  4427.8640    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 191         2  4443.4503    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 192         2  4592.3509    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 193         2  4631.3513    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 194         2  4703.2820    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 195         2  4778.0570    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 196         2  4807.1763    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 197         2  4940.0382    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 198         2  5096.0842    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 199         2  5365.4437    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 200         2  5372.2570    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 201         2  5432.5425    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 202         2  5577.3849    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 203         2  5682.7839    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 204         2  5688.0523    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 205         2  5718.7936    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 206         2  5910.7442    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 207         2  6048.7189    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 208         2  6130.4166    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 209         2  6349.9997    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 210         2  6425.3210    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 211         2  6506.7830    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 212         2  6572.6640    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 213         2  6586.1021    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 214         2  6784.9064    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 215         2  6800.7430    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 216         2  6872.1424    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 217         2  6930.0081    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 218         2  7141.8191    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 219         2  7243.3190    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 220         2  7352.9040    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 221         2  7354.0977    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 222         2  7398.4770    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 223         2  7496.9767    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 224         2  7572.9416    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 225         2  7635.8752    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 226         2  7729.2110    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 227         2  7748.9087    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 228         2  7807.0321    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 229         2  7821.2778    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 230         2  7855.2292    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 231         2  8049.9707    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 232         2  8111.1810    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 233         2  8165.1570    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 234         2  8289.7631    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 235         2  8320.9267    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 236         2  8432.9956    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 237         2  8483.3301    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 238         2  8680.5709    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 239         2  8882.3747    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 240         2  9181.1285    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 241         2  9186.0610    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 242         2  9263.9555    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 243         2  9269.8875    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 244         2  9281.1404    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 245         2  9310.2995    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 246         2  9321.7835    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 247         2  9351.1467    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 248         2  9364.9909    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 249         2  9415.6142    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 250         2  9530.2165    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 251         2  9556.7782    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 252         2  9666.0535    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 253         2  9704.6674    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 254         2  9773.5464    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 255         2  9848.4262    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 256         2 10056.3384    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 257         2 10285.2922    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 258         2 10294.8585    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 259         2 10330.9952    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 260         2 10398.9500    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 261         2 10402.5365    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 262         2 10426.6804    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 263         2 10461.3726    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 264         2 10515.8574    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 265         2 10538.3812    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 266         2 10578.5037    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 267         2 10619.7660    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 268         2 10835.3669    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 269         2 10837.6253    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 270         2 10893.6843    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 271         2 10959.6185    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 272         2 11006.8909    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 273         2 11029.1307    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 274         2 11192.1073    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 275         2 11255.5176    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 276         2 11287.0616    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 277         2 11331.7036    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 278         2 11486.5331    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 279         2 11503.8151    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 280         2 11675.1681    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 281         2 11691.3073    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 282         2 11780.5489    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 283         2 11792.4402    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 284         2 11793.7679    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 285         2 11855.7885    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 286         2 11939.9152    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 287         2 12009.2665    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 288         2 12015.3754    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 289         2 12022.3184    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 290         2 12150.9301    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 291         2 12181.1741    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 292         2 12284.9338    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 293         2 12313.8349    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 294         2 12358.7083    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 295         2 12388.1120    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 296         2 12397.8754    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 297         2 12485.5591    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 298         2 12528.1637    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 299         2 12649.2750    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 300         2 12661.4331    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 301         2 12665.4148    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 302         3   155.7178    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 303         3   193.8827    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 304         3   198.8823    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 305         3   334.7364    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 306         3   339.4752    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 307         3   383.9209    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 308         3   537.5278    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 309         3   562.2469    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 310         3   921.6241    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 311         3   960.1357    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 312         3  1346.5930    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 313         3  1523.0592    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 314         3  1762.7628    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 315         3  1826.2814    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 316         3  1867.6049    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 317         3  2031.3598    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 318         3  2048.9229    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 319         3  2089.7161    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 320         3  2121.0194    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 321         3  2125.9011    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 322         3  2231.3755    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 323         3  2425.7877    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 324         3  2435.1934    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 325         3  2460.0795    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 326         3  2484.3196    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 327         3  2662.3748    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 328         3  2765.7353    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 329         3  2834.8390    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 330         3  2845.9274    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 331         3  2863.8990    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 332         3  2942.2235    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 333         3  3105.6478    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 334         3  3160.3356    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 335         3  3163.2174    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 336         3  3227.5468    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 337         3  3367.9410    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 338         3  3368.1766    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 339         3  3377.5907    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 340         3  3398.9228    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 341         3  3493.9042    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 342         3  3591.6655    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 343         3  3614.7473    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 344         3  3944.1708    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 345         3  3960.5114    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 346         3  3965.5005    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 347         3  4038.3088    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 348         3  4060.5184    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 349         3  4106.0463    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 350         3  4208.9015    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 351         3  4368.2082    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 352         3  4418.3022    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 353         3  4780.7164    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 354         3  4810.8559    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 355         3  4849.7158    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 356         3  4988.2681    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 357         3  5058.4802    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 358         3  5193.0267    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 359         3  5206.2424    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 360         3  5381.1219    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 361         3  5420.4016    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 362         3  5631.9827    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 363         3  5797.2623    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 364         3  5841.9415    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 365         3  5933.9986    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 366         3  5974.8873    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 367         3  6000.2834    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 368         3  6107.3078    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 369         3  6332.0390    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 370         3  6487.8059    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 371         3  6794.3180    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 372         3  6798.8060    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 373         3  6846.9087    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 374         3  6893.8417    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 375         3  7045.7982    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 376         3  7113.5309    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 377         3  7302.2290    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 378         3  7481.4339    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 379         3  7740.4727    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 380         3  7880.2549    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 381         3  7978.3170    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 382         3  8013.5910    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 383         3  8067.0609    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 384         3  8116.7020    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 385         3  8324.5649    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 386         3  8398.6764    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 387         3  8581.4503    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 388         3  8623.9674    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 389         3  8652.9426    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 390         3  8680.9608    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 391         3  8714.9963    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 392         3  8734.4351    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 393         3  8756.4415    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 394         3  8779.9915    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 395         3  8841.0847    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 396         3  8907.0258    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 397         3  8960.7670    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 398         3  9100.2278    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 399         3  9148.8716    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 400         3  9260.3827    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 401         3  9445.5819    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 402         3  9461.3438    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 403         3  9499.8874    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 404         3  9525.3776    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 405         3  9571.3079    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 406         3  9651.3803    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 407         3  9685.9553    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 408         3  9756.1710    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 409         3  9847.4343    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 410         3 10050.0060    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 411         3 10051.7798    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 412         3 10334.4962    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 413         3 10390.2329    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 414         3 11174.2808    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 415         3 11177.7385    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 416         3 11254.0607    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 417         3 11425.8466    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 418         3 11573.3197    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 419         3 11597.4206    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 420         3 11781.5131    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 421         3 11804.2088    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 422         3 11875.9285    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 423         3 11951.4247    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 424         3 11999.7391    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 425         3 12029.2570    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 426         3 12037.5722    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 427         3 12418.5675    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 428         3 12444.1163    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 429         3 12464.3175    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 430         3 12482.7694    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 431         3 12597.8750    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 432         3 12718.4910    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 433         3 13169.1349    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 434         3 13330.4794    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 435         3 13339.2446    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 436         3 13532.5844    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 437         3 13591.1442    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 438         3 13601.6034    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 439         3 13659.2478    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 440         3 13891.0138    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 441         3 13935.3991    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 442         3 14349.8417    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 443         3 14545.6029    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 444         3 14627.2554    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 445         3 14636.8145    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 446         3 14639.2416    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 447         3 14665.9031    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 448         3 14718.3569    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 449         3 14874.7654    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 450         3 15081.0408    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 451         3 15189.6279    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#>     n_total   dunn_min_p mcpmod_min_p
#> 1       150 2.614327e-01 9.360871e-02
#> 2       150 2.614396e-01 9.275725e-02
#> 3       150 2.613955e-01 9.319514e-02
#> 4       150 2.614196e-01 9.270148e-02
#> 5       150 2.614895e-01 9.276225e-02
#> 6       150 2.614892e-01 9.345740e-02
#> 7       150 2.614015e-01 9.306070e-02
#> 8       150 2.615511e-01 9.272428e-02
#> 9       150 2.614361e-01 9.244525e-02
#> 10      150 2.614332e-01 9.290695e-02
#> 11      150 2.614584e-01 9.308355e-02
#> 12      150 2.614969e-01 9.253591e-02
#> 13      150 2.614508e-01 9.292989e-02
#> 14      150 2.614040e-01 9.363291e-02
#> 15      150 2.614044e-01 9.262052e-02
#> 16      150 2.615193e-01 9.286933e-02
#> 17      150 2.614002e-01 9.284280e-02
#> 18      150 2.614238e-01 9.340432e-02
#> 19      150 2.614416e-01 9.290941e-02
#> 20      150 2.616167e-01 9.323894e-02
#> 21      150 2.614347e-01 9.300686e-02
#> 22      150 2.613603e-01 9.293349e-02
#> 23      150 2.614027e-01 9.298311e-02
#> 24      150 2.614211e-01 9.312740e-02
#> 25      150 2.613984e-01 9.349511e-02
#> 26      150 2.614601e-01 9.277710e-02
#> 27      150 2.614322e-01 9.265344e-02
#> 28      150 2.615036e-01 9.326216e-02
#> 29      150 2.613462e-01 9.300571e-02
#> 30      150 2.614575e-01 9.273741e-02
#> 31      150 2.614319e-01 9.316567e-02
#> 32      150 2.613480e-01 9.312085e-02
#> 33      150 2.614617e-01 9.343095e-02
#> 34      150 2.614056e-01 9.312477e-02
#> 35      150 2.615702e-01 9.300466e-02
#> 36      150 2.614963e-01 9.327955e-02
#> 37      150 2.614518e-01 9.329042e-02
#> 38      150 2.614567e-01 9.273575e-02
#> 39      150 2.615624e-01 9.326930e-02
#> 40      150 2.613972e-01 9.303780e-02
#> 41      150 2.614363e-01 9.299402e-02
#> 42      150 2.614751e-01 9.284153e-02
#> 43      150 2.614512e-01 9.295286e-02
#> 44      150 2.614149e-01 9.237586e-02
#> 45      150 2.614151e-01 9.278536e-02
#> 46      150 2.614348e-01 9.250590e-02
#> 47      150 2.614754e-01 9.298851e-02
#> 48      150 2.613954e-01 9.314493e-02
#> 49      150 2.615101e-01 9.309354e-02
#> 50      150 2.615213e-01 9.347243e-02
#> 51      150 2.614228e-01 9.281253e-02
#> 52      150 2.614546e-01 9.329139e-02
#> 53      150 2.615069e-01 9.290822e-02
#> 54      150 2.614887e-01 9.343875e-02
#> 55      150 2.614446e-01 9.310967e-02
#> 56      150 2.614339e-01 9.302447e-02
#> 57      150 2.615390e-01 9.288955e-02
#> 58      150 2.614182e-01 9.292151e-02
#> 59      150 2.614378e-01 9.322694e-02
#> 60      150 2.614630e-01 9.310439e-02
#> 61      150 2.614645e-01 9.317999e-02
#> 62      150 2.614052e-01 9.274885e-02
#> 63      150 2.614066e-01 9.316136e-02
#> 64      150 2.613996e-01 9.305733e-02
#> 65      150 2.615829e-01 9.311495e-02
#> 66      150 2.613643e-01 9.326818e-02
#> 67      150 2.613852e-01 9.333006e-02
#> 68      150 2.615034e-01 9.350892e-02
#> 69      150 2.614114e-01 9.306649e-02
#> 70      150 2.613583e-01 9.301616e-02
#> 71      150 2.613885e-01 9.298870e-02
#> 72      150 2.614391e-01 9.286583e-02
#> 73      150 2.614707e-01 9.323755e-02
#> 74      150 2.615086e-01 9.326162e-02
#> 75      150 2.614686e-01 9.281220e-02
#> 76      150 2.613521e-01 9.321166e-02
#> 77      150 2.614946e-01 9.328788e-02
#> 78      150 2.614557e-01 9.268851e-02
#> 79      150 2.614772e-01 9.309326e-02
#> 80      150 2.614987e-01 9.313837e-02
#> 81      150 2.614306e-01 9.313236e-02
#> 82      150 2.614124e-01 9.361133e-02
#> 83      150 2.615002e-01 9.320870e-02
#> 84      150 2.614384e-01 9.335957e-02
#> 85      150 2.613176e-01 9.335430e-02
#> 86      150 2.614249e-01 9.288760e-02
#> 87      150 2.615418e-01 9.320970e-02
#> 88      150 2.614237e-01 9.285773e-02
#> 89      150 2.614312e-01 9.320397e-02
#> 90      150 2.614070e-01 9.296978e-02
#> 91      150 2.614820e-01 9.337207e-02
#> 92      150 2.614042e-01 9.321349e-02
#> 93      150 2.614723e-01 9.309716e-02
#> 94      150 2.614290e-01 9.382504e-02
#> 95      150 2.613988e-01 9.268601e-02
#> 96      150 2.615152e-01 9.282868e-02
#> 97      150 2.614500e-01 9.241847e-02
#> 98      150 2.614159e-01 9.328178e-02
#> 99      150 2.614869e-01 9.309507e-02
#> 100     150 2.614638e-01 9.291479e-02
#> 101     150 2.615125e-01 9.244476e-02
#> 102     150 2.614196e-01 9.135931e-02
#> 103     150 2.613786e-01 9.295149e-02
#> 104     150 2.614485e-01 9.297123e-02
#> 105     150 2.614764e-01 9.345298e-02
#> 106     150 2.613265e-01 9.319463e-02
#> 107     150 2.615018e-01 9.310204e-02
#> 108     150 2.614128e-01 9.269482e-02
#> 109     150 2.614334e-01 9.318676e-02
#> 110     150 2.613965e-01 9.323312e-02
#> 111     150 2.614217e-01 9.297889e-02
#> 112     150 2.613930e-01 9.354575e-02
#> 113     150 2.614780e-01 9.297310e-02
#> 114     150 2.614366e-01 9.295070e-02
#> 115     150 2.613386e-01 9.309881e-02
#> 116     150 2.614538e-01 9.306822e-02
#> 117     150 2.614093e-01 9.320874e-02
#> 118     150 2.613804e-01 9.338952e-02
#> 119     150 2.614565e-01 9.331719e-02
#> 120     150 2.614653e-01 9.296052e-02
#> 121     150 2.613838e-01 9.416520e-02
#> 122     150 2.613849e-01 9.308690e-02
#> 123     150 2.615031e-01 9.303579e-02
#> 124     150 2.615844e-01 9.244391e-02
#> 125     150 2.614941e-01 9.304041e-02
#> 126     150 2.615478e-01 9.272218e-02
#> 127     150 2.615309e-01 9.294953e-02
#> 128     150 2.614171e-01 9.315036e-02
#> 129     150 2.614899e-01 9.329863e-02
#> 130     150 2.613846e-01 9.289970e-02
#> 131     150 2.614287e-01 9.300815e-02
#> 132     150 2.613020e-01 9.335157e-02
#> 133     150 2.614981e-01 9.296093e-02
#> 134     150 2.615437e-01 9.284297e-02
#> 135     150 2.615497e-01 9.323686e-02
#> 136     150 2.613633e-01 9.267196e-02
#> 137     150 2.614964e-01 9.336882e-02
#> 138     150 2.614169e-01 9.309549e-02
#> 139     150 2.613611e-01 9.278984e-02
#> 140     150 2.613671e-01 9.266914e-02
#> 141     150 2.613687e-01 9.337392e-02
#> 142     150 2.614785e-01 9.289227e-02
#> 143     150 2.614767e-01 9.311363e-02
#> 144     150 2.615267e-01 9.311493e-02
#> 145     150 2.613882e-01 9.299885e-02
#> 146     150 2.614772e-01 9.313427e-02
#> 147     150 2.615392e-01 9.319678e-02
#> 148     150 2.614104e-01 9.313141e-02
#> 149     150 2.615117e-01 9.273186e-02
#> 150     150 2.615516e-01 9.237583e-02
#> 151     150 2.613845e-01 9.217027e-02
#> 152     150 8.981219e-05 1.123651e-05
#> 153     150 7.961448e-05 3.660842e-06
#> 154     150 7.860573e-05 3.293457e-06
#> 155     150 9.792402e-05 4.038984e-06
#> 156     150 9.217714e-05 3.253847e-06
#> 157     150 9.473143e-05 3.105401e-06
#> 158     150 9.080252e-05 3.706229e-06
#> 159     150 1.052573e-04 4.266154e-06
#> 160     150 8.708270e-05 4.546306e-06
#> 161     150 8.849276e-05 5.623681e-06
#> 162     150 8.256923e-05 4.426698e-06
#> 163     150 8.535531e-05 4.387355e-06
#> 164     150 9.959180e-05 3.542223e-06
#> 165     150 9.606843e-05 3.174445e-06
#> 166     150 8.546491e-05 3.907943e-06
#> 167     150 7.335827e-05 4.771371e-06
#> 168     150 7.789603e-05 3.058025e-06
#> 169     150 8.127999e-05 6.561985e-06
#> 170     150 7.613346e-05 7.149066e-06
#> 171     150 8.592316e-05 5.534546e-06
#> 172     150 8.550120e-05 3.156882e-06
#> 173     150 9.074101e-05 1.225210e-05
#> 174     150 7.406792e-05 3.494359e-06
#> 175     150 8.354890e-05 4.180008e-06
#> 176     150 8.195689e-05 3.841123e-06
#> 177     150 9.993704e-05 5.589129e-06
#> 178     150 9.753898e-05 2.961110e-06
#> 179     150 7.338543e-05 1.169763e-05
#> 180     150 8.659199e-05 3.351433e-06
#> 181     150 8.271319e-05 3.048907e-06
#> 182     150 9.544456e-05 4.022911e-06
#> 183     150 9.664125e-05 2.924540e-06
#> 184     150 9.100926e-05 3.382574e-06
#> 185     150 9.213875e-05 3.625582e-06
#> 186     150 8.541382e-05 3.410981e-06
#> 187     150 8.638013e-05 3.197464e-06
#> 188     150 7.753415e-05 3.212647e-06
#> 189     150 8.848871e-05 1.015368e-05
#> 190     150 9.241812e-05 3.391715e-06
#> 191     150 8.070419e-05 4.568120e-06
#> 192     150 8.982265e-05 4.710213e-06
#> 193     150 7.849511e-05 4.366608e-06
#> 194     150 9.924447e-05 9.715528e-06
#> 195     150 8.329449e-05 3.880861e-06
#> 196     150 8.000101e-05 3.295047e-06
#> 197     150 8.994803e-05 1.298347e-05
#> 198     150 8.038019e-05 3.136372e-06
#> 199     150 8.575402e-05 2.842082e-06
#> 200     150 1.034467e-04 3.987165e-06
#> 201     150 7.850827e-05 6.031850e-06
#> 202     150 8.247827e-05 1.061563e-05
#> 203     150 7.637356e-05 2.952576e-06
#> 204     150 9.407388e-05 3.127547e-06
#> 205     150 9.673583e-05 3.317754e-06
#> 206     150 8.014738e-05 5.329364e-06
#> 207     150 8.566631e-05 3.649790e-06
#> 208     150 8.196794e-05 5.677936e-06
#> 209     150 8.112970e-05 5.267769e-06
#> 210     150 8.375053e-05 3.101002e-06
#> 211     150 9.351520e-05 3.065850e-06
#> 212     150 7.845088e-05 4.538759e-06
#> 213     150 8.673866e-05 3.077422e-06
#> 214     150 8.343704e-05 3.406318e-06
#> 215     150 1.085818e-04 4.300808e-06
#> 216     150 8.687553e-05 5.950416e-06
#> 217     150 9.085260e-05 5.127386e-06
#> 218     150 8.159712e-05 3.848200e-06
#> 219     150 8.710774e-05 4.249256e-06
#> 220     150 9.324710e-05 3.163779e-06
#> 221     150 8.445660e-05 3.742343e-06
#> 222     150 9.823264e-05 5.100253e-06
#> 223     150 7.781721e-05 3.298102e-06
#> 224     150 8.093244e-05 3.623732e-06
#> 225     150 8.155167e-05 3.016929e-06
#> 226     150 9.046981e-05 3.523167e-06
#> 227     150 9.854356e-05 3.057922e-06
#> 228     150 7.610552e-05 4.073044e-06
#> 229     150 9.194991e-05 3.324073e-06
#> 230     150 8.039539e-05 3.238537e-06
#> 231     150 7.908109e-05 3.412122e-06
#> 232     150 7.865262e-05 3.049822e-06
#> 233     150 8.774424e-05 3.380072e-06
#> 234     150 9.936927e-05 4.110546e-06
#> 235     150 8.256012e-05 3.157050e-06
#> 236     150 9.565554e-05 4.699858e-06
#> 237     150 8.394185e-05 2.986234e-06
#> 238     150 8.053769e-05 4.120730e-06
#> 239     150 9.178960e-05 4.732087e-06
#> 240     150 8.146128e-05 4.419217e-06
#> 241     150 8.896865e-05 3.323733e-06
#> 242     150 8.324750e-05 3.573718e-06
#> 243     150 7.727374e-05 3.156861e-06
#> 244     150 8.597741e-05 3.333968e-06
#> 245     150 9.901049e-05 3.434634e-06
#> 246     150 8.585519e-05 3.768781e-06
#> 247     150 9.592965e-05 3.144893e-06
#> 248     150 8.946895e-05 5.456388e-06
#> 249     150 1.036129e-04 4.441485e-06
#> 250     150 8.639761e-05 3.125239e-06
#> 251     150 8.862373e-05 3.395543e-06
#> 252     150 8.100858e-05 5.381122e-06
#> 253     150 8.249337e-05 4.294650e-06
#> 254     150 9.258027e-05 3.218535e-06
#> 255     150 7.995027e-05 5.870885e-06
#> 256     150 8.035010e-05 1.063567e-05
#> 257     150 7.883884e-05 4.726642e-06
#> 258     150 1.335054e-04 3.415576e-06
#> 259     150 9.059259e-05 3.031212e-06
#> 260     150 9.759776e-05 3.513580e-06
#> 261     150 8.427311e-05 9.716059e-06
#> 262     150 1.050883e-04 1.429249e-05
#> 263     150 8.671316e-05 4.649378e-06
#> 264     150 8.030417e-05 3.155409e-06
#> 265     150 1.092399e-04 3.029741e-06
#> 266     150 1.044078e-04 3.376731e-06
#> 267     150 1.049721e-04 5.558451e-06
#> 268     150 1.149923e-04 3.884535e-06
#> 269     150 8.230897e-05 2.907439e-06
#> 270     150 8.922891e-05 4.698365e-06
#> 271     150 8.542828e-05 3.246686e-06
#> 272     150 8.528368e-05 3.627737e-06
#> 273     150 8.781407e-05 3.899593e-06
#> 274     150 7.956818e-05 3.145309e-06
#> 275     150 1.063663e-04 4.189047e-06
#> 276     150 8.127860e-05 5.753535e-06
#> 277     150 8.841512e-05 2.889243e-06
#> 278     150 1.071952e-04 3.230267e-06
#> 279     150 7.531926e-05 4.160609e-06
#> 280     150 8.166653e-05 3.250126e-06
#> 281     150 8.731245e-05 4.065649e-06
#> 282     150 7.861762e-05 3.086610e-06
#> 283     150 7.869225e-05 3.788437e-06
#> 284     150 7.903125e-05 4.672738e-06
#> 285     150 1.021289e-04 3.142116e-06
#> 286     150 8.853428e-05 4.476558e-06
#> 287     150 7.863446e-05 4.278731e-06
#> 288     150 7.751441e-05 4.601287e-06
#> 289     150 9.068121e-05 3.025168e-06
#> 290     150 8.383178e-05 6.816018e-06
#> 291     150 8.212471e-05 6.802007e-06
#> 292     150 9.848770e-05 3.188885e-06
#> 293     150 9.288432e-05 5.531478e-06
#> 294     150 1.028726e-04 5.493731e-06
#> 295     150 9.532545e-05 3.416792e-06
#> 296     150 8.175583e-05 6.115948e-06
#> 297     150 9.853731e-05 4.144229e-06
#> 298     150 8.443076e-05 3.171534e-06
#> 299     150 8.534431e-05 3.931003e-06
#> 300     150 8.463621e-05 3.234322e-06
#> 301     150 8.525217e-05 5.311228e-06
#> 302     150 1.540937e-04 1.831518e-05
#> 303     150 1.486112e-04 2.105323e-05
#> 304     150 1.538273e-04 2.127895e-05
#> 305     150 1.422745e-04 1.702716e-05
#> 306     150 2.005675e-04 2.343646e-05
#> 307     150 1.485906e-04 3.806392e-05
#> 308     150 1.623462e-04 2.441022e-05
#> 309     150 1.503481e-04 3.733663e-05
#> 310     150 1.797116e-04 2.915798e-05
#> 311     150 1.614737e-04 2.068826e-05
#> 312     150 1.843080e-04 1.853905e-05
#> 313     150 1.907370e-04 2.714468e-05
#> 314     150 1.592947e-04 4.017872e-05
#> 315     150 1.531499e-04 2.694495e-05
#> 316     150 1.443367e-04 1.904073e-05
#> 317     150 1.509031e-04 1.738920e-05
#> 318     150 1.481064e-04 3.571296e-05
#> 319     150 1.653353e-04 2.341770e-05
#> 320     150 1.538056e-04 3.135866e-05
#> 321     150 1.553463e-04 2.817980e-05
#> 322     150 1.560652e-04 1.951457e-05
#> 323     150 1.438423e-04 2.858206e-05
#> 324     150 1.448520e-04 1.890106e-05
#> 325     150 1.483413e-04 2.692030e-05
#> 326     150 1.736704e-04 3.171798e-05
#> 327     150 1.797214e-04 1.754470e-05
#> 328     150 1.439721e-04 1.609126e-05
#> 329     150 1.620651e-04 1.818050e-05
#> 330     150 1.571943e-04 1.631641e-05
#> 331     150 1.542269e-04 1.734375e-05
#> 332     150 1.581470e-04 2.295099e-05
#> 333     150 1.666742e-04 2.268851e-05
#> 334     150 1.520787e-04 1.790566e-05
#> 335     150 1.704419e-04 3.194841e-05
#> 336     150 1.601266e-04 1.922853e-05
#> 337     150 1.539982e-04 1.182779e-04
#> 338     150 1.598785e-04 1.659750e-05
#> 339     150 1.424631e-04 1.545492e-05
#> 340     150 1.411589e-04 1.682579e-05
#> 341     150 1.677200e-04 1.570373e-05
#> 342     150 1.562820e-04 2.466652e-05
#> 343     150 1.720883e-04 2.126923e-05
#> 344     150 1.866546e-04 1.766676e-05
#> 345     150 1.551548e-04 1.730262e-05
#> 346     150 1.566000e-04 2.129428e-05
#> 347     150 1.640063e-04 2.122334e-05
#> 348     150 1.483623e-04 2.619601e-05
#> 349     150 1.458905e-04 1.893657e-05
#> 350     150 1.760133e-04 2.420377e-05
#> 351     150 1.739001e-04 1.656794e-05
#> 352     150 1.835615e-04 2.148656e-05
#> 353     150 1.453044e-04 3.058325e-05
#> 354     150 1.873121e-04 2.288007e-05
#> 355     150 1.620287e-04 2.044667e-05
#> 356     150 1.819241e-04 2.102506e-05
#> 357     150 1.591919e-04 1.656483e-05
#> 358     150 1.500343e-04 1.855311e-05
#> 359     150 1.701330e-04 2.657025e-05
#> 360     150 1.559658e-04 2.700287e-05
#> 361     150 1.470612e-04 2.069972e-05
#> 362     150 1.846304e-04 2.968998e-05
#> 363     150 2.044847e-04 1.750060e-05
#> 364     150 1.784955e-04 2.650080e-05
#> 365     150 1.927798e-04 2.904097e-05
#> 366     150 1.685548e-04 3.940434e-05
#> 367     150 1.466113e-04 2.428654e-05
#> 368     150 1.723505e-04 2.653453e-05
#> 369     150 1.463968e-04 1.725588e-05
#> 370     150 1.542035e-04 1.970338e-05
#> 371     150 1.540188e-04 2.063157e-05
#> 372     150 1.691332e-04 3.565796e-05
#> 373     150 1.865557e-04 1.797510e-05
#> 374     150 1.461914e-04 1.732894e-05
#> 375     150 1.867721e-04 1.958258e-05
#> 376     150 1.688553e-04 6.769882e-05
#> 377     150 1.934134e-04 1.666127e-05
#> 378     150 1.843593e-04 1.960807e-05
#> 379     150 1.669163e-04 2.291828e-05
#> 380     150 1.681233e-04 1.785072e-05
#> 381     150 1.646755e-04 1.905143e-05
#> 382     150 1.873625e-04 1.833666e-05
#> 383     150 1.584607e-04 1.909523e-05
#> 384     150 1.494487e-04 2.744092e-05
#> 385     150 1.539640e-04 1.900844e-05
#> 386     150 1.604525e-04 2.128285e-05
#> 387     150 1.571142e-04 2.582744e-05
#> 388     150 1.467875e-04 2.025071e-05
#> 389     150 1.817136e-04 2.586465e-05
#> 390     150 1.601060e-04 3.486711e-05
#> 391     150 1.527500e-04 2.630914e-05
#> 392     150 1.825167e-04 1.531044e-05
#> 393     150 1.505055e-04 1.991110e-05
#> 394     150 1.837598e-04 1.644515e-05
#> 395     150 1.614936e-04 4.544858e-05
#> 396     150 1.555598e-04 2.474531e-05
#> 397     150 1.751303e-04 1.776269e-05
#> 398     150 1.651111e-04 2.589872e-05
#> 399     150 1.637869e-04 3.064272e-05
#> 400     150 1.694460e-04 1.520250e-05
#> 401     150 1.738804e-04 1.685538e-05
#> 402     150 1.722240e-04 3.123791e-05
#> 403     150 1.471121e-04 1.658194e-05
#> 404     150 1.687523e-04 1.959535e-05
#> 405     150 2.945914e-04 2.241390e-05
#> 406     150 1.631506e-04 1.788922e-05
#> 407     150 1.674377e-04 2.819270e-05
#> 408     150 1.661221e-04 1.755670e-05
#> 409     150 1.504291e-04 1.841208e-05
#> 410     150 1.501874e-04 8.222209e-05
#> 411     150 1.514275e-04 1.984766e-05
#> 412     150 1.760528e-04 2.726822e-04
#> 413     150 1.593024e-04 2.251545e-05
#> 414     150 1.568605e-04 3.848548e-05
#> 415     150 1.500970e-04 2.020868e-05
#> 416     150 1.655930e-04 1.724270e-05
#> 417     150 1.579336e-04 2.109224e-05
#> 418     150 1.709845e-04 2.863219e-05
#> 419     150 1.592823e-04 1.923521e-05
#> 420     150 1.713582e-04 1.710489e-05
#> 421     150 1.427187e-04 1.657991e-05
#> 422     150 1.569930e-04 2.155653e-05
#> 423     150 1.501364e-04 1.433629e-05
#> 424     150 1.805139e-04 2.670733e-05
#> 425     150 1.535945e-04 1.788749e-05
#> 426     150 1.665864e-04 1.584645e-05
#> 427     150 1.753771e-04 2.525791e-05
#> 428     150 1.636859e-04 3.648080e-05
#> 429     150 1.621157e-04 2.506888e-05
#> 430     150 1.577081e-04 1.975706e-05
#> 431     150 1.904207e-04 1.974486e-05
#> 432     150 1.756263e-04 1.669109e-05
#> 433     150 1.425341e-04 1.946775e-05
#> 434     150 1.520832e-04 1.785669e-05
#> 435     150 1.495511e-04 3.419877e-05
#> 436     150 1.705446e-04 1.683549e-05
#> 437     150 1.612790e-04 1.598924e-05
#> 438     150 1.749557e-04 1.723506e-05
#> 439     150 1.420741e-04 2.662785e-05
#> 440     150 1.509071e-04 2.191179e-05
#> 441     150 1.645254e-04 1.650869e-05
#> 442     150 1.546173e-04 3.886625e-05
#> 443     150 1.710297e-04 1.777615e-05
#> 444     150 1.629038e-04 1.883680e-05
#> 445     150 1.629410e-04 3.222050e-05
#> 446     150 1.532248e-04 2.501029e-05
#> 447     150 1.405828e-04 1.555058e-05
#> 448     150 1.786101e-04 4.265987e-05
#> 449     150 1.475729e-04 1.816474e-05
#> 450     150 1.562633e-04 2.517472e-05
#> 451     150 1.546718e-04 2.178772e-05
```

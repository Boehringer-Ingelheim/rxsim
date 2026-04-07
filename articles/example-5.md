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
#>   replicate timepoint analysis sample_size    allocation n_arms alpha delta
#> 1         1  143.0045    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 2         2  148.2639    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#> 3         3  155.7178    final         150 1, 1, 1, 1, 1      5  0.05   0.1
#>   n_total   dunn_min_p mcpmod_min_p
#> 1     150 2.614327e-01 9.360871e-02
#> 2     150 8.302926e-05 3.381185e-06
#> 3     150 1.502310e-04 3.383530e-05
```

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

## Scenario

A **multi-arm, fixed design** (placebo + 4 active doses) with a **single
continuous endpoint**.  
At the **final analysis**, we perform:

- **Dunnett test**: all active doses vs placebo using a normal-theory
  linear model.
- **MCP-Mod**: model-based multiple contrast test + model fitting on a
  candidate set of dose–response shapes.

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

Create a `Timer` with discrete time points to reach the target sample
size.

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

Define population generators for each dose level using the Emax model.

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
#>   replicate sample_size    allocation n_arms alpha delta n_total   dunn_min_p
#> 1         1         150 1, 1, 1, 1, 1      5  0.05   0.1     150 2.614327e-01
#> 2         2         150 1, 1, 1, 1, 1      5  0.05   0.1     150 8.981219e-05
#> 3         3         150 1, 1, 1, 1, 1      5  0.05   0.1     150 1.540937e-04
#>   mcpmod_min_p
#> 1 9.360871e-02
#> 2 1.123651e-05
#> 3 1.831518e-05
```

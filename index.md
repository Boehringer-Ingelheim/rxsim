# rxsim

> **Reduce friction in randomized controlled trial simulation.** rxsim
> lets you prototype Phase II/III trial designs at high speed — without
> writing hundreds of lines of low-level plumbing code.

------------------------------------------------------------------------

## Why rxsim?

Simulating a randomized controlled trial from scratch means manually
wiring together enrollment schedules, randomization, data snapshots,
analyses, and result collection — across potentially thousands of
replicates. Done in base R or even with the tidyverse, that
orchestration routinely reaches **1,000+ lines of fragile,
hard-to-reuse, -maintain, and -review code**.

**rxsim replaces that scaffolding with \<100 lines of high-level,
declarative code**, so you can stay focused on the design decisions that
matter: endpoints, allocation, analysis rules, and adaptive strategies.

## Key Features

- 👌 **Very easy to use** — Add/remove arms, and interims as list
  elements
- 🏗️ **Declarative trial design** — define arms, populations, and
  analyses as plain R functions; let rxsim handle the orchestration
- 📅 **Flexible enrollment & dropout** — stochastic or deterministic
  schedules
- ⚡ **Trigger-based analyses** — fire any analysis at a calendar time,
  an enrollment fraction, or any custom data condition
- 🔁 **Scalable replication** — generate and run thousands of
  independent trial replicates
- 🧩 **Any endpoint type** — continuous, time-to-event, binary,
  correlated multi-endpoint; supply any analysis function
- 🔬 **Adaptive designs** — interims, Go/No-Go rules, Bayesian
  borrowing; all composable from the same building blocks
- 📦 **Integrates cleanly** — works with `survival`, `RBesT`,
  `DoseFinding`, `multcomp`, `dplyr`, and more

## Core Concepts

rxsim is built on three composable classes:

| Class            | Role                                                                                                              |
|------------------|-------------------------------------------------------------------------------------------------------------------|
| **`Population`** | Holds subject-level data (endpoints, covariates) and tracks each subject’s enrollment and dropout times           |
| **`Timer`**      | Drives the trial clock — stores discrete timepoints per arm and evaluates condition-triggered analyses            |
| **`Trial`**      | Orchestrates the simulation — iterates over timepoints, updates populations, snapshots data, and collects results |

The high-level helpers
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
and
[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md)
wrap these three classes into an ergonomic one-call workflow for the
common case.

## Quick Example

A two-arm, fixed-design trial with a continuous endpoint in under 30
lines:

``` r
library(rxsim)

# 1. Define populations — one generator function per arm
population_generators <- list(
  placebo   = function(n) data.frame(id = 1:n, y = rnorm(n, 0.0), readout_time = 1),
  treatment = function(n) data.frame(id = 1:n, y = rnorm(n, 0.3), readout_time = 1)
)

# 2. Define what triggers an analysis and what to compute
sample_size <- 100L
analysis_generators <- list(
  final = list(
    trigger  = rlang::exprs(sum(!is.na(enroll_time)) >= !!sample_size),
    analysis = function(df, timer) {
      enrolled <- subset(df, !is.na(enroll_time))
      data.frame(n = nrow(enrolled), p_value = t.test(y ~ arm, data = enrolled)$p.value)
    }
  )
)

# 3. Create 500 independent replicates and run them all
trials <- replicate_trial(
  trial_name            = "phase2_example",
  sample_size           = sample_size,
  arms                  = c("placebo", "treatment"),
  allocation            = c(1, 1),
  enrollment            = function(n) rexp(n, rate = 1),
  dropout               = function(n) rexp(n, rate = 0.01),
  analysis_generators   = analysis_generators,
  population_generators = population_generators,
  n                     = 500
)

run_trials(trials)
```

Collect results across all replicates in one pass — see the [Getting
Started
vignette](https://boehringer-ingelheim.github.io/rxsim/articles/getting-started.html)
for the full walkthrough.

## Documentation

Full documentation lives on the **[rxsim package
website](https://boehringer-ingelheim.github.io/rxsim)**.

| Vignette                                                                                                        | What it covers                                                                                                                                                                                                          |
|-----------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [Getting Started](https://boehringer-ingelheim.github.io/rxsim/articles/getting-started.html)                   | End-to-end walkthrough of a complete trial simulation                                                                                                                                                                   |
| [Core Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.html)                            | Population, Timer, and Trial in depth; trigger expressions explained                                                                                                                                                    |
| [Enrollment & Dropout Modeling](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment-dropout.html)  | [`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md) vs [`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md) — stochastic vs piecewise-constant |
| [Example 1: Fixed design, continuous](https://boehringer-ingelheim.github.io/rxsim/articles/example-1.html)     | Two-arm fixed design, t-test                                                                                                                                                                                            |
| [Example 2: Interim analysis, continuous](https://boehringer-ingelheim.github.io/rxsim/articles/example-2.html) | Interim + final analyses                                                                                                                                                                                                |
| [Example 3: Correlated endpoints](https://boehringer-ingelheim.github.io/rxsim/articles/example-3.html)         | Two correlated continuous endpoints, Holm correction                                                                                                                                                                    |
| [Example 4: Continuous + time-to-event](https://boehringer-ingelheim.github.io/rxsim/articles/example-4.html)   | Mixed endpoint — continuous + TTE, Cox PH                                                                                                                                                                               |
| [Example 5: Multi-arm, Dunnett + MCP-Mod](https://boehringer-ingelheim.github.io/rxsim/articles/example-5.html) | Dose-finding with Dunnett test and MCP-Mod                                                                                                                                                                              |
| [Example 6: Subgroup analysis](https://boehringer-ingelheim.github.io/rxsim/articles/example-6.html)            | Overall + exploratory subgroup treatment effects                                                                                                                                                                        |
| [Example 7: Bayesian Go/No-Go](https://boehringer-ingelheim.github.io/rxsim/articles/example-7.html)            | Bayesian decision rule with historical placebo borrowing (RBesT)                                                                                                                                                        |

## Installation

Install the development version from
[GitHub](https://github.com/Boehringer-Ingelheim/rxsim):

``` r
install.packages("pak")
pak::pak("Boehringer-Ingelheim/rxsim")
```

------------------------------------------------------------------------

## FAQ

**Why use rxsim instead of writing simulations in base R?**

  

Building a clinical trial simulation from scratch typically balloons to
**1,000+ lines of low-level code**. Even with the tidyverse, the
orchestration is entirely on you: you write the clock loop, you decide
when analyses fire, you manage state across replicates, and handle
result aggregation.

rxsim provides that orchestration layer out of the box. The same
simulation fits in **\<100 lines of high-level, readable code** where
each section has a clear purpose: define trial scenario, populations,
analysis rules, replicate, and run. This dramatically reduces the
cognitive burden and makes it faster to: - prototype a new trial
design - swap out an endpoint model - add/remove arms - add/remove looks
or analyses - hand the code to a colleague

**Is rxsim only useful for simple two-arm trials?**

  

No. rxsim supports any number of arms, any endpoint type, and arbitrary
analysis functions. The built-in examples cover continuous endpoints,
correlated multi-endpoints, time-to-event endpoints, multi-arm
dose-finding (Dunnett + MCP-Mod), subgroup analyses, and Bayesian
Go/No-Go rules with historical borrowing. If your analysis can be
written as an R function, rxsim can run it.

**Can rxsim handle adaptive designs with interim analyses?**

  

Yes. Trigger conditions fire whenever a boolean expression on the
accumulating data is satisfied — at a calendar time, an enrollment
fraction, or any custom condition. Multiple triggers can be registered
on the same trial, allowing interim and final analyses to coexist
naturally. Adaptive rules (e.g., stopping early, adjusting sample size)
can be implemented inside analysis functions and propagated across the
trial state.

See [Core
Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.html)
for a detailed explanation of the trigger system.

**Is rxsim appropriate for regulatory submissions?**

  

rxsim is designed for **rapid trial design prototyping and operating
characteristic evaluation** — estimating power, type I error, and
decision probabilities under a range of scenarios. It is not a validated
submission package, yet. Simulation outputs it produces can feed into
submission-supporting analyses when the user applies appropriate
validation and documentation.

**How do I scale to thousands of replicates efficiently?**

  

`replicate_trial(n = N)` creates `N` fully independent `Trial` objects
(each with its own population data and timer state).
[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md)
then executes them sequentially. For very large `N`, parallelization is
straightforward — replace
[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md)
with any R parallel back-end:

``` r
# Example using parallel::mclapply (Unix/macOS)
parallel::mclapply(trials, function(tr) tr$run(), mc.cores = parallel::detectCores())
```

**What is the difference between `gen_plan()` and `gen_timepoints()`?**

  

Both functions produce an enrollment/dropout schedule that feeds into a
`Timer`, but they model timing differently:

- **[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md)**
  is **stochastic**: you supply a inter-event time probability
  distribution function (e.g., `function(n) rexp(n, rate = 1)`). Each
  replicate gets a different realization of the enrollment process,
  capturing natural variability in study timelines.
- **[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md)**
  is **deterministic piecewise-constant rate**: you specify enrollment
  rates per time period (e.g., 5 patients/month for months 1–4, then
  10/month). Every replicate follows the same fixed schedule.

See the [Enrollment & Dropout Modeling
vignette](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment-dropout.html)
for a side-by-side comparison.

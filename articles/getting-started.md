# Getting Started

This vignette walks through a complete rxsim simulation from start to
finish: defining a trial scenario, building arm populations, registering
an analysis trigger, running replicates, and collecting results. By the
end you will have seen the full rxsim workflow in one place and will be
ready to adapt it to your own design. For deeper explanations of each
building block, see the [Core
Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.md)
vignette.

## Load the package

``` r
library(rxsim)
```

## Step 1 — Define the scenario

The `scenario` tibble captures the design parameters you want to track
in your results — sample size, allocation ratios, or any other factors
you might vary across simulation runs. Wrapping them in
[`tidyr::expand_grid()`](https://tidyr.tidyverse.org/reference/expand_grid.html)
produces a one-row data frame that gets embedded in every analysis
result. When you later stack results across many scenarios or
sensitivity analyses, this metadata keeps each row traceable back to its
design.

``` r
sample_size <- 30
arms        <- c("control", "treatment")
allocation  <- c(1, 1)
true_delta  <- 0.5

scenario <- tidyr::expand_grid(
  sample_size = sample_size,
  allocation  = list(allocation),
  true_delta  = true_delta
)
```

## Step 2 — Define arm populations

Each trial arm is represented by a **generator function** — a plain R
function that takes `n` (the number of subjects planned for that arm)
and returns a data frame of subject-level endpoint data. rxsim calls
these functions once per replicate to draw a fresh population, so every
replicate gets independent data.

The data frame must contain at least three columns:

- `id` — a unique integer identifier per subject
- `readout_time` — the time (in trial-time units) after enrollment when
  the endpoint is observed (use `0` for baseline or
  immediately-available data)
- at least one endpoint column

[`vector_to_dataframe()`](https://boehringer-ingelheim.github.io/rxsim/reference/vector_to_dataframe.md)
is a convenience helper that wraps a numeric vector into this standard
format with `id`, `data`, and `readout_time = 0`.

Here we simulate a simple continuous endpoint: the control arm is
standard-normal, the treatment arm has a mean shift of `true_delta`.

``` r
population_generators <- list(
  control   = function(n) vector_to_dataframe(rnorm(n)),
  treatment = function(n) vector_to_dataframe(rnorm(n, mean = true_delta))
)
```

## Step 3 — Define enrollment and dropout

Enrollment and dropout are modelled as random processes. Each function
takes `n` and returns a vector of **inter-event times** — the waiting
times between successive enrollments (or dropouts). These times are
drawn independently for every replicate, giving each trial its own
enrollment trajectory.

`rexp(n, rate = 1)` generates exponentially-distributed inter-arrival
times with a mean gap of 1 time unit between enrolments — a common
approximation for Poisson arrivals. A lower dropout rate (`0.05`)
reflects a trial where most subjects complete the study.

``` r
enrollment <- function(n) rexp(n, rate = 1.0)
dropout    <- function(n) rexp(n, rate = 0.05)
```

## Step 4 — Define analysis triggers

Analysis triggers are the heart of rxsim. Each trigger pairs a
**condition** (when to fire) with an **analysis function** (what to
compute when it fires).

The condition is written as a `dplyr`-style boolean expression inside
[`rlang::exprs()`](https://rlang.r-lib.org/reference/defusing-advanced.html).
It is evaluated against a snapshot of all currently enrolled subjects at
each timepoint. Here the condition fires once full enrollment is
reached, i.e., `sample_size` subjects have accumulated a non-`NA`
`enroll_time`.

The `!!` operator (pronounced “bang-bang”) **injects the current value**
of `sample_size` into the expression at definition time, rather than
looking it up at evaluation time. This is necessary because the
expression is stored and evaluated later inside the simulation loop.

When the condition is met, rxsim calls the analysis function with two
arguments:

- `df`: a data frame snapshot of all enrolled subjects at the triggering
  timepoint, with columns from the population data plus `enroll_time`,
  `drop_time`, and `arm`
- `time`: the current trial clock time at which the trigger fired

The function should return a data frame (one row per trigger event is
the conventional pattern).

``` r
analysis_generators <- list(
  final = list(
    trigger = rlang::exprs(
      sum(!is.na(enroll_time)) >= !!sample_size
    ),
    analysis = function(df, timer) {
      enrolled <- subset(df, !is.na(enroll_time))
      data.frame(
        scenario,
        n_enrolled  = nrow(enrolled),
        mean_ctrl   = mean(enrolled$data[enrolled$arm == "control"]),
        mean_trt    = mean(enrolled$data[enrolled$arm == "treatment"]),
        stringsAsFactors = FALSE
      )
    }
  )
)
```

## Step 5 — Create replicates and run

[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md)
generates `n` fully independent `Trial` objects. For each replicate it:

1.  Samples a fresh enrollment/dropout timeline from your functions
2.  Calls each population generator to draw new subject-level data
3.  Registers your analysis triggers on the trial’s internal timer

[`run_trials()`](https://boehringer-ingelheim.github.io/rxsim/reference/run_trials.md)
then executes every replicate’s simulation loop in sequence.

``` r
set.seed(42)

trials <- replicate_trial(
  trial_name            = "getting_started",
  sample_size           = sample_size,
  arms                  = arms,
  allocation            = allocation,
  enrollment            = enrollment,
  dropout               = dropout,
  analysis_generators   = analysis_generators,
  population_generators = population_generators,
  n                     = 5
)

run_trials(trials)
#> [[1]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: getting_started_1
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
#>     name: getting_started_2
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
#>     name: getting_started_3
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[4]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: getting_started_4
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
#> 
#> [[5]]
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: getting_started_5
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6
```

## Step 6 — Collect and interpret results

### Analysis results

After running, each `Trial` object exposes a `results` list. It is
indexed first by the **timepoint** at which an analysis fired (e.g.,
`"time_30.4"`), then by the **analysis name** you gave it (here
`"final"`). The value is whatever your analysis function returned.

The
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
helper gathers every analysis from every timepoint across all replicates
into one tidy data frame:

``` r
replicate_results <- collect_results(trials)
replicate_results
#>     replicate timepoint analysis sample_size allocation true_delta n_enrolled
#> 1           1  29.71727    final          30       1, 1        0.5         30
#> 2           1  71.75884    final          30       1, 1        0.5         30
#> 3           1  81.71066    final          30       1, 1        0.5         30
#> 4           1  88.80173    final          30       1, 1        0.5         30
#> 5           1 123.92991    final          30       1, 1        0.5         30
#> 6           1 138.67989    final          30       1, 1        0.5         30
#> 7           1 139.25641    final          30       1, 1        0.5         30
#> 8           1 146.28070    final          30       1, 1        0.5         30
#> 9           1 154.13389    final          30       1, 1        0.5         30
#> 10          1 192.17561    final          30       1, 1        0.5         30
#> 11          1 207.74148    final          30       1, 1        0.5         30
#> 12          1 220.86782    final          30       1, 1        0.5         30
#> 13          1 228.59602    final          30       1, 1        0.5         30
#> 14          1 256.91786    final          30       1, 1        0.5         30
#> 15          1 287.38043    final          30       1, 1        0.5         30
#> 16          1 296.15466    final          30       1, 1        0.5         30
#> 17          1 379.51726    final          30       1, 1        0.5         30
#> 18          1 403.41940    final          30       1, 1        0.5         30
#> 19          1 403.99571    final          30       1, 1        0.5         30
#> 20          1 540.92715    final          30       1, 1        0.5         30
#> 21          1 544.19131    final          30       1, 1        0.5         30
#> 22          1 577.18203    final          30       1, 1        0.5         30
#> 23          1 599.76724    final          30       1, 1        0.5         30
#> 24          1 605.59251    final          30       1, 1        0.5         30
#> 25          1 616.62545    final          30       1, 1        0.5         30
#> 26          1 619.17132    final          30       1, 1        0.5         30
#> 27          1 648.08404    final          30       1, 1        0.5         30
#> 28          1 666.36446    final          30       1, 1        0.5         30
#> 29          1 673.06152    final          30       1, 1        0.5         30
#> 30          1 842.30823    final          30       1, 1        0.5         30
#> 31          2  27.31408    final          30       1, 1        0.5         30
#> 32          2  48.00359    final          30       1, 1        0.5         30
#> 33          2 109.71186    final          30       1, 1        0.5         30
#> 34          2 131.84283    final          30       1, 1        0.5         30
#> 35          2 133.73988    final          30       1, 1        0.5         30
#> 36          2 136.11920    final          30       1, 1        0.5         30
#> 37          2 137.84044    final          30       1, 1        0.5         30
#> 38          2 142.18902    final          30       1, 1        0.5         30
#> 39          2 160.61632    final          30       1, 1        0.5         30
#> 40          2 198.01478    final          30       1, 1        0.5         30
#> 41          2 264.40849    final          30       1, 1        0.5         30
#> 42          2 294.33890    final          30       1, 1        0.5         30
#> 43          2 366.11547    final          30       1, 1        0.5         30
#> 44          2 367.31573    final          30       1, 1        0.5         30
#> 45          2 394.88895    final          30       1, 1        0.5         30
#> 46          2 438.12156    final          30       1, 1        0.5         30
#> 47          2 440.57908    final          30       1, 1        0.5         30
#> 48          2 484.79892    final          30       1, 1        0.5         30
#> 49          2 526.34752    final          30       1, 1        0.5         30
#> 50          2 528.33233    final          30       1, 1        0.5         30
#> 51          2 544.26888    final          30       1, 1        0.5         30
#> 52          2 560.74408    final          30       1, 1        0.5         30
#> 53          2 561.97384    final          30       1, 1        0.5         30
#> 54          2 638.74485    final          30       1, 1        0.5         30
#> 55          2 650.69527    final          30       1, 1        0.5         30
#> 56          2 699.76729    final          30       1, 1        0.5         30
#> 57          2 704.52535    final          30       1, 1        0.5         30
#> 58          3  35.29848    final          30       1, 1        0.5         30
#> 59          3  69.62658    final          30       1, 1        0.5         30
#> 60          3  70.81995    final          30       1, 1        0.5         30
#> 61          3  78.09509    final          30       1, 1        0.5         30
#> 62          3 112.62422    final          30       1, 1        0.5         30
#> 63          3 120.26618    final          30       1, 1        0.5         30
#> 64          3 166.04672    final          30       1, 1        0.5         30
#> 65          3 244.41745    final          30       1, 1        0.5         30
#> 66          3 244.94703    final          30       1, 1        0.5         30
#> 67          3 250.17608    final          30       1, 1        0.5         30
#> 68          3 277.54075    final          30       1, 1        0.5         30
#> 69          3 279.03777    final          30       1, 1        0.5         30
#> 70          3 294.09805    final          30       1, 1        0.5         30
#> 71          3 342.65569    final          30       1, 1        0.5         30
#> 72          3 367.39120    final          30       1, 1        0.5         30
#> 73          3 411.00304    final          30       1, 1        0.5         30
#> 74          3 443.18313    final          30       1, 1        0.5         30
#> 75          3 457.78388    final          30       1, 1        0.5         30
#> 76          3 468.89933    final          30       1, 1        0.5         30
#> 77          3 483.59472    final          30       1, 1        0.5         30
#> 78          3 508.70768    final          30       1, 1        0.5         30
#> 79          3 510.42419    final          30       1, 1        0.5         30
#> 80          3 555.17049    final          30       1, 1        0.5         30
#> 81          3 573.74553    final          30       1, 1        0.5         30
#> 82          3 574.31182    final          30       1, 1        0.5         30
#> 83          3 574.90874    final          30       1, 1        0.5         30
#> 84          3 576.72932    final          30       1, 1        0.5         30
#> 85          3 591.75515    final          30       1, 1        0.5         30
#> 86          3 662.12355    final          30       1, 1        0.5         30
#> 87          4  23.61616    final          30       1, 1        0.5         30
#> 88          4  29.72873    final          30       1, 1        0.5         30
#> 89          4  72.60486    final          30       1, 1        0.5         30
#> 90          4  88.63264    final          30       1, 1        0.5         30
#> 91          4  99.34682    final          30       1, 1        0.5         30
#> 92          4 139.98048    final          30       1, 1        0.5         30
#> 93          4 206.27598    final          30       1, 1        0.5         30
#> 94          4 266.56996    final          30       1, 1        0.5         30
#> 95          4 283.19883    final          30       1, 1        0.5         30
#> 96          4 304.96656    final          30       1, 1        0.5         30
#> 97          4 314.46469    final          30       1, 1        0.5         30
#> 98          4 328.46072    final          30       1, 1        0.5         30
#> 99          4 329.15820    final          30       1, 1        0.5         30
#> 100         4 339.53599    final          30       1, 1        0.5         30
#> 101         4 344.97937    final          30       1, 1        0.5         30
#> 102         4 385.33577    final          30       1, 1        0.5         30
#> 103         4 385.34138    final          30       1, 1        0.5         30
#> 104         4 456.78964    final          30       1, 1        0.5         30
#> 105         4 534.01374    final          30       1, 1        0.5         30
#> 106         4 570.68747    final          30       1, 1        0.5         30
#> 107         4 611.73903    final          30       1, 1        0.5         30
#> 108         4 639.27803    final          30       1, 1        0.5         30
#> 109         4 657.52944    final          30       1, 1        0.5         30
#> 110         4 670.87695    final          30       1, 1        0.5         30
#> 111         4 695.91081    final          30       1, 1        0.5         30
#> 112         4 710.49484    final          30       1, 1        0.5         30
#> 113         4 728.27773    final          30       1, 1        0.5         30
#> 114         4 739.33553    final          30       1, 1        0.5         30
#> 115         4 757.50669    final          30       1, 1        0.5         30
#> 116         4 758.12874    final          30       1, 1        0.5         30
#> 117         4 774.26422    final          30       1, 1        0.5         30
#> 118         5  25.83584    final          30       1, 1        0.5         30
#> 119         5  42.69936    final          30       1, 1        0.5         30
#> 120         5  55.30675    final          30       1, 1        0.5         30
#> 121         5  78.63514    final          30       1, 1        0.5         30
#> 122         5  98.58472    final          30       1, 1        0.5         30
#> 123         5 124.76022    final          30       1, 1        0.5         30
#> 124         5 137.43827    final          30       1, 1        0.5         30
#> 125         5 185.54669    final          30       1, 1        0.5         30
#> 126         5 204.11105    final          30       1, 1        0.5         30
#> 127         5 209.57048    final          30       1, 1        0.5         30
#> 128         5 217.52736    final          30       1, 1        0.5         30
#> 129         5 218.93812    final          30       1, 1        0.5         30
#> 130         5 224.62617    final          30       1, 1        0.5         30
#> 131         5 292.09763    final          30       1, 1        0.5         30
#> 132         5 303.22840    final          30       1, 1        0.5         30
#> 133         5 306.62372    final          30       1, 1        0.5         30
#> 134         5 356.25605    final          30       1, 1        0.5         30
#> 135         5 381.40767    final          30       1, 1        0.5         30
#> 136         5 418.34338    final          30       1, 1        0.5         30
#> 137         5 447.44924    final          30       1, 1        0.5         30
#> 138         5 461.90921    final          30       1, 1        0.5         30
#> 139         5 464.03910    final          30       1, 1        0.5         30
#> 140         5 495.45607    final          30       1, 1        0.5         30
#> 141         5 508.54993    final          30       1, 1        0.5         30
#> 142         5 516.79836    final          30       1, 1        0.5         30
#> 143         5 547.04878    final          30       1, 1        0.5         30
#> 144         5 633.22196    final          30       1, 1        0.5         30
#> 145         5 640.68178    final          30       1, 1        0.5         30
#> 146         5 647.93051    final          30       1, 1        0.5         30
#> 147         5 653.59243    final          30       1, 1        0.5         30
#>      mean_ctrl   mean_trt
#> 1   -0.1657292 0.61748974
#> 2   -0.1657292 0.61748974
#> 3   -0.1657292 0.61748974
#> 4   -0.1657292 0.61748974
#> 5   -0.1657292 0.61748974
#> 6   -0.1657292 0.61748974
#> 7   -0.1657292 0.61748974
#> 8   -0.1657292 0.61748974
#> 9   -0.1657292 0.61748974
#> 10  -0.1657292 0.61748974
#> 11  -0.1657292 0.61748974
#> 12  -0.1657292 0.61748974
#> 13  -0.1657292 0.61748974
#> 14  -0.1657292 0.61748974
#> 15  -0.1657292 0.61748974
#> 16  -0.1657292 0.61748974
#> 17  -0.1657292 0.61748974
#> 18  -0.1657292 0.61748974
#> 19  -0.1657292 0.61748974
#> 20  -0.1657292 0.61748974
#> 21  -0.1657292 0.61748974
#> 22  -0.1657292 0.61748974
#> 23  -0.1657292 0.61748974
#> 24  -0.1657292 0.61748974
#> 25  -0.1657292 0.61748974
#> 26  -0.1657292 0.61748974
#> 27  -0.1657292 0.61748974
#> 28  -0.1657292 0.61748974
#> 29  -0.1657292 0.61748974
#> 30  -0.1657292 0.61748974
#> 31  -0.0978443 0.64303730
#> 32  -0.0978443 0.64303730
#> 33  -0.0978443 0.64303730
#> 34  -0.0978443 0.64303730
#> 35  -0.0978443 0.64303730
#> 36  -0.0978443 0.64303730
#> 37  -0.0978443 0.64303730
#> 38  -0.0978443 0.64303730
#> 39  -0.0978443 0.64303730
#> 40  -0.0978443 0.64303730
#> 41  -0.0978443 0.64303730
#> 42  -0.0978443 0.64303730
#> 43  -0.0978443 0.64303730
#> 44  -0.0978443 0.64303730
#> 45  -0.0978443 0.64303730
#> 46  -0.0978443 0.64303730
#> 47  -0.0978443 0.64303730
#> 48  -0.0978443 0.64303730
#> 49  -0.0978443 0.64303730
#> 50  -0.0978443 0.64303730
#> 51  -0.0978443 0.64303730
#> 52  -0.0978443 0.64303730
#> 53  -0.0978443 0.64303730
#> 54  -0.0978443 0.64303730
#> 55  -0.0978443 0.64303730
#> 56  -0.0978443 0.64303730
#> 57  -0.0978443 0.64303730
#> 58  -0.2573593 0.07050019
#> 59  -0.2573593 0.07050019
#> 60  -0.2573593 0.07050019
#> 61  -0.2573593 0.07050019
#> 62  -0.2573593 0.07050019
#> 63  -0.2573593 0.07050019
#> 64  -0.2573593 0.07050019
#> 65  -0.2573593 0.07050019
#> 66  -0.2573593 0.07050019
#> 67  -0.2573593 0.07050019
#> 68  -0.2573593 0.07050019
#> 69  -0.2573593 0.07050019
#> 70  -0.2573593 0.07050019
#> 71  -0.2573593 0.07050019
#> 72  -0.2573593 0.07050019
#> 73  -0.2573593 0.07050019
#> 74  -0.2573593 0.07050019
#> 75  -0.2573593 0.07050019
#> 76  -0.2573593 0.07050019
#> 77  -0.2573593 0.07050019
#> 78  -0.2573593 0.07050019
#> 79  -0.2573593 0.07050019
#> 80  -0.2573593 0.07050019
#> 81  -0.2573593 0.07050019
#> 82  -0.2573593 0.07050019
#> 83  -0.2573593 0.07050019
#> 84  -0.2573593 0.07050019
#> 85  -0.2573593 0.07050019
#> 86  -0.2573593 0.07050019
#> 87  -0.4199802 0.25120329
#> 88  -0.4199802 0.25120329
#> 89  -0.4199802 0.25120329
#> 90  -0.4199802 0.25120329
#> 91  -0.4199802 0.25120329
#> 92  -0.4199802 0.25120329
#> 93  -0.4199802 0.25120329
#> 94  -0.4199802 0.25120329
#> 95  -0.4199802 0.25120329
#> 96  -0.4199802 0.25120329
#> 97  -0.4199802 0.25120329
#> 98  -0.4199802 0.25120329
#> 99  -0.4199802 0.25120329
#> 100 -0.4199802 0.25120329
#> 101 -0.4199802 0.25120329
#> 102 -0.4199802 0.25120329
#> 103 -0.4199802 0.25120329
#> 104 -0.4199802 0.25120329
#> 105 -0.4199802 0.25120329
#> 106 -0.4199802 0.25120329
#> 107 -0.4199802 0.25120329
#> 108 -0.4199802 0.25120329
#> 109 -0.4199802 0.25120329
#> 110 -0.4199802 0.25120329
#> 111 -0.4199802 0.25120329
#> 112 -0.4199802 0.25120329
#> 113 -0.4199802 0.25120329
#> 114 -0.4199802 0.25120329
#> 115 -0.4199802 0.25120329
#> 116 -0.4199802 0.25120329
#> 117 -0.4199802 0.25120329
#> 118 -0.4170091 0.59879107
#> 119 -0.4170091 0.59879107
#> 120 -0.4170091 0.59879107
#> 121 -0.4170091 0.59879107
#> 122 -0.4170091 0.59879107
#> 123 -0.4170091 0.59879107
#> 124 -0.4170091 0.59879107
#> 125 -0.4170091 0.59879107
#> 126 -0.4170091 0.59879107
#> 127 -0.4170091 0.59879107
#> 128 -0.4170091 0.59879107
#> 129 -0.4170091 0.59879107
#> 130 -0.4170091 0.59879107
#> 131 -0.4170091 0.59879107
#> 132 -0.4170091 0.59879107
#> 133 -0.4170091 0.59879107
#> 134 -0.4170091 0.59879107
#> 135 -0.4170091 0.59879107
#> 136 -0.4170091 0.59879107
#> 137 -0.4170091 0.59879107
#> 138 -0.4170091 0.59879107
#> 139 -0.4170091 0.59879107
#> 140 -0.4170091 0.59879107
#> 141 -0.4170091 0.59879107
#> 142 -0.4170091 0.59879107
#> 143 -0.4170091 0.59879107
#> 144 -0.4170091 0.59879107
#> 145 -0.4170091 0.59879107
#> 146 -0.4170091 0.59879107
#> 147 -0.4170091 0.59879107
```

Each row is one replicate. The `mean_ctrl` and `mean_trt` columns show
the arm-level sample means at the moment the final analysis fired. The
`timepoint` and `analysis` columns identify when and which analysis
produced each row — essential when a trial has both interim and final
analyses. Variation across replicates reflects both the stochastic
endpoint data and the randomness in enrollment timing; exactly what
operating characteristic simulations are designed to characterise.

### Locked data snapshots

Every time an analysis fires, rxsim also saves a **locked data
snapshot**, the full subject-level data frame at that timepoint. You can
inspect it directly to debug your analysis function, audit which
subjects were enrolled, or compute additional statistics post-hoc.

``` r
# One snapshot per timepoint that fired in replicate 1
names(trials[[1]]$locked_data)
#>  [1] "time_29.7172711929598" "time_71.7588397920052" "time_81.7106552388439"
#>  [4] "time_88.8017284464011" "time_123.929908531333" "time_138.679889248391"
#>  [7] "time_139.256406634753" "time_146.280697615524" "time_154.133894694986"
#> [10] "time_192.175611581853" "time_207.741484369442" "time_220.867823778942"
#> [13] "time_228.5960165945"   "time_256.917855072035" "time_287.380430737777"
#> [16] "time_296.154664244769" "time_379.517257070982" "time_403.419397851426"
#> [19] "time_403.995706183809" "time_540.927149055388" "time_544.191309155759"
#> [22] "time_577.182029686913" "time_599.767237766782" "time_605.592512902939"
#> [25] "time_616.625447408432" "time_619.171321070889" "time_648.084043521793"
#> [28] "time_666.364456690367" "time_673.061517277252" "time_842.308234163735"

# First six rows of the snapshot
head(trials[[1]]$locked_data[[1]])
#>   id       data readout_time     arm enroll_time drop_time subject_id
#> 1  1 -2.6882473            0 control   0.8592321        NA          1
#> 2  2  0.8666502            0 control  11.6939492        NA          2
#> 3  3  0.1687397            0 control  26.6106187        NA          3
#> 4  4 -1.0908241            0 control   1.6540916        NA          4
#> 5  5 -0.3803481            0 control   9.5016932        NA          5
#> 6  6 -0.9480918            0 control  13.9319792        NA          6
#>   measurement_time     time
#> 1        0.8592321 29.71727
#> 2       11.6939492 29.71727
#> 3       26.6106187 29.71727
#> 4        1.6540916 29.71727
#> 5        9.5016932 29.71727
#> 6       13.9319792 29.71727
```

The locked data contains the population columns (`id`, `data`, `arm`,
`readout_time`) plus three columns added by rxsim:

| Column        | Meaning                                                           |
|---------------|-------------------------------------------------------------------|
| `enroll_time` | Calendar time the subject was enrolled (`NA` if not yet enrolled) |
| `drop_time`   | Calendar time the subject dropped out (`NA` if still active)      |
| `time`        | The trial clock time at which this snapshot was taken             |

## Next steps

Now that you have seen the full workflow, here is where to go next:

- **[Core
  Concepts](https://boehringer-ingelheim.github.io/rxsim/articles/concepts.md)**
  — understand how `Population`, `Timer`, and `Trial` compose, and how
  to write more advanced trigger expressions
- **[Enrollment & Dropout
  Modeling](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment-dropout.md)**
  — choose between stochastic (`gen_plan`) and piecewise-constant
  (`gen_timepoints`) schedules
- **[Example
  1](https://boehringer-ingelheim.github.io/rxsim/articles/example-1.md)**
  through **[Example
  7](https://boehringer-ingelheim.github.io/rxsim/articles/example-7.md)**
  — progressively complex designs: correlated endpoints, time-to-event,
  multi-arm dose-finding, subgroup analyses, and Bayesian Go/No-Go rules

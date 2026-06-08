# Conditions and Triggers

``` r
library(rxsim)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## What a Condition does

A `Condition` pairs two things:

1.  a trigger (`where`) that decides **when** the condition fires
2.  an analysis function (`analysis`) that decides **what** to compute
    when it fires

You can think of one `Condition` as one planned look. If your design has
two looks (for example, interim + final), define two `Condition`
objects. More looks means more `Condition` instances.

## Trigger primitives

rxsim provides trigger primitives that return `trigger` objects.

``` r
sample_size <- 120

trig_enroll_50  <- enroll_trigger(0.5, sample_size)     # 50% enrolled
trig_enroll_100 <- enroll_trigger(1.0, sample_size)     # full enrollment
trig_calendar   <- calendar_trigger(8)                         # single look
trig_count      <- count_trigger("enroll_time", ">=", 40)       # count-based
trig_value      <- value_trigger("analysis_flag", "==", TRUE)   # value-based

class(trig_enroll_50)
#> [1] "trigger"
class(trig_calendar)
#> [1] "trigger"
```

Use cases:

- [`enroll_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md):
  percentage-based interims tied to accrued sample size
- [`calendar_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md):
  fixed-time looks (e.g., month 4, 8, 12)
- [`count_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  /
  [`value_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md):
  custom state-driven logic

## Composing triggers with `&` and `|`

Triggers are composable:

- `a & b` means both conditions must hold
- `a | b` means either condition can fire the look

``` r
combined_and <- enroll_trigger(0.5, sample_size) & calendar_trigger(6)
combined_or  <- calendar_trigger(12) | enroll_trigger(1.0, sample_size)

class(combined_and)
#> [1] "trigger"
class(combined_or)
#> [1] "trigger"
```

This is useful when an interim should happen only after both enough
enrollment and sufficient follow-up time.

## Registering composed triggers with `Condition$new()`

Helper constructors are great for common single-trigger cases. When your
rule is composed (for example `a & b`), register it directly with
`Condition$new()`.

``` r
cond_composed <- Condition$new(
  where = enroll_trigger(0.5, sample_size) & calendar_trigger(6),
  analysis = function(df, current_time) {
    enrolled <- subset(df, !is.na(enroll_time))
    data.frame(
      time = current_time,
      n = nrow(enrolled)
    )
  },
  name = "interim_enroll_and_time"
)

class(cond_composed)
#> [1] "Condition" "R6"
length(cond_composed$where)
#> [1] 2
```

`where` accepts any `trigger` object (primitive or composed) and
internally converts it to the filter predicates used during condition
evaluation.

## Helper constructors for common conditions

Two helper constructors cover common workflows and return `Condition`
objects directly.

``` r
cond_calendar <- condition_calendar_time(
  cal_time = 8,
  analysis = function(df, current_time) {
    enrolled <- subset(df, !is.na(enroll_time))
    data.frame(time = current_time, n = nrow(enrolled))
  },
  name = "cal_8"
)

cond_fraction <- condition_enrollment_fraction(
  fraction = 0.5,
  sample_size = sample_size,
  analysis = function(df, current_time) {
    enrolled <- subset(df, !is.na(enroll_time))
    data.frame(time = current_time, frac = nrow(enrolled) / sample_size)
  },
  name = "frac_50"
)

class(cond_calendar)
#> [1] "Condition" "R6"
class(cond_fraction)
#> [1] "Condition" "R6"
```

## Interim + final example (two Conditions)

This pattern is the core of multi-look designs: define one interim
condition and one final condition, then pass both into `Trial$new()`.

``` r
set.seed(111)

arms <- c("pbo", "trt")
allocation <- c(1, 1)
sample_size <- 80

plan <- deterministic_schedule(
  sample_size = sample_size,
  arms = arms,
  allocation = allocation,
  enrollment = list(end_time = c(4, 8), rate = c(8, 2)),
  dropout = list(end_time = c(4, 8), rate = c(0, 1))
)

tmr <- Timer$new(name = "conditions_timer")
add_timepoints(tmr, plan)

n_by_arm <- c(
  pbo = sum(plan$enroll[plan$arm == "pbo"]),
  trt = sum(plan$enroll[plan$arm == "trt"])
)

pbo_data <- data.frame(
  id = seq_len(n_by_arm[["pbo"]]),
  y = rnorm(n_by_arm[["pbo"]], mean = 0),
  readout_time = 1
)

trt_data <- data.frame(
  id = seq_len(n_by_arm[["trt"]]),
  y = rnorm(n_by_arm[["trt"]], mean = 0.4),
  readout_time = 1
)

populations <- list(
  Population$new("pbo", pbo_data),
  Population$new("trt", trt_data)
)
```

``` r
interim <- condition_calendar_time(
  cal_time = 4,
  analysis = function(df, current_time) {
    dat <- df |> filter(!is.na(enroll_time))
    data.frame(
      look = "interim",
      time = current_time,
      n_total = nrow(dat),
      p_value = t.test(y ~ arm, data = dat)$p.value
    )
  },
  name = "interim"
)

final <- condition_enrollment_fraction(
  fraction = 1.0,
  sample_size = sample_size,
  analysis = function(df, current_time) {
    dat <- df |> filter(!is.na(enroll_time))
    data.frame(
      look = "final",
      time = current_time,
      n_total = nrow(dat),
      p_value = t.test(y ~ arm, data = dat)$p.value
    )
  },
  name = "final"
)
```

``` r
trial <- Trial$new(
  name = "conditions_demo",
  seed = 111,
  timer = tmr,
  population = populations,
  conditions = list(interim, final)
)

trial$run()
```

``` r
trial$results
#> $time_4
#> $time_4$interim
#>      look time n_total      p_value
#> 1 interim    4      32 0.0007027072
#> 
#> 
#> $time_9
#> $time_9$final
#>    look time n_total      p_value
#> 1 final    9      80 9.633389e-06
collect_results(trial)
#>   replicate timepoint analysis    look time n_total      p_value
#> 1         1         4  interim interim    4      32 7.027072e-04
#> 2         1         9    final   final    9      80 9.633389e-06
```

Each look writes to a different entry under `trial$results` and appears
as a separate row in
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
via the `analysis` column.

## Cooldown and max triggers

For recurrent conditions, you can cap firing frequency with:

- `cooldown`: minimum time gap between triggers
- `max_triggers`: maximum total number of triggers

These settings are useful for frequent data-driven checks where you want
guardrails against over-triggering.

## Pattern summary

1.  Define one `Condition` per planned look.
2.  Keep trigger logic explicit and readable.
3.  Use composable triggers (`&`, `|`) for richer rules.
4.  Aggregate with
    [`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md)
    for downstream operating characteristics.

## Next steps

- [Enrollment and
  Dropout](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment.md) -
  building timers for fixed or stochastic schedules
- [Trial
  reference](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md) -
  how Trial orchestrates timer, populations, and conditions
- [Example
  5](https://boehringer-ingelheim.github.io/rxsim/articles/example-5.md) -
  seamless Ph2a/2b with interim and final looks

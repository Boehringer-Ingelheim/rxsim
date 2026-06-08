# Population: Subject-Level Data and Enrollment Tracking

## What is a Population?

A `Population` is the subject-level state for one trial arm. It owns the
arm’s endpoint data and two tracking vectors: one for enrollment times
and one for dropout times. In practice, it is the object that answers
three questions at once: who belongs to this arm, what data does each
subject contribute, and what has happened to each subject so far?

The example below starts with a single-readout arm. No subject has
enrolled or dropped yet, so both state vectors start as `NA`.

``` r
set.seed(1001)

subject_data <- data.frame(
  id = 1:4,
  readout_time = 0,
  response = c(11.2, 9.8, 10.5, 12.1)
)

pop <- Population$new(name = "active", data = subject_data)

list(
  arm = pop$name,
  subjects = pop$n,
  readouts_per_subject = pop$n_readouts,
  enrolled = pop$enrolled,
  dropped = pop$dropped
)
#> $arm
#> [1] "active"
#> 
#> $subjects
#> [1] 4
#> 
#> $readouts_per_subject
#> [1] 1
#> 
#> $enrolled
#> [1] NA NA NA NA
#> 
#> $dropped
#> [1] NA NA NA NA
```

## Data frame format

`Population$data` must be a data frame with three structural columns
plus at least one endpoint column.

| Column         | Type      | Meaning                                         |
|----------------|-----------|-------------------------------------------------|
| `id`           | integer   | Unique subject identifier within the arm        |
| `arm`          | character | Arm label; auto-filled from `name` if absent    |
| `readout_time` | numeric   | Time from enrollment to the readout             |
| endpoint       | any       | At least one endpoint column such as `response` |

`readout_time` is relative to enrollment, not calendar time. A subject
can have one row or several rows, depending on how many scheduled
readouts you need to store.

Here `subject_data` omitted `arm`, so `Population$new()` filled it in
from the object name.

``` r
pop$data
#>   id readout_time response    arm
#> 1  1            0     11.2 active
#> 2  2            0      9.8 active
#> 3  3            0     10.5 active
#> 4  4            0     12.1 active
```

## Creating a Population

You can construct a `Population` directly from a data frame when you
already have subject-level rows. For quick single-endpoint examples,
[`as_population_data()`](https://boehringer-ingelheim.github.io/rxsim/reference/as_population_data.md)
is a compact shortcut that creates `id`, `data`, and `readout_time` so
that `Population$new()` can fill the arm label.

The first object reuses the explicit data frame from above. The second
shows the vector shortcut for a placebo arm.

``` r
pop_df <- Population$new(name = "active", data = subject_data)

pop_vector <- Population$new(
  name = "placebo",
  data = as_population_data(c(8.7, 9.1, 8.9, 9.4))
)

pop_df$data
#>   id readout_time response    arm
#> 1  1            0     11.2 active
#> 2  2            0      9.8 active
#> 3  3            0     10.5 active
#> 4  4            0     12.1 active
pop_vector$data
#>   id data readout_time     arm
#> 1  1  8.7            0 placebo
#> 2  2  9.1            0 placebo
#> 3  3  8.9            0 placebo
#> 4  4  9.4            0 placebo
```

## Enrollment and dropout tracking

The `enrolled` and `dropped` fields are numeric vectors with one entry
per subject. `NA` means the event has not happened yet.
`set_enrolled(n, time)` randomly chooses `n` currently unenrolled
subjects and stamps them with an enrollment time. `set_dropped(n, time)`
randomly chooses from subjects who are enrolled and not yet dropped.

A step-by-step example makes the bookkeeping visible.

``` r
set.seed(2024)

state_0 <- data.frame(
  id = seq_len(pop$n),
  enrolled = pop$enrolled,
  dropped = pop$dropped
)

pop$set_enrolled(n = 2, time = 1)
state_1 <- data.frame(
  id = seq_len(pop$n),
  enrolled = pop$enrolled,
  dropped = pop$dropped
)

pop$set_enrolled(n = 1, time = 3)
state_2 <- data.frame(
  id = seq_len(pop$n),
  enrolled = pop$enrolled,
  dropped = pop$dropped
)

pop$set_dropped(n = 1, time = 5)
state_3 <- data.frame(
  id = seq_len(pop$n),
  enrolled = pop$enrolled,
  dropped = pop$dropped
)

state_0
#>   id enrolled dropped
#> 1  1       NA      NA
#> 2  2       NA      NA
#> 3  3       NA      NA
#> 4  4       NA      NA
state_1
#>   id enrolled dropped
#> 1  1        1      NA
#> 2  2        1      NA
#> 3  3       NA      NA
#> 4  4       NA      NA
state_2
#>   id enrolled dropped
#> 1  1        1      NA
#> 2  2        1      NA
#> 3  3        3      NA
#> 4  4       NA      NA
state_3
#>   id enrolled dropped
#> 1  1        1      NA
#> 2  2        1      NA
#> 3  3        3       5
#> 4  4       NA      NA
```

## Multiple readout times

When a subject contributes more than one scheduled measurement, store
multiple rows for the same `id`. `Population$n` still counts unique
subjects, while `Population$n_readouts` records how many rows belong to
each subject.

The data below uses two readouts per subject: one at enrollment and one
eight units later.

``` r
multi_readout_data <- data.frame(
  id = rep(1:4, each = 2),
  readout_time = rep(c(0, 8), times = 4),
  response = c(11.2, 10.8, 9.8, 9.4, 10.5, 10.0, 12.1, 11.6)
)

multi_pop <- Population$new(name = "active", data = multi_readout_data)

multi_pop$n
#> [1] 4
multi_pop$n_readouts
#> [1] 2
multi_pop$data
#>   id readout_time response    arm
#> 1  1            0     11.2 active
#> 2  1            8     10.8 active
#> 3  2            0      9.8 active
#> 4  2            8      9.4 active
#> 5  3            0     10.5 active
#> 6  3            8     10.0 active
#> 7  4            0     12.1 active
#> 8  4            8     11.6 active
```

## Endpoint types

`Population` itself is endpoint-agnostic: it stores rows and metadata,
not a special endpoint class. Endpoint type is determined by which
endpoint columns you include in the data frame.

| Endpoint type | Typical data construction                                                            | Example column |
|---------------|--------------------------------------------------------------------------------------|----------------|
| Continuous    | [`rnorm()`](https://rdrr.io/r/stats/Normal.html) draws in a numeric column           | `response`     |
| Binary        | `rbinom(..., size = 1, prob = p)` as 0/1                                             | `response_bin` |
| Time-to-event | [`rexp()`](https://rdrr.io/r/stats/Exponential.html) / Weibull draws for event times | `tte`          |

This matters because the analysis function in each `Condition` must
match the endpoint scale: t-test or linear model for continuous
endpoints, proportion or logistic methods for binary endpoints, and
survival methods for time-to-event outcomes.

``` r
set.seed(4040)

pop_cont <- Population$new(
  name = "cont",
  data = data.frame(
    id = 1:5,
    response = rnorm(5, mean = 0, sd = 1),
    readout_time = 1
  )
)

pop_bin <- Population$new(
  name = "bin",
  data = data.frame(
    id = 1:5,
    response_bin = rbinom(5, size = 1, prob = 0.35),
    readout_time = 1
  )
)

pop_tte <- Population$new(
  name = "tte",
  data = data.frame(
    id = 1:5,
    tte = rexp(5, rate = 0.2),
    readout_time = 1
  )
)

names(pop_cont$data)
#> [1] "id"           "response"     "readout_time" "arm"
names(pop_bin$data)
#> [1] "id"           "response_bin" "readout_time" "arm"
names(pop_tte$data)
#> [1] "id"           "tte"          "readout_time" "arm"
```

## Role in Trial

`Population` is one arm-level building block inside `Trial`. At each
timer update, `Trial$run()` calls `set_enrolled()` and `set_dropped()`
for each arm, then builds a locked snapshot from the currently enrolled
subjects. That snapshot carries forward the original population columns
and adds `enroll_time`, `drop_time`, `measurement_time`, and the current
trial `time`.

This example defines one simple analysis so that the trial stores a
`locked_data` snapshot at calendar time 5.

``` r
set.seed(303)

trial_pop <- Population$new(name = "active", data = subject_data)

timer <- Timer$new(name = "population_timer")
timer$add_timepoint(time = 1, arm = "active", enroll = 2L, drop = 0L)
timer$add_timepoint(time = 3, arm = "active", enroll = 1L, drop = 0L)
timer$add_timepoint(time = 5, arm = "active", enroll = 0L, drop = 1L)

snapshot_condition <- Condition$new(
  where = calendar_trigger(5),
  analysis = function(df, current_time) {
    data.frame(
      current_time = current_time,
      n_subjects = length(unique(df$id))
    )
  },
  name = "snapshot"
)

trial <- Trial$new(
  name = "population_trial",
  seed = 303,
  timer = timer,
  population = list(trial_pop),
  conditions = list(snapshot_condition)
)

trial$run()

trial$locked_data$time_5[
  , c(
    "id",
    "arm",
    "readout_time",
    "response",
    "enroll_time",
    "drop_time",
    "measurement_time",
    "time"
  )
]
#>   id    arm readout_time response enroll_time drop_time measurement_time time
#> 1  1 active            0     11.2           3        NA                3    5
#> 3  3 active            0     10.5           1         5                1    5
#> 4  4 active            0     12.1           1        NA                1    5

trial$results$time_5$snapshot
#>   current_time n_subjects
#> 1            5          3
```

## Tips and gotchas

A few rules matter when you build populations by hand. `id` should
identify the subject, not the row, so repeated readouts reuse the same
`id`. Omitting `arm` is fine because `Population$new()` fills it from
`name`. Finally, the data frame must have at least four columns after
`arm` is added: the three structural columns and at least one endpoint
column.

The checks below summarise those rules on the objects already created
above.

``` r
data.frame(
  check = c(
    "Unique subjects in multi_pop",
    "Rows per subject in multi_pop",
    "Arm auto-filled in pop",
    "Minimum column count met"
  ),
  value = c(
    multi_pop$n,
    paste(unique(as.integer(table(multi_pop$data$id))), collapse = ", "),
    all(pop$data$arm == pop$name),
    ncol(pop$data) >= 4
  )
)
#>                           check value
#> 1  Unique subjects in multi_pop     4
#> 2 Rows per subject in multi_pop     2
#> 3        Arm auto-filled in pop  TRUE
#> 4      Minimum column count met  TRUE
```

## Next steps

To see how `Population` fits into the broader workflow, continue with:

- [Enrollment and
  Dropout](https://boehringer-ingelheim.github.io/rxsim/articles/enrollment.md):
  stochastic and deterministic enrollment, Timer reference
- [Trial
  reference](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md):
  orchestration across population, timer, and conditions

# Condition: Stateful trigger and analysis unit

A `Condition` encapsulates a single trigger rule that is evaluated
against a data snapshot at each simulated timepoint. It combines three
concerns:

1.  **Filtering** - a
    [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
    expression selects the rows relevant to this condition (e.g. "only
    enrolled subjects in arm A").

2.  **Analysis** - an optional function transforms the filtered snapshot
    into a result (e.g. a t-test, a subject count, a Go/No-Go decision).

3.  **Trigger bookkeeping** - the condition fires only when the filtered
    data is non-empty, the cooldown period has elapsed since the last
    trigger, and the maximum trigger count has not been reached.

`Condition` objects are stored in `trial$conditions` and evaluated by
[`Trial`](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)`$run()`
at each timepoint.

## Details

**Three-gate logic.** A trigger fires only when all three gates pass:

1.  The filtered snapshot contains at least one row.

2.  `trigger_count < max_triggers`.

3.  `current_time - last_trigger_time >= cooldown` (or the condition has
    never fired before).

If any gate fails, `check_conditions()` returns an empty list and state
is not updated.

On a successful trigger, the condition calls
`analysis(df, current_time, ...)` and stores the result under `name` (or
`1L` when no name is set). Any values in `analysis_args` are appended as
additional named arguments. If no analysis function is provided, the
filtered data frame is returned as-is with a warning.

## See also

- [`Timer`](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md)
  for managing trial timepoints

- [`Trial`](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)
  for running the simulation and iterating over conditions

- [`value_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md),
  [`count_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md),
  [`enroll_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md),
  [`calendar_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)
  for building trigger specifications

- [`condition_calendar_time()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_calendar_time.md)
  and
  [`condition_enrollment_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_enrollment_fraction.md)
  for convenient `Condition` constructors

## Public fields

- `where`:

  `list` of quosures
  ([`rlang::quos()`](https://rlang.r-lib.org/reference/defusing-advanced.html))
  used as
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  predicates, or a `trigger` object (converted automatically). `NULL` or
  empty list passes the full snapshot.

- `trigger_spec`:

  The original `trigger` object passed to `where` (before quosure
  conversion), or `NULL`. Used by the fixed fast path to compute the
  earliest possible firing time and skip non-firing timepoints. `NULL`
  means "evaluate at every timepoint" (safe fallback).

- `analysis`:

  `function` or `NULL`. Called as `analysis(df, current_time, ...)` on a
  successful trigger, where `...` are any values from `analysis_args`.

- `analysis_args`:

  `list` or `NULL`. Named list of extra values injected into the
  analysis function call as additional named arguments.

- `name`:

  `character` or `NULL`. Key labelling the result in the output list.
  Falls back to `1L` when `NULL`.

- `cooldown`:

  `numeric`. Minimum time units between consecutive triggers. Default
  `0`.

- `max_triggers`:

  `integer` or `Inf`. Maximum number of times this condition may fire.
  Default `1L`.

- `trigger_count`:

  `integer`. Number of successful triggers so far. Initialised to `0L`.

- `last_trigger_time`:

  `numeric`. Calendar time of the most recent successful trigger.
  `NA_real_` until first trigger.

## Methods

### Public methods

- [`Condition$new()`](#method-Condition-initialize)

- [`Condition$check_conditions()`](#method-Condition-check_conditions)

- [`Condition$clone()`](#method-Condition-clone)

------------------------------------------------------------------------

### `Condition$new()`

Create a new `Condition` instance.

#### Usage

    Condition$new(
      where = NULL,
      analysis = NULL,
      analysis_args = NULL,
      name = NULL,
      cooldown = 0,
      max_triggers = 1L
    )

#### Arguments

- `where`:

  `trigger` object (from
  [`value_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md),
  [`count_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md),
  [`enroll_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md),
  or
  [`calendar_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md)),
  a `list` of quosures from
  [`rlang::quos()`](https://rlang.r-lib.org/reference/defusing-advanced.html),
  or `NULL` to use the full snapshot.

- `analysis`:

  `function` or `NULL`. Called as `analysis(df, current_time, ...)` on a
  successful trigger, where `...` are the values from `analysis_args`.

- `analysis_args`:

  `list` or `NULL`. Named list of extra arguments passed to the analysis
  function after `df` and `current_time`.

- `name`:

  `character` or `NULL`. Result key. Defaults to `1L`.

- `cooldown`:

  `numeric`. Minimum time between triggers. Default `0`.

- `max_triggers`:

  `integer`. Maximum trigger count. Default `1L`. Use `Inf` for
  unlimited.

#### Returns

A new `Condition` instance.

------------------------------------------------------------------------

### `Condition$check_conditions()`

Evaluate this condition against a data snapshot.

Applies the three-gate logic: non-empty filter result, cooldown elapsed,
and trigger count below `max_triggers`. Returns the analysis result (or
filtered data) on a successful trigger, or an empty list otherwise.

#### Usage

    Condition$check_conditions(locked_data, current_time)

#### Arguments

- `locked_data`:

  `data.frame` The trial snapshot at the current time.

- `current_time`:

  `numeric` Calendar time of the current timepoint.

#### Returns

Named `list` with one entry (the analysis result) on success, or an
empty `list` if the condition did not fire.

------------------------------------------------------------------------

### `Condition$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Condition$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# Build a snapshot data frame (enroll_time = NA means not yet enrolled)
snapshot <- data.frame(
  arm         = c("A", "A", "A", "B"),
  status      = c("active", "active", "active", "active"),
  enroll_time = c(1, 2, 3, NA_real_),
  stringsAsFactors = FALSE
)

# Analysis function: count enrolled subjects and record fire time
count_fn <- function(df, current_time) {
  data.frame(n_active = nrow(df), fired_at = current_time)
}

# Condition fires once when 3+ of 4 subjects are enrolled (max_triggers = 1)
cond <- Condition$new(
  where        = enroll_trigger(fraction = 0.75, sample_size = 4),
  analysis     = count_fn,
  name         = "interim_A",
  cooldown     = 0,
  max_triggers = 1L
)

# First call: fires and returns analysis result
res <- cond$check_conditions(snapshot, current_time = 5)
res[["interim_A"]]  # data.frame(n_active = 4, fired_at = 5)
#>   n_active fired_at
#> 1        3        5

# Second call: does not fire (max_triggers already reached)
res2 <- cond$check_conditions(snapshot, current_time = 6)
length(res2)  # 0
#> [1] 0
```

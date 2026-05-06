# Condition: Stateful trigger and analysis unit

A `Condition` encapsulates a single trigger rule that is evaluated
against a data snapshot at each simulated timepoint. It combines three
concerns:

1.  **Filtering** — a
    [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
    expression selects the rows relevant to this condition (e.g. "only
    enrolled subjects in arm A").

2.  **Analysis** — an optional function transforms the filtered snapshot
    into a result (e.g. a t-test, a subject count, a Go/No-Go decision).

3.  **Trigger bookkeeping** — the condition fires only when the filtered
    data is non-empty, the cooldown period has elapsed since the last
    trigger, and the maximum trigger count has not been reached.

`Condition` objects are stored in `trial$conditions` and evaluated by
[`Trial`](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)`$run()`
at each timepoint.

## Details

**Three-gate logic.** A trigger fires only when all three gates pass:

1.  The filtered snapshot contains at least one row.

2.  `current_time - last_trigger_time >= cooldown` (or the condition has
    never fired before).

3.  `trigger_count < max_triggers`.

If any gate fails, `check_conditions()` returns an empty list and state
is not updated.

On a successful trigger, the condition calls
`analysis(filtered_data, current_time)` and stores the result under
`name` (or `1L` when no name is set). If no analysis function is
provided, the filtered data frame is returned as-is with a warning.

## Fields

- `where`:

  `list` of quosures (from
  [`rlang::quos()`](https://rlang.r-lib.org/reference/defusing-advanced.html))
  used as
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  predicates. Pass `NULL` or an empty list to skip filtering and pass
  the full snapshot to the analysis.

- `analysis`:

  `function` or `NULL`. Called as
  `analysis(filtered_data, current_time)` on a successful trigger.
  Should return a `data.frame` or named list. If `NULL`, the filtered
  data frame is returned with a warning.

- `name`:

  `character` or `NULL`. Key used to label the result in the returned
  list. Falls back to `1L` when `NULL`.

- `cooldown`:

  `numeric`. Minimum time units that must elapse between consecutive
  triggers. Default `0` (no cooldown).

- `max_triggers`:

  `integer`. Maximum number of times this condition may fire. Use `Inf`
  for unlimited. Default `1L`.

- `trigger_count`:

  `integer`. Number of successful triggers so far. Initialised to `0L`.

- `last_trigger_time`:

  `numeric`. Calendar time of the most recent successful trigger.
  Initialised to `NA_real_`.

## Methods

- `$new(where, analysis, name, cooldown, max_triggers)`:

  Construct a new `Condition`. All arguments except `where` are
  optional. `cooldown` must be a single non-negative number;
  `max_triggers` must be a single non-negative integer or `Inf`.

- `$check_conditions(locked_data, current_time)`:

  Evaluate the condition against `locked_data` at `current_time`.
  Returns a named `list` containing the analysis result (or filtered
  data frame) if the condition fires, or an empty `list` otherwise. On a
  successful trigger, `trigger_count` is incremented and
  `last_trigger_time` is updated.

## See also

- [`Timer`](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md)
  for managing trial timepoints

- [`Trial`](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)
  for running the simulation and iterating over conditions

- [`trigger_by_calendar()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_calendar.md)
  and
  [`trigger_by_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_fraction.md)
  for convenient `Condition` constructors

- [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  for predicate syntax

## Public fields

- `where`:

  `list` of quosures
  ([`rlang::quos()`](https://rlang.r-lib.org/reference/defusing-advanced.html))
  used as
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  predicates. `NULL` or empty list passes the full snapshot.

- `analysis`:

  `function` or `NULL`. Called as
  `analysis(filtered_data, current_time)` on a successful trigger.

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

- [`Condition$new()`](#method-Condition-new)

- [`Condition$check_conditions()`](#method-Condition-check_conditions)

- [`Condition$clone()`](#method-Condition-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `Condition` instance.

#### Usage

    Condition$new(
      where = NULL,
      analysis = NULL,
      name = NULL,
      cooldown = 0,
      max_triggers = 1L
    )

#### Arguments

- `where`:

  `list` of quosures (from
  [`rlang::quos()`](https://rlang.r-lib.org/reference/defusing-advanced.html))
  used as filter predicates. Pass `NULL` or omit to use the full
  snapshot.

- `analysis`:

  `function` or `NULL`. Called as
  `analysis(filtered_data, current_time)` on a successful trigger.

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

### Method `check_conditions()`

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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Condition$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# Build a snapshot data frame
snapshot <- data.frame(
  arm    = c("A", "A", "A", "B"),
  status = c("active", "active", "active", "active"),
  stringsAsFactors = FALSE
)

# Analysis function: count active subjects per arm
count_fn <- function(dat, current_time) {
  data.frame(n_active = nrow(dat), fired_at = current_time)
}

# Condition fires once when arm A has active subjects (max_triggers = 1)
cond <- Condition$new(
  where        = rlang::quos(arm == "A", status == "active"),
  analysis     = count_fn,
  name         = "interim_A",
  cooldown     = 0,
  max_triggers = 1L
)

# First call: fires and returns analysis result
res <- cond$check_conditions(snapshot, current_time = 5)
res[["interim_A"]]  # data.frame(n_active = 3, fired_at = 5)
#>   n_active fired_at
#> 1        3        5

# Second call: does not fire (max_triggers already reached)
res2 <- cond$check_conditions(snapshot, current_time = 6)
length(res2)  # 0
#> [1] 0
```

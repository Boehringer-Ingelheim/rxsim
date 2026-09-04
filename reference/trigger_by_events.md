# Trigger Analysis When Enough TTE Events Have Occurred

Builds a
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
that fires once a target number of time-to-event (TTE) endpoint events
have been observed by the current trial time. Enrollment guards are
always included: only subjects who are enrolled before the current time
are eligible.

## Usage

``` r
trigger_by_events(
  event_col,
  n_events,
  time_col = "time",
  enroll_col = "enroll_time",
  analysis = NULL,
  name = NULL,
  op = ">="
)
```

## Arguments

- event_col:

  `character` Column holding each subject's event calendar time. `NA`
  indicates no event has occurred yet.

- n_events:

  `numeric` Number of events required to fire.

- time_col:

  `character` Column holding the current trial time. Defaults to
  `"time"` (the column added by `Trial$run()`).

- enroll_col:

  `character` Column holding each subject's enrollment time. Defaults to
  `"enroll_time"`.

- analysis:

  `function` or `NULL` Optional analysis function called as
  `analysis(df, current_time, ...)`. If `NULL`, the filtered snapshot is
  returned as-is with a warning.

- name:

  `character` or `NULL` Result key. Defaults to `"events_<n_events>"`.

- op:

  `character` Comparison operator for the event-count threshold. Must be
  one of `c(">=", "<=", ">", "<", "==", "!=", "%in%")`. Defaults to
  `">="`.

## Value

A
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
object.

## See also

[Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md),
[`condition_calendar_time()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_calendar_time.md),
[`condition_enrollment_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_enrollment_fraction.md),
[`timed_count_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md),
[`notna_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md),
[`col_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md).

## Examples

``` r
cond <- trigger_by_events(
  event_col = "pfs_event_time",
  n_events  = 100,
  analysis  = function(df, current_time) {
    data.frame(n_events = sum(!is.na(df$pfs_event_time)), fired_at = current_time)
  }
)
```

# Build a Condition that Fires at a Calendar Time

Convenience constructor: wraps `Condition$new()` with a
[`calendar_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md).
The returned `Condition` should be passed to
`Trial$new(conditions = list(...))`.

## Usage

``` r
condition_calendar_time(cal_time, analysis = NULL, name = NULL)
```

## Arguments

- cal_time:

  `numeric` Calendar time at or after which to fire.

- analysis:

  `function` or `NULL` Called as `analysis(df, current_time, ...)`. If
  `NULL`, the filtered snapshot is returned as-is with a warning.

- name:

  `character` or `NULL` Result key. Defaults to `"cal_time_<cal_time>"`.

## Value

A
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
object.

## See also

[Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md),
[`condition_enrollment_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_enrollment_fraction.md),
[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md),
[`calendar_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md).

## Examples

``` r
cond <- condition_calendar_time(
  cal_time = 12,
  analysis = function(df, current_time) {
    data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
  }
)
```

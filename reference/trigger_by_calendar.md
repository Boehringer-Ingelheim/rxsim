# Trigger Analysis at a Calendar Time

Builds a
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
that fires when the trial clock reaches a specified calendar time. The
returned `Condition` should be passed to
`Trial$new(conditions = list(...))`.

## Usage

``` r
trigger_by_calendar(cal_time, analysis = NULL, name = NULL)
```

## Arguments

- cal_time:

  `numeric` Calendar time(s) at which to trigger.

- analysis:

  `function` or `NULL` Optional analysis function called as
  `analysis(df, current_time, ...)`. If `NULL`, the filtered snapshot is
  returned as-is with a warning.

- name:

  `character` or `NULL` Result key. Defaults to `"cal_time_<cal_time>"`.

## Value

A
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
object.

## See also

[Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md),
[`trigger_by_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_fraction.md),
[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md),
[`calendar_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md).

## Examples

``` r
cond <- trigger_by_calendar(
  cal_time = 12,
  analysis = function(df, current_time) {
    data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
  }
)
```

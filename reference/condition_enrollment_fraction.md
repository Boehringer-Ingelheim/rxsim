# Build a Condition that Fires at an Enrollment Fraction

Convenience constructor: wraps `Condition$new()` with an
[`enroll_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md).
The returned `Condition` should be passed to
`Trial$new(conditions = list(...))`.

## Usage

``` r
condition_enrollment_fraction(
  fraction,
  sample_size,
  analysis = NULL,
  name = NULL
)
```

## Arguments

- fraction:

  `numeric` Sample fraction (0 \< fraction \<= 1).

- sample_size:

  `integer` Target sample size.

- analysis:

  `function` or `NULL` Called as `analysis(df, current_time, ...)`. If
  `NULL`, the filtered snapshot is returned as-is with a warning.

- name:

  `character` or `NULL` Result key. Defaults to `"frac_<fraction>"`.

## Value

A
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
object.

## See also

[Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md),
[`condition_calendar_time()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_calendar_time.md),
[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md),
[`enroll_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md).

## Examples

``` r
cond <- condition_enrollment_fraction(
  fraction    = 0.5,
  sample_size = 100,
  analysis    = function(df, current_time) {
    data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
  }
)
```

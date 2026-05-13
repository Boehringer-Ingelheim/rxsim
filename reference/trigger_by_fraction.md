# Trigger Analysis at a Sample Fraction

Builds a
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
that fires when a given fraction of the target sample has been enrolled.
The returned `Condition` should be passed to
`Trial$new(conditions = list(...))`.

## Usage

``` r
trigger_by_fraction(fraction, sample_size, analysis = NULL, name = NULL)
```

## Arguments

- fraction:

  `numeric` Sample fraction (0 \< fraction \<= 1).

- sample_size:

  `integer` Target sample size.

- analysis:

  `function` or `NULL` Optional analysis function called as
  `analysis(df, current_time, ...)`. If `NULL`, the filtered snapshot is
  returned as-is with a warning.

- name:

  `character` or `NULL` Result key. Defaults to `"frac_<fraction>"`.

## Value

A
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
object.

## See also

[Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md),
[`trigger_by_calendar()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_calendar.md),
[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md),
[`enroll_trigger()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_primitives.md).

## Examples

``` r
cond <- trigger_by_fraction(
  fraction    = 0.5,
  sample_size = 100,
  analysis    = function(df, current_time) {
    data.frame(n_enrolled = sum(!is.na(df$enroll_time)))
  }
)
```

# Trigger Analysis at a Sample Fraction

Adds an analysis trigger when a fraction of the target sample enrolls.

## Usage

``` r
trigger_by_fraction(fraction, timer, sample_size, analysis = NULL)
```

## Arguments

- fraction:

  `numeric` Sample fraction (0 \< fraction \<= 1).

- timer:

  [`Timer`](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md)
  instance to update.

- sample_size:

  `integer` Target sample size.

- analysis:

  `function` or `NULL` Optional function to apply.

## Value

Invisible
[`Timer`](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md).

## See also

[Timer](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md),
[`trigger_by_calendar()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_calendar.md).

## Examples

``` r
t <- Timer$new("timer")
trigger_by_fraction(0.5, t, sample_size = 100, analysis = function(df, current_time) nrow(df))
```

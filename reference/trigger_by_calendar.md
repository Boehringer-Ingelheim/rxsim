# Trigger Analysis at a Calendar Time

Adds an analysis trigger at a specified calendar time.

## Usage

``` r
trigger_by_calendar(cal_time, timer, analysis = NULL)
```

## Arguments

- cal_time:

  `numeric` Calendar time(s) to trigger.

- timer:

  [`Timer`](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md)
  instance to update.

- analysis:

  `function` or `NULL` Optional function to apply.

## Value

Invisible
[`Timer`](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md).

## See also

[Timer](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md),
[`trigger_by_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_fraction.md).

## Examples

``` r
t <- Timer$new("timer")
trigger_by_calendar(2, t, analysis = function(df, current_time) nrow(df))
```

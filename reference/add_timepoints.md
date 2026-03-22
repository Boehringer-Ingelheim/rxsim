# Add Timepoints to a Timer

Adds multiple enrollment and dropout events from a data frame.

## Usage

``` r
add_timepoints(timer, df)
```

## Arguments

- timer:

  [`Timer`](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md)
  instance.

- df:

  `data.frame` with columns: `time` (numeric), `arm` (character),
  `enroller` (integer), `dropper` (integer).

## See also

[Timer](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md),
[`gen_plan()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_plan.md),
[`gen_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/gen_timepoints.md).

## Examples

``` r
t <- Timer$new(name = "Timer")

timepoints <- data.frame(
  time = c(1, 2, 3.1, 4, 5, 6),
  arm = rep("Arm A", 6),
  dropper = c(2L, rep(1L, 5)),
  enroller = rep(3L, 6)
)

add_timepoints(t, timepoints)
```

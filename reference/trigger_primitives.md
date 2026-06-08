# Build Trial Triggers

Create trigger specifications that can be passed to
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
or composed with `&` and `|`.

## Usage

``` r
value_trigger(col, op, rhs)

count_trigger(col, op, rhs)

enroll_trigger(fraction, sample_size)

calendar_trigger(cal_time)

# S3 method for class 'trigger'
e1 & e2

# S3 method for class 'trigger'
e1 | e2
```

## Arguments

- col:

  `character` Column name referenced by the trigger.

- op:

  `character` Comparison operator. Must be one of
  `c(">=", "<=", ">", "<", "==", "!=", "%in%")`.

- rhs:

  Right-hand side value. Must be atomic for `value_trigger()` and
  numeric for `count_trigger()`.

- fraction:

  `numeric` Sample fraction (0 \< fraction \<= 1).

- sample_size:

  `numeric` Target sample size.

- cal_time:

  `numeric` Calendar time at or after which to trigger. Uses `>=` so it
  works with both stochastic (continuous) and deterministic (integer)
  schedules.

- e1, e2:

  `trigger` objects to combine.

## Value

A `trigger` object.

## See also

[Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md),
`enroll_trigger()`, `calendar_trigger()`.

## Examples

``` r
t1 <- value_trigger("time", ">=", 52)
t2 <- count_trigger("enroll_time", ">=", 100)
t3 <- enroll_trigger(0.5, 200) & calendar_trigger(52)
```

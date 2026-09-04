# Timer: Track timed events across arms

A class to collect and query *timepoints* - time-based enrollment and
dropout events - across trial arms.

Use `add_schedule()` to register events and `get_end_timepoint()` /
`get_n_arms()` / `get_unique_times()` for summary queries. The full
event table is the public `timelist` field.

## Details

Trigger conditions (filtering + analysis) are now managed by the
separate
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
class. `Condition` objects are stored in `trial$conditions` and
evaluated by
[`Trial`](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)`$run()`
at each timepoint.

Helper functions
[`condition_calendar_time()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_calendar_time.md)
and
[`condition_enrollment_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_enrollment_fraction.md)
provide convenient shortcuts for building `Condition` objects; both
return a
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
that you pass to `Trial$new(conditions = list(...))`.

## See also

[`Trial`](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)
to coordinate simulations with populations,
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
for trigger/analysis logic,
[`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)
/
[`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md)
to build a schedule.

## Public fields

- `name`:

  `character` Unique identifier for the `Timer` instance.

- `timelist`:

  `data.frame` A data.frame of timepoints with columns:

  - `time` `numeric` Calendar time

  - `arm` `character` Unique identifier of the arm

  - `drop` `integer` \# of subjects dropped at `time`

  - `enroll` `integer` \# of subjects enrolled at `time`

## Methods

### Public methods

- [`Timer$new()`](#method-Timer-initialize)

- [`Timer$add_schedule()`](#method-Timer-add_schedule)

- [`Timer$get_end_timepoint()`](#method-Timer-get_end_timepoint)

- [`Timer$get_n_arms()`](#method-Timer-get_n_arms)

- [`Timer$get_unique_times()`](#method-Timer-get_unique_times)

- [`Timer$clone()`](#method-Timer-clone)

------------------------------------------------------------------------

### `Timer$new()`

Create a new `Timer` instance.

#### Usage

    Timer$new(name, timelist = NULL)

#### Arguments

- `name`:

  `character` Unique identifier.

- `timelist`:

  `data.frame` Optional data.frame of timepoints with columns `time`,
  `arm`, `drop`, `enroll`. If `NULL`, an empty frame is created.

#### Returns

A new `Timer` instance.

#### Examples

    t <- Timer$new(name = "Timer")

------------------------------------------------------------------------

### `Timer$add_schedule()`

Add a schedule of timepoints to the timer.

#### Usage

    Timer$add_schedule(schedule)

#### Arguments

- `schedule`:

  `data.frame` with columns `time` (numeric), `arm` (character),
  `enroll` (integer), `drop` (integer). One row per event; a single
  event is a one-row data frame. Typically the output of
  [`stochastic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/stochastic_schedule.md)
  or
  [`deterministic_schedule()`](https://boehringer-ingelheim.github.io/rxsim/reference/deterministic_schedule.md).

  `enroll` and `drop` are subject counts and must be integer (`3L`, not
  `3`) - fractional counts are silently truncated downstream, so they
  are rejected here.

#### Examples

    t <- Timer$new(name = "Timer")

    # single event
    t$add_schedule(data.frame(time = 1, arm = "A", drop = 1L, enroll = 3L))

    # whole schedule (data.frame() recycles the constant columns)
    t$add_schedule(data.frame(time = 2:3, arm = "A", enroll = 2L, drop = 0L))

------------------------------------------------------------------------

### `Timer$get_end_timepoint()`

Determine the last timepoint for a given instance of `Timer` class.

#### Usage

    Timer$get_end_timepoint()

#### Examples

    t <- Timer$new(name = "Timer")
    t$add_schedule(data.frame(time = 3.14, arm = "A", drop = 7L, enroll = 22L))
    t$get_end_timepoint()

------------------------------------------------------------------------

### `Timer$get_n_arms()`

Get number of unique arms.

#### Usage

    Timer$get_n_arms()

#### Returns

`integer` Number of unique arms.

#### Examples

    t <- Timer$new(name = "Timer")
    t$add_schedule(data.frame(time = 3.14, arm = "A", drop = 7L, enroll = 22L))
    t$add_schedule(data.frame(time = 3.28, arm = "B", drop = 6L, enroll = 23L))
    t$get_n_arms()

------------------------------------------------------------------------

### `Timer$get_unique_times()`

Get unique timepoints.

#### Usage

    Timer$get_unique_times()

#### Returns

`numeric` vector of unique times.

#### Examples

    t <- Timer$new(name = "Timer")
    t$add_schedule(data.frame(time = 3.14, arm = "A", drop = 7L, enroll = 22L))
    t$add_schedule(data.frame(time = 3.28, arm = "B", drop = 6L, enroll = 23L))
    t$get_unique_times()

------------------------------------------------------------------------

### `Timer$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Timer$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# Basic construction
t <- Timer$new(name = "Timer")

# Add timepoints
t$add_schedule(data.frame(
  time   = c(1, 2, 1),
  arm    = c("A", "A", "B"),
  drop   = c(2L, 1L, 0L),
  enroll = c(10L, 12L, 8L)
))

# Query
t$get_end_timepoint() # max time => 2
#> [1] 2
t$get_n_arms()        # unique arms => 2
#> [1] 2
t$get_unique_times()  # unique times => c(1, 2)
#> [1] 1 2
t$timelist            # the full event table
#>   time arm drop enroll
#> 1    1   A    2     10
#> 2    2   A    1     12
#> 3    1   B    0      8


## ------------------------------------------------
## Method `Timer$new()`
## ------------------------------------------------

t <- Timer$new(name = "Timer")

## ------------------------------------------------
## Method `Timer$add_schedule()`
## ------------------------------------------------

t <- Timer$new(name = "Timer")

# single event
t$add_schedule(data.frame(time = 1, arm = "A", drop = 1L, enroll = 3L))

# whole schedule (data.frame() recycles the constant columns)
t$add_schedule(data.frame(time = 2:3, arm = "A", enroll = 2L, drop = 0L))

## ------------------------------------------------
## Method `Timer$get_end_timepoint()`
## ------------------------------------------------

t <- Timer$new(name = "Timer")
t$add_schedule(data.frame(time = 3.14, arm = "A", drop = 7L, enroll = 22L))
t$get_end_timepoint()
#> [1] 3.14

## ------------------------------------------------
## Method `Timer$get_n_arms()`
## ------------------------------------------------

t <- Timer$new(name = "Timer")
t$add_schedule(data.frame(time = 3.14, arm = "A", drop = 7L, enroll = 22L))
t$add_schedule(data.frame(time = 3.28, arm = "B", drop = 6L, enroll = 23L))
t$get_n_arms()
#> [1] 2

## ------------------------------------------------
## Method `Timer$get_unique_times()`
## ------------------------------------------------

t <- Timer$new(name = "Timer")
t$add_schedule(data.frame(time = 3.14, arm = "A", drop = 7L, enroll = 22L))
t$add_schedule(data.frame(time = 3.28, arm = "B", drop = 6L, enroll = 23L))
t$get_unique_times()
#> [1] 3.14 3.28
```

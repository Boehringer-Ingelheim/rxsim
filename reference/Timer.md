# Timer: Track timed events across arms

A class to collect and query *timepoints* — time-based enrollment and
dropout events — across trial arms.

Use `add_timepoint()` to register events, `get_timepoint()` for lookup,
`get_end_timepoint()` / `get_n_arms()` / `get_unique_times()` for
summary queries.

## Details

Trigger conditions (filtering + analysis) are now managed by the
separate
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
class. `Condition` objects are stored in `trial$conditions` and
evaluated by
[`Trial`](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)`$run()`
at each timepoint.

Helper functions
[`trigger_by_calendar()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_calendar.md)
and
[`trigger_by_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_fraction.md)
provide convenient shortcuts for building `Condition` objects; both
return a
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
that you pass to `Trial$new(conditions = list(...))`.

## See also

[`Trial`](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)
to coordinate simulations with populations,
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
for trigger/analysis logic,
[`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md)
to attach multiple timepoints.

## Public fields

- `name`:

  `character` Unique identifier for the `Timer` instance.

- `timelist`:

  `list` A list of timepoints. Each timepoint is a list with keys:

  - `time` `numeric` Calendar time

  - `arm` `character` Unique identifier of the arm

  - `dropper` `integer` \# of subjects dropped at `time`

  - `enroller` `integer` \# of subjects enrolled at `time`

## Methods

### Public methods

- [`Timer$new()`](#method-Timer-new)

- [`Timer$add_timepoint()`](#method-Timer-add_timepoint)

- [`Timer$get_end_timepoint()`](#method-Timer-get_end_timepoint)

- [`Timer$get_n_arms()`](#method-Timer-get_n_arms)

- [`Timer$get_unique_times()`](#method-Timer-get_unique_times)

- [`Timer$get_timepoint()`](#method-Timer-get_timepoint)

- [`Timer$clone()`](#method-Timer-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `Timer` instance.

#### Usage

    Timer$new(name, timelist = NULL)

#### Arguments

- `name`:

  `character` Unique identifier.

- `timelist`:

  `list` Optional list of timepoints.

#### Returns

A new `Timer` instance.

#### Examples

    t <- Timer$new(name = "Timer")

------------------------------------------------------------------------

### Method `add_timepoint()`

Add a timepoint to a timer.

#### Usage

    Timer$add_timepoint(time, arm, dropper, enroller)

#### Arguments

- `time`:

  `numeric` Calendar time.

- `arm`:

  `character` Arm identifier.

- `dropper`:

  `integer` Count of subjects to drop.

- `enroller`:

  `integer` Count of subjects to enroll.

#### Examples

    t <- Timer$new(name = "Timer")
    t$add_timepoint(
      time = 1,
      arm = "A",
      dropper = 1L,
      enroller = 3L
    )

------------------------------------------------------------------------

### Method `get_end_timepoint()`

Determine the last timepoint for a given instance of `Timer` class.

#### Usage

    Timer$get_end_timepoint()

#### Examples

    t <- Timer$new(name = "Timer")
    t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
    t$get_end_timepoint()

------------------------------------------------------------------------

### Method `get_n_arms()`

Get number of unique arms.

#### Usage

    Timer$get_n_arms()

#### Returns

`integer` Number of unique arms.

#### Examples

    t <- Timer$new(name = "Timer")
    t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
    t$add_timepoint(time = 3.28, arm = "B", dropper = 6L, enroller = 23L)
    t$get_n_arms()

------------------------------------------------------------------------

### Method `get_unique_times()`

Get unique timepoints.

#### Usage

    Timer$get_unique_times()

#### Returns

`numeric` vector of unique times.

#### Examples

    t <- Timer$new(name = "Timer")
    t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
    t$add_timepoint(time = 3.28, arm = "B", dropper = 6L, enroller = 23L)
    t$get_unique_times()

------------------------------------------------------------------------

### Method `get_timepoint()`

Get a timepoint by arm and index.

#### Usage

    Timer$get_timepoint(arm, i)

#### Arguments

- `arm`:

  `character` Arm identifier.

- `i`:

  `integer` Timepoint index.

#### Returns

`list` timepoint or `NULL` if not found.

#### Examples

    t <- Timer$new(name = "Timer")
    t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
    t$add_timepoint(time = 3.28, arm = "B", dropper = 6L, enroller = 23L)

    t$get_timepoint("A", 1)

------------------------------------------------------------------------

### Method `clone()`

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
t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 10L)
t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 12L)
t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 8L)

# Query
t$get_end_timepoint() # max time => 2
#> [1] 2
t$get_n_arms()        # unique arms => 2
#> [1] 2
t$get_unique_times()  # unique times => c(1, 2)
#> [1] 1 2
t$get_timepoint("A", 1) # returns a single timepoint
#> $time
#> [1] 1
#> 
#> $arm
#> [1] "A"
#> 
#> $dropper
#> [1] 2
#> 
#> $enroller
#> [1] 10
#> 


## ------------------------------------------------
## Method `Timer$new`
## ------------------------------------------------

t <- Timer$new(name = "Timer")

## ------------------------------------------------
## Method `Timer$add_timepoint`
## ------------------------------------------------

t <- Timer$new(name = "Timer")
t$add_timepoint(
  time = 1,
  arm = "A",
  dropper = 1L,
  enroller = 3L
)

## ------------------------------------------------
## Method `Timer$get_end_timepoint`
## ------------------------------------------------

t <- Timer$new(name = "Timer")
t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
t$get_end_timepoint()
#> [1] 3.14

## ------------------------------------------------
## Method `Timer$get_n_arms`
## ------------------------------------------------

t <- Timer$new(name = "Timer")
t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
t$add_timepoint(time = 3.28, arm = "B", dropper = 6L, enroller = 23L)
t$get_n_arms()
#> [1] 2

## ------------------------------------------------
## Method `Timer$get_unique_times`
## ------------------------------------------------

t <- Timer$new(name = "Timer")
t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
t$add_timepoint(time = 3.28, arm = "B", dropper = 6L, enroller = 23L)
t$get_unique_times()
#> [1] 3.14 3.28

## ------------------------------------------------
## Method `Timer$get_timepoint`
## ------------------------------------------------

t <- Timer$new(name = "Timer")
t$add_timepoint(time = 3.14, arm = "A", dropper = 7L, enroller = 22L)
t$add_timepoint(time = 3.28, arm = "B", dropper = 6L, enroller = 23L)

t$get_timepoint("A", 1)
#> NULL
```

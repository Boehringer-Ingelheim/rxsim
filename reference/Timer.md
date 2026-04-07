# Timer: Track timed events and apply condition-triggered analyses

A class to collect and query *timepoints*, time-based events, across
arms. Timer class also supports conditions that filter data using
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
and apply custom analyses.

Use `add_timepoint()` to append timepoints, `get_timepoint()` for a
lookup, and `check_conditions()` to filter a data frame based on a
trigger condition and return either analysis results or the filtered
data.

## Details

Helper functions
[`trigger_by_calendar()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_calendar.md)
and
[`trigger_by_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_fraction.md)
provide convenient shortcuts for common trigger patterns.

## See also

[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)
to coordinate simulations with populations,
[`add_timepoints()`](https://boehringer-ingelheim.github.io/rxsim/reference/add_timepoints.md)
to attach multiple timepoints,
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
for condition syntax.

## Public fields

- `name`:

  `character` Unique identifier for the `Timer` instance.

- `timelist`:

  `list` A list of timepoints. Each timepoint is a list with keys:

  - `time` `numeric` Calendar time

  - `arm` `character` Unique identifier of the arm

  - `dropper` `integer` \# of subjects dropper at `time`

  - `enroller` `integer` \# of subjects enrolled at `time`

- `conditions`:

  `list` A list of condition entries. Each entry is a list with keys:

  - `where` `expr` filter conditions in
    [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
    style

  - `analysis` `function` or `NULL` analysis applied to filtered data

  - `name` `character` or `NULL` unique key for the condition

  - `cooldown` `numeric` minimum time between consecutive triggers

  - `max_triggers` `integer` maximum number of times this condition can
    trigger

## Methods

### Public methods

- [`Timer$new()`](#method-Timer-new)

- [`Timer$add_timepoint()`](#method-Timer-add_timepoint)

- [`Timer$add_condition()`](#method-Timer-add_condition)

- [`Timer$get_end_timepoint()`](#method-Timer-get_end_timepoint)

- [`Timer$get_n_arms()`](#method-Timer-get_n_arms)

- [`Timer$get_unique_times()`](#method-Timer-get_unique_times)

- [`Timer$get_timepoint()`](#method-Timer-get_timepoint)

- [`Timer$check_conditions()`](#method-Timer-check_conditions)

- [`Timer$clone()`](#method-Timer-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `Timer` instance.

#### Usage

    Timer$new(name, timelist = NULL, conditions = NULL)

#### Arguments

- `name`:

  `character` Unique identifier.

- `timelist`:

  `list` Optional list of timepoints.

- `conditions`:

  `list` Optional list of condition entries.

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

### Method `add_condition()`

Add a trigger condition to a timer.

#### Usage

    Timer$add_condition(
      ...,
      analysis = NULL,
      name = NULL,
      cooldown = 0,
      max_triggers = 1L
    )

#### Arguments

- `...`:

  `expression` Boolean expression(s) for
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).

- `analysis`:

  `function` or `NULL` Optional function to apply.

- `name`:

  `character` Unique condition identifier.

- `cooldown`:

  `numeric` Minimum time between consecutive triggers (default: 0, no
  cooldown).

- `max_triggers`:

  `integer` Maximum number of times this condition can trigger (default:
  1, single trigger).

#### Examples

    #' t <- Timer$new(name = "Timer")

    # Add timepoints
    t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 10L)

    # Add conditions using `dplyr` style
    # Suppose you have a data.frame:
    df <- data.frame(
      id = 1:6,
      arm = c("A","A","B","B","A","B"),
      status = c("active","inactive","active","active","inactive","active"),
      visit = c(1,2,1,3,3,2)
    )

    # Analysis function: count rows at/after a given visit, per arm
    my_analysis <- function(dat, current_time) {
      out <- aggregate(id ~ arm, dat, length)
      out$current_time <- current_time
      out
    }

    # Condition 1: active only
    t$add_condition(
      status == "active",
      analysis = my_analysis,
      name = "active_only"
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

### Method `check_conditions()`

Check conditions and return filtered data or analysis results.

#### Usage

    Timer$check_conditions(locked_data, current_time)

#### Arguments

- `locked_data`:

  `data.frame` Trial data.

- `current_time`:

  `numeric` Calendar time.

#### Returns

`list` of filtered data or analysis results per condition.

#### Examples

    #' t <- Timer$new(name = "Timer")

    # Add timepoints
    t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 10L)
    t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 12L)
    t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 8L)

    # Query
    t$get_end_timepoint()     # max time => 2
    t$get_n_arms()            # unique arms => 2
    t$get_unique_times()      # unique times => c(1, 2)
    t$get_timepoint("A", 1)   # returns a single timepoint

    # Add conditions using dplyr style
    # Suppose you have a data.frame:
    df <- data.frame(
      id = 1:6,
      arm = c("A","A","B","B","A","B"),
      status = c("active", "inactive", "active", "active", "inactive", "active"),
      visit = c(1,2,1,3,3,2)
    )

    # Analysis function: count rows at/after a given visit, per arm
    my_analysis <- function(dat, current_time) {
      out <- aggregate(id ~ arm, dat, length)
      out$current_time <- current_time
      out
    }

    # Condition: active only
    t$add_condition(
      status == "active",
      analysis = my_analysis,
      name = "active_only"
    )

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
t$get_n_arms() # unique arms => 2
#> [1] 2
t$get_unique_times() # unique times => c(1, 2)
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

# Add conditions using trigger helpers or dplyr style
# Suppose you have a data.frame:
df <- data.frame(
  id = 1:6,
  arm = c("A", "A", "B", "B", "A", "B"),
  status = c("active", "inactive", "active", "active", "inactive", "active"),
  visit = c(1, 2, 1, 3, 3, 2)
)

# Analysis function: count rows at/after a given visit, per arm
my_analysis <- function(dat, current_time) {
  out <- aggregate(id ~ arm, dat, length)
  out$current_time <- current_time
  out
}

# Or add conditions manually with dplyr style
# Condition: arm A, visit >= 2, no analysis -> returns filtered df
t$add_condition(
  arm == "A", visit >= 2,
  name = "armA_visit2plus"
)

# Run checks
res <- t$check_conditions(locked_data = df, current_time = 3)
#> Warning:  returning filtered data as is because condition 'armA_visit2plus' has no applicable analysis 
names(res)
#> [1] "armA_visit2plus"


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
## Method `Timer$add_condition`
## ------------------------------------------------

#' t <- Timer$new(name = "Timer")

# Add timepoints
t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 10L)

# Add conditions using `dplyr` style
# Suppose you have a data.frame:
df <- data.frame(
  id = 1:6,
  arm = c("A","A","B","B","A","B"),
  status = c("active","inactive","active","active","inactive","active"),
  visit = c(1,2,1,3,3,2)
)

# Analysis function: count rows at/after a given visit, per arm
my_analysis <- function(dat, current_time) {
  out <- aggregate(id ~ arm, dat, length)
  out$current_time <- current_time
  out
}

# Condition 1: active only
t$add_condition(
  status == "active",
  analysis = my_analysis,
  name = "active_only"
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

## ------------------------------------------------
## Method `Timer$check_conditions`
## ------------------------------------------------

#' t <- Timer$new(name = "Timer")

# Add timepoints
t$add_timepoint(time = 1, arm = "A", dropper = 2L, enroller = 10L)
t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 12L)
t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 8L)

# Query
t$get_end_timepoint()     # max time => 2
#> [1] 3.28
t$get_n_arms()            # unique arms => 2
#> [1] 2
t$get_unique_times()      # unique times => c(1, 2)
#> [1] 3.14 3.28 1.00 2.00
t$get_timepoint("A", 1)   # returns a single timepoint
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

# Add conditions using dplyr style
# Suppose you have a data.frame:
df <- data.frame(
  id = 1:6,
  arm = c("A","A","B","B","A","B"),
  status = c("active", "inactive", "active", "active", "inactive", "active"),
  visit = c(1,2,1,3,3,2)
)

# Analysis function: count rows at/after a given visit, per arm
my_analysis <- function(dat, current_time) {
  out <- aggregate(id ~ arm, dat, length)
  out$current_time <- current_time
  out
}

# Condition: active only
t$add_condition(
  status == "active",
  analysis = my_analysis,
  name = "active_only"
)
```

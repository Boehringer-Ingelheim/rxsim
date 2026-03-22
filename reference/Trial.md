# Trial: Simulate a multi‑arm clinical trial

The `Trial` class coordinates one or more `Population` objects and a
`Timer` to simulate a clinical trial.

At each unique time defined in the trial's `Timer`, the `Trial`:

- applies enrollment and dropout updates to each `Population`

- builds a snapshot of all currently enrolled subjects

- evaluates all conditions in the `Timer`

- stores both the snapshot `locked_data` and the analysis outputs
  `results`

Use `run()` to execute the simulation. Trigger conditions are best added
with helper functions
[`trigger_by_calendar()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_calendar.md)
or
[`trigger_by_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/trigger_by_fraction.md).

## See also

[Population](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md),
[Timer](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md),
[`prettify_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/prettify_results.md),
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md),
[`clone_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/clone_trial.md).

[Timer](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md),
[`prettify_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/prettify_results.md).

## Public fields

- `name`:

  `character` Unique trial identifier.

- `seed`:

  `numeric` or `NULL` Random seed for reproducibility.

- `timer`:

  `Timer` object with timepoints and conditions.

- `population`:

  `list` of
  [Population](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md)
  objects, one per arm.

- `locked_data`:

  `list` Snapshots at each timepoint.

- `results`:

  `list` Analysis outputs per condition.

## Methods

### Public methods

- [`Trial$new()`](#method-Trial-new)

- [`Trial$run()`](#method-Trial-run)

- [`Trial$clone()`](#method-Trial-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `Trial` instance.

#### Usage

    Trial$new(
      name,
      seed = NULL,
      timer = NULL,
      population = list(),
      locked_data = list(),
      results = list()
    )

#### Arguments

- `name`:

  `character` Unique identifier for the trial.

- `seed`:

  `numeric` or `NULL` Optional random seed for reproducibility.

- `timer`:

  `Timer` object defining timepoints and conditions.

- `population`:

  `list` of
  [Population](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md)
  objects, one per arm.

- `locked_data`:

  `list` Generated at each `$run()` call.

- `results`:

  `list` Analysis outputs generated at each `$run()` call.

#### Returns

A new `Trial` instance.

#### Examples

    t <- Timer$new(name="simple_timer")
    pop <- Population$new(
      name = "simple_pop",
      data = vector_to_dataframe(rnorm(5))
    )
    pop$set_enrolled(5, 1)
    Trial$new(name = "simple_trial", timer=t, population = list(pop))

------------------------------------------------------------------------

### Method `run()`

Execute a trial simulation.

At each unique time defined by the trial's `Timer`:

- Apply enrollment and dropout actions to each `Population`

- Build a combined snapshot of all currently enrolled subjects

- Attach a `time` column to the snapshot

- Evaluate all condition readers via `Timer$check_conditions()`

- Store snapshots and condition outputs under time‑indexed list keys

#### Usage

    Trial$run()

#### Returns

Updates `locked_data` and `results` fields.

#### Examples

    # Create two populations
    popA <- Population$new("A", data = vector_to_dataframe(rnorm(10)))
    popB <- Population$new("B", data = vector_to_dataframe(rnorm(12)))

    # Create a timer and add timepoints
    t <- Timer$new("Timer")
    t$add_timepoint(time = 1, arm = "A", dropper = 0L, enroller = 4L)
    t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 5L)
    t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 2L)
    t$add_timepoint(time = 2, arm = "B", dropper = 2L, enroller = 3L)

    # Create a trial
    trial <- Trial$new(
      name = "ExampleTrial",
      seed = 123,
      timer = t,
      population = list(popA, popB)
    )

    # Run the simulation
    trial$run()

    prettify_results(trial$results)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Trial$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# Create two populations
popA <- Population$new("A", data = vector_to_dataframe(rnorm(10)))
popB <- Population$new("B", data = vector_to_dataframe(rnorm(12)))

# Create a timer and add timepoints
t <- Timer$new("Timer")
t$add_timepoint(time = 1, arm = "A", dropper = 0L, enroller = 4L)
t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 5L)
t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 2L)
t$add_timepoint(time = 2, arm = "B", dropper = 2L, enroller = 3L)

# Add trigger for final analysis at time 2
trigger_by_calendar(2, t, analysis = function(df, current_time) {
  nrow(df)
})

# Create a trial
trial <- Trial$new(
  name = "ExampleTrial",
  seed = 123,
  timer = t,
  population = list(popA, popB)
)

# Run the simulation
trial$run()

prettify_results(trial$results)
#>   time cal_time_2
#> 1    2         14


## ------------------------------------------------
## Method `Trial$new`
## ------------------------------------------------

t <- Timer$new(name="simple_timer")
pop <- Population$new(
  name = "simple_pop",
  data = vector_to_dataframe(rnorm(5))
)
pop$set_enrolled(5, 1)
Trial$new(name = "simple_trial", timer=t, population = list(pop))
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     initialize: function (name, seed = NULL, timer = NULL, population = list(), 
#>     locked_data: list
#>     name: simple_trial
#>     population: list
#>     results: list
#>     run: function () 
#>     seed: NULL
#>     timer: Timer, R6

## ------------------------------------------------
## Method `Trial$run`
## ------------------------------------------------

# Create two populations
popA <- Population$new("A", data = vector_to_dataframe(rnorm(10)))
popB <- Population$new("B", data = vector_to_dataframe(rnorm(12)))

# Create a timer and add timepoints
t <- Timer$new("Timer")
t$add_timepoint(time = 1, arm = "A", dropper = 0L, enroller = 4L)
t$add_timepoint(time = 1, arm = "B", dropper = 0L, enroller = 5L)
t$add_timepoint(time = 2, arm = "A", dropper = 1L, enroller = 2L)
t$add_timepoint(time = 2, arm = "B", dropper = 2L, enroller = 3L)

# Create a trial
trial <- Trial$new(
  name = "ExampleTrial",
  seed = 123,
  timer = t,
  population = list(popA, popB)
)

# Run the simulation
trial$run()

prettify_results(trial$results)
#> NULL
```

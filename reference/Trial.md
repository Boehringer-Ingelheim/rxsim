# Trial: Simulate a multi‑arm clinical trial

The `Trial` class coordinates one or more `Population` objects, a
`Timer`, and a list of `Condition` objects to simulate a clinical trial.

At each unique time defined in the trial's `Timer`, the `Trial`:

- applies enrollment and dropout updates to each `Population`

- builds a snapshot of all currently enrolled subjects

- evaluates each
  [`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
  in `self$conditions` against the snapshot

- stores both the snapshot (`locked_data`) and the analysis outputs
  (`results`)

Use `run()` to execute the simulation. Trigger conditions are built with
[`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)`$new()`
(or helpers
[`condition_calendar_time()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_calendar_time.md)
/
[`condition_enrollment_fraction()`](https://boehringer-ingelheim.github.io/rxsim/reference/condition_enrollment_fraction.md))
and stored in `trial$conditions`.

## See also

[Population](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md),
[Timer](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md),
[Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md),
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md),
[`replicate_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/replicate_trial.md),
[`clone_trial()`](https://boehringer-ingelheim.github.io/rxsim/reference/clone_trial.md).

[Timer](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md),
[Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md),
[`collect_results()`](https://boehringer-ingelheim.github.io/rxsim/reference/collect_results.md).

## Public fields

- `name`:

  `character` Unique trial identifier.

- `seed`:

  `numeric` or `NULL` Random seed for reproducibility.

- `timer`:

  `Timer` object with timepoints.

- `population`:

  `list` of
  [Population](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md)
  objects, one per arm.

- `conditions`:

  `list` of
  [Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
  objects evaluated at each timepoint.

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
      conditions = list(),
      results = list()
    )

#### Arguments

- `name`:

  `character` Unique identifier for the trial.

- `seed`:

  `numeric` or `NULL` Optional random seed for reproducibility.

- `timer`:

  `Timer` object defining timepoints.

- `population`:

  `list` of
  [Population](https://boehringer-ingelheim.github.io/rxsim/reference/Population.md)
  objects, one per arm.

- `locked_data`:

  `list` Generated at each `$run()` call.

- `conditions`:

  `list` of
  [Condition](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
  objects to evaluate at each timepoint.

- `results`:

  `list` Analysis outputs generated at each `$run()` call.

#### Returns

A new `Trial` instance.

#### Examples

    t <- Timer$new(name="simple_timer")
    pop <- Population$new(
      name = "simple_pop",
      data = as_population_data(rnorm(5))
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

- Evaluate each
  [`Condition`](https://boehringer-ingelheim.github.io/rxsim/reference/Condition.md)
  in `self$conditions` against the snapshot

- Store snapshots and condition outputs under time‑indexed list keys

#### Usage

    Trial$run()

#### Returns

Updates `locked_data` and `results` fields.

#### Examples

    # Create two populations
    popA <- Population$new("A", data = as_population_data(rnorm(10)))
    popB <- Population$new("B", data = as_population_data(rnorm(12)))

    # Create a timer and add timepoints
    t <- Timer$new("Timer")
    t$add_timepoint(time = 1, arm = "A", drop = 0L, enroll = 4L)
    t$add_timepoint(time = 1, arm = "B", drop = 0L, enroll = 5L)
    t$add_timepoint(time = 2, arm = "A", drop = 1L, enroll = 2L)
    t$add_timepoint(time = 2, arm = "B", drop = 2L, enroll = 3L)

    # Create a trial
    trial <- Trial$new(
      name = "ExampleTrial",
      seed = 123,
      timer = t,
      population = list(popA, popB)
    )

    # Run the simulation
    trial$run()

    collect_results(trial)

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
popA <- Population$new("A", data = as_population_data(rnorm(10)))
popB <- Population$new("B", data = as_population_data(rnorm(12)))

# Create a timer and add timepoints
t <- Timer$new("Timer")
t$add_timepoint(time = 1, arm = "A", drop = 0L, enroll = 4L)
t$add_timepoint(time = 1, arm = "B", drop = 0L, enroll = 5L)
t$add_timepoint(time = 2, arm = "A", drop = 1L, enroll = 2L)
t$add_timepoint(time = 2, arm = "B", drop = 2L, enroll = 3L)

# Build a condition: fire at time >= 2 and count enrolled rows
cond <- Condition$new(
  where    = calendar_trigger(2),
  analysis = function(df, current_time) nrow(df),
  name     = "final"
)

# Create a trial
trial <- Trial$new(
  name       = "ExampleTrial",
  seed       = 123,
  timer      = t,
  population = list(popA, popB),
  conditions = list(cond)
)

# Run the simulation
trial$run()

collect_results(trial)
#>   replicate timepoint analysis X14L
#> 1         1         2    final   14


## ------------------------------------------------
## Method `Trial$new`
## ------------------------------------------------

t <- Timer$new(name="simple_timer")
pop <- Population$new(
  name = "simple_pop",
  data = as_population_data(rnorm(5))
)
pop$set_enrolled(5, 1)
Trial$new(name = "simple_trial", timer=t, population = list(pop))
#> <Trial>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     conditions: list
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
popA <- Population$new("A", data = as_population_data(rnorm(10)))
popB <- Population$new("B", data = as_population_data(rnorm(12)))

# Create a timer and add timepoints
t <- Timer$new("Timer")
t$add_timepoint(time = 1, arm = "A", drop = 0L, enroll = 4L)
t$add_timepoint(time = 1, arm = "B", drop = 0L, enroll = 5L)
t$add_timepoint(time = 2, arm = "A", drop = 1L, enroll = 2L)
t$add_timepoint(time = 2, arm = "B", drop = 2L, enroll = 3L)

# Create a trial
trial <- Trial$new(
  name = "ExampleTrial",
  seed = 123,
  timer = t,
  population = list(popA, popB)
)

# Run the simulation
trial$run()

collect_results(trial)
#> # A tibble: 0 × 0
```

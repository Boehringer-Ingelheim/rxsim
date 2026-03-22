# Population: Manage a patient population

The `Population` class stores subject-level data and manages enrollment
and dropout times for each subject.

Use `set_enrolled()` and `set_dropped()` to assign enrollment and
dropout times to random subsets of subjects.

## See also

[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md)
for multi-arm simulations,
[Timer](https://boehringer-ingelheim.github.io/rxsim/reference/Timer.md)
for managing timepoints.

Population,
[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md).

Population,
[Trial](https://boehringer-ingelheim.github.io/rxsim/reference/Trial.md).

Population.

## Public fields

- `name`:

  `character` Unique identifier for the population.

- `data`:

  `data.frame` Subject-level data frame with columns:

  - `id` `integer`

  - `arm` `character`

  - `readout_time` `numeric`

  - `data` `numeric`

  - may contain more columns

- `enrolled`:

  `numeric` Vector of enrollment times for each subject.

- `dropped`:

  `numeric` Vector of dropout times for each subject.

- `n`:

  `integer` Number of unique subjects in the population.

- `n_readouts`:

  `integer` Number of readout_times in the population.

## Methods

### Public methods

- [`Population$new()`](#method-Population-new)

- [`Population$set_enrolled()`](#method-Population-set_enrolled)

- [`Population$set_dropped()`](#method-Population-set_dropped)

- [`Population$set_data()`](#method-Population-set_data)

- [`Population$clone()`](#method-Population-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `Population` instance.

#### Usage

    Population$new(
      name,
      data = NULL,
      enrolled = NULL,
      dropped = NULL,
      n = NULL,
      n_readouts = NULL
    )

#### Arguments

- `name`:

  `character` Unique identifier for the population.

- `data`:

  `data.frame` with columns: `id`, `arm`, `readout_time`, `data`, and
  optionally more columns.

- `enrolled`:

  `numeric` Optional enrollment times (auto-initialized if `NULL`).

- `dropped`:

  `numeric` Optional dropout times (auto-initialized if `NULL`).

- `n`:

  `integer` Auto-computed from data (optional).

- `n_readouts`:

  `integer` Auto-computed from data (optional).

#### Returns

A new `Population` instance.

#### Examples

    Population$new(name = "Intervention", data = vector_to_dataframe(rnorm(5)))

------------------------------------------------------------------------

### Method `set_enrolled()`

Mark subjects as enrolled at a given time.

Enrollment applies only to unenrolled subjects (`NA`).

#### Usage

    Population$set_enrolled(n, time)

#### Arguments

- `n`:

  `integer` Number of subjects to enroll.

- `time`:

  `numeric` Enrollment time.

#### Examples

    pop <- Population$new("Test", vector_to_dataframe(rnorm(10)))
    pop$set_enrolled(n = 4, time = 2)

------------------------------------------------------------------------

### Method `set_dropped()`

Mark subjects as dropped at a given time.

Dropout applies only to enrolled, not-yet-dropped subjects.

#### Usage

    Population$set_dropped(n, time)

#### Arguments

- `n`:

  `integer` Number of subjects to drop.

- `time`:

  `numeric` Dropout time.

#### Examples

    pop <- Population$new("Test", vector_to_dataframe(rnorm(10)))
    pop$set_enrolled(n = 5, time = 1)
    pop$set_dropped(n = 2, time = 3)

------------------------------------------------------------------------

### Method `set_data()`

Replace underlying subject data and reset enrollment/dropout status.

#### Usage

    Population$set_data(data)

#### Arguments

- `data`:

  `data.frame` with columns: `id`, `data`, `arm`, `readout_time`, and
  optionally more columns.

#### Examples

    pop <- Population$new("ResetDemo", vector_to_dataframe(rnorm(5)))
    pop$set_data(
      data.frame(
        id = 1:8,
        data = rnorm(8),
        arm = "ResetDemo",
        readout_time = 0
      )
    )

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Population$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# Basic example: vector input
pop <- Population$new(name = "Control", data = vector_to_dataframe(rnorm(10)))

pop$n # number of subjects
#> [1] 10
head(pop$data) # generated subject-level data
#>   id         data readout_time     arm
#> 1  1 -1.400043517            0 Control
#> 2  2  0.255317055            0 Control
#> 3  3 -2.437263611            0 Control
#> 4  4 -0.005571287            0 Control
#> 5  5  0.621552721            0 Control
#> 6  6  1.148411606            0 Control

# Set enrollment for 5 subjects at time = 1
pop$set_enrolled(n = 5, time = 1)

# Drop 2 subjects at time = 3
pop$set_dropped(n = 2, time = 3)

# Reset underlying data
pop$set_data(vector_to_dataframe(rnorm(8)))


## ------------------------------------------------
## Method `Population$new`
## ------------------------------------------------

Population$new(name = "Intervention", data = vector_to_dataframe(rnorm(5)))
#> <Population>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     data: data.frame
#>     dropped: NA NA NA NA NA
#>     enrolled: NA NA NA NA NA
#>     initialize: function (name, data = NULL, enrolled = NULL, dropped = NULL, 
#>     n: 5
#>     n_readouts: 1
#>     name: Intervention
#>     set_data: function (data) 
#>     set_dropped: function (n, time) 
#>     set_enrolled: function (n, time) 

## ------------------------------------------------
## Method `Population$set_enrolled`
## ------------------------------------------------

pop <- Population$new("Test", vector_to_dataframe(rnorm(10)))
pop$set_enrolled(n = 4, time = 2)

## ------------------------------------------------
## Method `Population$set_dropped`
## ------------------------------------------------

pop <- Population$new("Test", vector_to_dataframe(rnorm(10)))
pop$set_enrolled(n = 5, time = 1)
pop$set_dropped(n = 2, time = 3)

## ------------------------------------------------
## Method `Population$set_data`
## ------------------------------------------------

pop <- Population$new("ResetDemo", vector_to_dataframe(rnorm(5)))
pop$set_data(
  data.frame(
    id = 1:8,
    data = rnorm(8),
    arm = "ResetDemo",
    readout_time = 0
  )
)
```

# rxsim 0.3.0

This release is mostly about paying down debt. The timer got a simpler internal
representation, fixed designs got a faster execution path, and two long-deprecated
helpers finally came out. If you are upgrading from 0.2.0, the removals below are
the only changes that should require edits on your side.

## Breaking changes

* `add_timepoints()` and `Timer$get_timepoint()` have been removed (#107). Both
  predate `Timer$add_schedule()`, which does the same job in one call and accepts
  the output of `stochastic_schedule()` and `deterministic_schedule()` without
  reshaping.

  ```r
  # before
  add_timepoints(timer, plan)

  # after
  timer$add_schedule(plan)
  ```

  If you built timepoints by hand, pass a data frame with `time`, `arm`,
  `enroll`, and `drop` columns straight to `add_schedule()`. Note that `enroll`
  and `drop` are subject counts and must be integers (`3L`, not `3`).

* `Timer$timelist` is now a plain data frame rather than a list of lists (#107).
  Code that reassembled it with `do.call(rbind, timer$timelist)` can read the
  field directly. Anything that indexed into it positionally will need updating.

* `stochastic_schedule()` and `deterministic_schedule()` are now genuine drop-in
  replacements for each other. Both return a plain `data.frame` sorted by `arm`
  and then `time`; previously `deterministic_schedule()` returned a tibble with a
  different row order, which meant swapping one for the other could quietly change
  downstream behaviour.

## New features

* Fixed designs (`adaptive = FALSE`, the default on `Trial$new()` and
  `replicate_trial()`) now run through a dedicated fast path that precomputes
  enrollment and dropout times and takes cheap prefix slices for each snapshot
  (#108). A companion optimisation skips leading timepoints where no condition
  can possibly fire (#109). Results are identical to the adaptive loop; only the
  runtime differs.

## Bug fixes

* `deterministic_schedule()` no longer fails when a sample-size correction lands
  before the first dropout period boundary. `findInterval()` returns zero in that
  case, and the resulting empty vector crashed the schedule assembly. Corrections
  are now attributed to the first period, which is what the right-closed period
  semantics imply.

* Trigger expressions combined with `|` now handle `NA` guards correctly, and
  quosure generation for composed triggers was tightened up (#106).

* `n_readouts > 1` is now handled consistently across both execution paths.

## Documentation

* A thorough accuracy pass over the README and all five vignettes. Several column
  tables, code examples, and cross-references had drifted behind earlier API
  changes; these now match the shipped behaviour. The README quick-start example
  and the assembly walkthrough in *Core Concepts* are executed at build time, so
  they cannot silently rot again.

* The trigger reference is complete: `notna_trigger()`, `col_trigger()`,
  `timed_count_trigger()`, and `trigger_by_events()` are now documented alongside
  the primitives they complement, and the fire-once default (`max_triggers = 1L`)
  is stated explicitly with a runnable example.

# rxsim 0.2.0.9000

## New features

* New event-driven trigger, `trigger_by_events()`, fires an analysis once a target
  number of time-to-event endpoints has been observed — the natural way to schedule
  looks in oncology and other event-driven settings (#99). It is built from three
  new lower-level primitives that are also exported for custom rules:
  `notna_trigger()`, `col_trigger()`, and `timed_count_trigger()`.

## Improvements

* Added statistical ground-truth tests covering the stochastic assignment paths,
  plus verification that the adaptive and fixed execution paths produce identical
  results.
* Unit test coverage extended to the helper functions.

# rxsim 0.1.7.9000

## Breaking changes

* The five worked examples have moved to a separate
  [rxsim-gallery](https://boehringer-ingelheim.github.io/rxsim-gallery/) site
  (#103). They were the only reason rxsim suggested `RBesT`, `DoseFinding`,
  `BayesianMCPMod`, `multcomp`, `survival`, `MASS`, and several others; dropping
  them removes a substantial installation burden from a package whose core has
  always been small. The examples themselves are unchanged and remain online.

# rxsim 0.1.6.9000

A large API consolidation. This is the release where the trigger system took its
current shape.

## Breaking changes

* `gen_plan()` and `gen_timepoints()` are now `stochastic_schedule()` and
  `deterministic_schedule()`. The new names say what the functions do rather than
  how they were implemented, and the signatures are aligned so the two can be
  swapped.
* `vector_to_dataframe()` is now `as_population_data()`.
* `prettify_results()` has been removed; `collect_results()` covers the same
  ground and returns a tidy data frame.
* `trigger_by_calendar()` and `trigger_by_fraction()` are replaced by
  `condition_calendar_time()` and `condition_enrollment_fraction()`.

## New features

* A new `Condition` class separates trigger logic from the `Timer`, which is now
  purely the trial clock. Conditions are passed to `Trial$new(conditions = ...)`.
* Composable trigger objects: `value_trigger()`, `count_trigger()`,
  `enroll_trigger()`, and `calendar_trigger()` can be combined with `&` and `|`
  to express rules such as "half enrolled *and* past month six" without writing
  filter expressions by hand.
* `Population` now validates that all arms declare the same number of readouts,
  catching a class of mismatch that previously surfaced as a confusing error
  deeper in the run loop.

## Improvements

* `collect_results()` rewritten for clarity and speed.
* Test coverage added for the `Trial` class.
* Vignettes reorganised around concepts rather than individual classes.

# rxsim 0.1.2.9000

## New features

* Conditions can now be capped so they do not fire repeatedly: `cooldown` sets a
  minimum time gap between triggers, and `max_triggers` sets a hard ceiling
  (#58). Useful for data-driven checks evaluated at every timepoint.

## Improvements

* A helper for building population data frames, and a polish pass across all
  worked examples and the getting-started guide.

# rxsim 0.1.0.9000

* First tagged release with the package website published via pkgdown, along with
  continuous integration running `R CMD check` on every push.

# rxsim 0.0.6.9000

* Initial version. Establishes the three building blocks the package is still
  organised around — `Population`, `Timer`, and `Trial` — together with
  `replicate_trial()` and `run_trials()` for generating and executing many
  independent replicates.

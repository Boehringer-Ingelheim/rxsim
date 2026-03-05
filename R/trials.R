#' Clone a Trial Object Multiple Times
#'
#' Creates deep copies of a `Trial` R6 object with independent
#' timer and population instances.
#'
#' @param trial `Trial` R6 object to clone.
#' @param n `integer` Number of clones to create. Defaults to `1`.
#'
#' @return `list` of `n` independently cloned `Trial` objects.
#'
#' @seealso [Trial], [replicate_trial()], [run_trials()].
#'
#' @examples
#' # Create 5 independent copies of a trial:
#' # clones <- clone_trial(trial, n = 5)
#'
#' @export
clone_trial <- function(trial, n = 1) {
  if (missing(trial) || !inherits(trial, "Trial")) stop("`trial` must be a Trial instance")
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) stop("`n` must be a single positive integer")

  lapply(seq_len(n), function(i) {
    Trial$new(
      name = paste(trial$name, i, sep="_"),
      timer = trial$timer$clone(),
      population = lapply(trial$population, function(x) x$clone())
    )
  })
}

#' Generate a Population Object
#'
#' Creates a `Population` R6 object by executing a generator function.
#'
#' @param name `character` Population name (e.g., arm label).
#' @param generator `function` that returns a data.frame-like object.
#' @param sample_size `integer` Number of subjects to generate.
#'
#' @return `Population` R6 object.
#'
#' @seealso [Population], [replicate_trial()], [vector_to_dataframe()].
#'
#' @examples
#' gen_control <- function(n) data.frame(id = 1:n, value = rnorm(n), readout_time = 0)
#' pop <- gen_population("control", gen_control, sample_size = 50)
#'
#' @export
gen_population <- function(name, generator, sample_size = 1) {
  if (!is.character(name) || length(name) != 1L) stop("`name` must be a single character value")
  if (!is.function(generator)) stop("`generator` must be a function")
  if (!is.numeric(sample_size) || length(sample_size) != 1L || sample_size < 0) stop("`sample_size` must be a single non-negative number")

  Population$new(
    name = name,
    data = generator(sample_size)
  )
}

#' Create Multiple Trials with Generated Populations
#'
#' Generates `n` independent `Trial` objects from a template,
#' timer, and population generators. Arms must match between plan and
#' population specifications.
#'
#' @param trial_name `character` Base name for generated trials (index appended).
#'   Defaults to `"name"`.
#' @param sample_size `integer` Total sample size across all arms.
#' @param arms `character` vector of arm identifiers.
#' @param allocation `numeric` vector of arm allocation ratios.
#' @param enrollment `function` that generates inter-enrollment times.
#' @param dropout `function` that generates inter-dropout times.
#' @param analysis_generators `list` (named) of analysis trigger specifications.
#' @param population_generators `list` (named) of population generator functions.
#' @param n `integer` Number of trials to create.
#'
#' @return `list` of `n` `Trial` objects with indexed names,
#'   cloned timers, and generated populations.
#'
#' @seealso [gen_plan()], [gen_population()], [clone_trial()], [run_trials()].
#'
#' @export
replicate_trial <- function(
  trial_name = "name",
  sample_size,
  arms,
  allocation,
  enrollment,
  dropout,
  analysis_generators,
  population_generators,
  n
) {

  timers <- lapply(seq_len(n), function(i) {
    t <- Timer$new(name = paste("timer", i, sep="_"))
    plan <- gen_plan(sample_size, arms, allocation, enrollment, dropout)
    add_timepoints(t, plan)
    lapply(names(analysis_generators), function(name) {
      t$add_condition(
        !!!analysis_generators[[name]]$trigger,
        analysis = analysis_generators[[name]]$analysis,
        name = name
      )
    })
    return(t)
  })

  pop_names <- names(population_generators)
  if (is.null(pop_names) || any(pop_names == "")) {
    stop("`population_generators` must be a named list.")
  }

  n_target <- lapply(timers, function(t) {
    plan_df <- dplyr::bind_rows(t$timelist)
    planned_arms <- unique(plan_df$arm)
    plan_missing_arms <- setdiff(pop_names, planned_arms)
    pop_missing_arms <- setdiff(planned_arms, pop_names)
    if (length(plan_missing_arms) > 0L) {
      stop(sprintf(
        "Population generator names not found in plan arms: %s",
        paste(plan_missing_arms, collapse = ", ")
      ))
    }
    if (length(pop_missing_arms) > 0L) {
      stop(sprintf(
        "Plan contains arms not found in population_generators: %s",
        paste(pop_missing_arms, collapse = ", ")
      ))
    }
    vapply(
      pop_names,
      function(a) as.integer(sum(plan_df$enroller[plan_df$arm == a], na.rm = TRUE)),
      integer(1)
    )
  })

  bad_trials <- which(vapply(n_target, function(x) sum(x) != sample_size, logical(1)))
  if (length(bad_trials) > 0L) {
    stop(sprintf(
      "Planned enrollment mismatch in trial(s): %s. Expected total %s enrolled subjects.",
      paste(bad_trials, collapse = ", "),
      as.character(sample_size)
    ))
  }

  trials <- lapply(seq_len(n), function(i) {
    Trial$new(
      name = paste(trial_name, i, sep="_"),
      timer = timers[[i]],
      population = lapply(names(population_generators), function(name) {
        gen_population(name, population_generators[[name]], n_target[[i]][[name]])
      })
    )
  })

  for (i in seq_along(trials)) {
    pop_actual <- vapply(trials[[i]]$population, function(p) as.integer(p$n), integer(1))
    names(pop_actual) <- vapply(trials[[i]]$population, function(p) p$name, character(1))

    pop_actual <- pop_actual[pop_names]
    pop_plan <- as.integer(n_target[[i]][pop_names])
    names(pop_plan) <- pop_names

    if (!identical(pop_actual, pop_plan)) {
      mismatch <- pop_names[pop_actual != pop_plan]
      detail <- paste(sprintf("%s planned=%s actual=%s", mismatch, pop_plan[mismatch], pop_actual[mismatch]), collapse = "; ")
      stop(sprintf("Population size mismatch in trial %s: %s", i, detail))
    }
  }

  return(trials)
}

#' Run Multiple Trial Objects
#'
#' Executes the `$run()` method for each trial in a list.
#'
#' @param trials `list` of `Trial` R6 objects.
#'
#' @return `list` of results from each trial's `$run()` method.
#'
#' @seealso [Trial], [replicate_trial()].
#'
#' @examples
#' # results <- run_trials(trials)
#'
#' @export


run_trials <- function(trials) {
  if (!is.list(trials)) stop("`trials` must be a list of Trial objects")
  lapply(trials, function(trial) trial$run())
}

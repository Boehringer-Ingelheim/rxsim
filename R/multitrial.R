# ======================================================================
# Trial replication and generation utilities
# ======================================================================

#' Replicate a Trial Object Multiple Times
#'
#' Creates `n` independent copies of a given `Trial` object. Each replicated trial receives:
#' - a new name (original name + index),
#' - a cloned timer object,
#' - a cloned population (each element cloned individually).
#'
#' This is useful for simulation studies where identical trial structures
#' must be evaluated independently.
#'
#' @param trial A `Trial` R6 object with fields:
#'   - `name`: character scalar,
#'   - `timer`: R6 object with `$clone()` method,
#'   - `population`: list of R6 objects, each supporting `$clone()`.
#' @param n Integer (length 1) number of replications to generate. Must be >= 1.
#'
#' @return A named list of length `n` containing cloned `Trial` objects.
#'   List elements are named `paste0(trial$name, "_", seq_len(n))`.
#'
#' @examples
#' # Assume `trial` is a fully constructed Trial R6 object
#' # trials <- replicate_trial(trial, n = 2)
#'
#' @export
replicate_trial <- function(trial, n = 1) {
  # ---- argument validation ----
  stopifnot("`n` must be a single positive integer" = is.numeric(n) && length(n) == 1 && n >= 1)
  n <- as.integer(n)

  if (!inherits(trial, "Trial")) {
    stop("`trial` must inherit from class 'Trial'.")
  }
  if (!is.character(trial$name) || length(trial$name) != 1) {
    stop("`trial$name` must be a character scalar.")
  }
  if (!is.list(trial$population)) {
    stop("`trial$population` must be a list of R6 objects with $clone().")
  }
  if (!is.function(trial$timer$clone)) {
    stop("`trial$timer` must provide a $clone() method.")
  }

  # ---- preallocate and build ----
  trials <- vector("list", n)
  names(trials) <- paste0(trial$name, "_", seq_len(n))

  for (i in seq_len(n)) {
    trials[[i]] <- Trial$new(
      name       = names(trials)[i],
      timer      = trial$timer$clone(deep = TRUE),
      population = lapply(trial$population, function(x) {
        if (!is.function(x$clone)) stop("Each population element must support $clone().")
        x$clone(deep = TRUE)
      })
    )
  }

  trials
}


#' Generate a Population Object from a Named Generator
#'
#' Helper to create a `Population` R6 object by running a provided generator function.
#'
#' @param name Character scalar; the population name (e.g., arm label).
#' @param generator A function with no required arguments that returns a data.frame-like object.
#'
#' @return A `Population` R6 object constructed as `Population$new(name = name, data = generator())`.
#'
#' @examples
#' gen_control <- function() data.frame(id = 1:5, value = rnorm(5))
#' pop <- gen_population(name = "control", generator = gen_control)
#'
#' @export
gen_population <- function(name = NULL, generator = NULL) {
  if (!is.character(name) || length(name) != 1 || is.na(name)) {
    stop("`name` must be a non-missing character scalar.")
  }
  if (!is.function(generator)) {
    stop("`generator` must be a function that returns data.")
  }

  dat <- generator()
  if (is.null(dat)) stop("`generator()` returned NULL; expected data.")

  Population$new(
    name = name,
    data = dat
  )
}


#' Create Multiple Trials from Data Generators
#'
#' Constructs `n` independent `Trial` objects given a trial base name, a timer,
#' and a list of `(name, generator)` pairs for populations. Each trial gets:
#' - a unique name `paste0(trial_name, "_", i)`,
#' - a deep-cloned timer,
#' - freshly generated `Population` objects by calling each generator.
#'
#' @param trial_name Character scalar base name for all generated trials.
#' @param data_gen_list A list where each element is a length-2 list or vector:
#'   - `[[1]]`: population name (character),
#'   - `[[2]]`: generator function (returns data for that population).
#'   Example: `list( list("pbo", gen_control), list("trt", gen_trt) )`
#' @param timer An R6 `Timer` object with `$clone()` method.
#' @param n Integer (length 1) number of trials to create. Must be >= 1.
#'
#' @return A named list of length `n` with constructed `Trial` objects.
#'
#' @examples
#' # gen_control <- function() data.frame(subject_id = 1:50, value = rnorm(50))
#' # gen_trt     <- function() data.frame(subject_id = 1:50, value = rnorm(50, 0.2))
#' # pops_gen <- list(list("pbo", gen_control), list("trt", gen_trt))
#' # t <- Timer$new(name = "trial_timer")
#' # trials <- create_multiple_trials_gen_pop("trial", pops_gen, t, n = 10)
#'
#' @export
create_multiple_trials_gen_pop <- function(trial_name = "trial",
                                           data_gen_list = NULL,
                                           timer,
                                           n = 1) {
  # ---- validation ----
  if (!is.character(trial_name) || length(trial_name) != 1) {
    stop("`trial_name` must be a character scalar.")
  }
  stopifnot("`n` must be a single positive integer" = is.numeric(n) && length(n) == 1 && n >= 1)
  n <- as.integer(n)

  if (!is.list(data_gen_list) || length(data_gen_list) == 0) {
    stop("`data_gen_list` must be a non-empty list of (name, generator) pairs.")
  }
  if (!inherits(timer, "Timer") || !is.function(timer$clone)) {
    stop("`timer` must be a `Timer` R6 object with a $clone() method.")
  }

  # Validate data_gen_list structure
  for (i in seq_along(data_gen_list)) {
    elt <- data_gen_list[[i]]
    if (length(elt) < 2) stop("Each element of `data_gen_list` must have at least two items: name, generator.")
    if (!is.character(elt[[1]]) || length(elt[[1]]) != 1) {
      stop("Population name (first element) must be a character scalar.")
    }
    if (!is.function(elt[[2]])) {
      stop("Generator (second element) must be a function.")
    }
  }

  # ---- preallocate and build ----
  trials <- vector("list", n)
  names(trials) <- paste0(trial_name, "_", seq_len(n))

  for (i in seq_len(n)) {
    trials[[i]] <- Trial$new(
      name       = names(trials)[i],
      timer      = timer$clone(deep = TRUE),
      population = lapply(data_gen_list, function(x) {
        gen_population(name = x[[1]], generator = x[[2]])
      })
    )
  }

  trials
}


#' Run Multiple Trials
#'
#' Applies `$run()` to each `Trial` in a list. Returns any values produced by
#' `$run()` (if it returns a value), otherwise returns `invisible(NULL)` for each.
#'
#' @param trials A list of `Trial` objects.
#' @param simplify Logical; if `TRUE` and results are list-like of equal length/shape,
#'   attempts to simplify. Defaults to `FALSE`.
#'
#' @return A list of results (or `NULL`s) corresponding to each trial’s `$run()`.
#'
#' @examples
#' # results <- run_multiple_trial(trials)
#'
#' @export
run_multiple_trial <- function(trials = NULL, simplify = FALSE) {
  if (!is.list(trials) || length(trials) == 0) {
    stop("`trials` must be a non-empty list of Trial objects.")
  }
  for (i in seq_along(trials)) {
    if (!inherits(trials[[i]], "Trial") || !is.function(trials[[i]]$run)) {
      stop("Each element of `trials` must be a `Trial` R6 object with a $run() method.")
    }
  }

  res <- lapply(trials, function(x) x$run())

  if (isTRUE(simplify)) {
    # A conservative attempt to simplify if possible
    tryCatch({
      return(simplify2array(res))
    }, error = function(e) {
      # Fallback to list if simplification fails
      return(res)
    })
  }

  res
}

# ======================================================================
# Example usage (can be placed in a vignette or examples section)
# ======================================================================

#--- Configuration / Generators ---
sample_size <- 100
arms        <- c("pbo", "trt")
allocation  <- c(1, 1)
enrollment  <- list(end_time = c(4, 8, 12, 16), rate = c(5, 10, 15, 20))
dropout     <- list(end_time = c(4, 8, 12, 16), rate = c(1, 2, 4, 8))

timepoints <- gen_timepoints(
  sample_size = sample_size,
  arms        = arms,
  allocation  = allocation,
  enrollment  = enrollment,
  dropout     = dropout
)

t <- Timer$new(name = "trial_timer")
add_timepoints(t, timepoints)

t$add_condition(
  time %in% c(6, 12),
  analysis = function(df, time) {
    df_enrolled <- df |>
      dplyr::filter(!is.na(enroll_time))
    stats::t.test(value ~ arm, data = df_enrolled)$p.value
  },
  name = "t_test_final"
)

gen_control <- function() {
  data.frame(
    subject_id = 1:50,
    value = stats::rnorm(50)
  )
}

gen_trt <- function() {
  data.frame(
    subject_id = 1:50,
    value = stats::rnorm(50, mean = 0.2)
  )
}

pops_gen <- list(
  list("pbo", gen_control),
  list("trt", gen_trt)
)

# Create 100 trials and run them
trials <- create_multiple_trials_gen_pop(
  trial_name    = "trial",
  data_gen_list = pops_gen,
  timer         = t,
  n             = 100
)

results <- run_multiple_trial(trials)

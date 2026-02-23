#' Population: Manage a patient population
#'
#' @description
#' The `Population` class stores subject-level data and manages enrollment and
#' dropout times for each subject.
#'
#' Use `set_enrolled()` and `set_dropped()` to assign enrollment and dropout
#' times to random subsets of subjects.
#'
#' @examples
#' # Basic example: vector input
#' pop <- Population$new(name = "Control", data = vector_to_dataframe(rnorm(10)))
#'
#' pop$n              # number of subjects
#' head(pop$data)     # generated subject-level data
#'
#' # Set enrollment for 5 subjects at time = 1
#' pop$set_enrolled(n = 5, time = 1)
#'
#' # Drop 2 subjects at time = 3
#' pop$set_dropped(n = 2, time = 3)
#'
#' # Reset underlying data
#' pop$set_data(vector_to_dataframe(rnorm(8)))
#'
#' @export
Population <- R6::R6Class(
  classname = "Population",
  public = list(
    # --- fields ---
    #' @field name `character` Unique identifier for the population.
    name = NULL,

    #' @field data `data.frame` Subject-level data frame with columns:
    #' - `subject_id` `integer`
    #' - `data` `numeric`
    #' - `population_name` `character`
    #' - may contain more columns
    data = NULL,

    #' @field enrolled `numeric` Vector of enrollment times for each subject.
    enrolled = NULL,

    #' @field dropped `numeric` Vector of dropout times for each subject.
    dropped = NULL,

    #' @field n `integer` Number of unique subjects in the population.
    n=NULL,
    #' @field m `integer` Number of readout_times in the population.
    m=NULL,
    # --- constructor ---
    #' @description
    #' Create a new `Population`` instance.
    #'
    #' @param name `character` Unique identifier for the population.
    #'
    #' @param data `data.frame` Underlying subject data frame with columns
    #' `subject_id`, `data`, and `population_name`.
    #'
    #' @param enrolled `numeric` enrollment times are always initialized
    #' automatically as `NA`.
    #'
    #' @param dropped `numeric` dropout times are always initialized
    #' automatically as `NA`.
    #'
    #' @param n `integer` `n` is automatically computed from data.
    #'
    #' @returns A new `Population` instance.
    #'
    #' @examples
    #' Population$new(name = "Intervention", data = vector_to_dataframe(rnorm(5)))
    initialize = function(
    name,
    data = NULL,
    enrolled = NULL,
    dropped = NULL,
    n=NULL,
    m=NULL
    ) {
      stopifnot(is.character(name))
      self$name <- name
      if(!("arm" %in% names(data))){
        data$arm <- name
      }
      self$data <- data
      self$n <- length(unique(self$data$subject_id))
      self$m <- (self$n/dim(self$data)[1])

      self$enrolled <- rep(NA, self$n)
      self$dropped <- rep(NA, self$n)
    },

    # --- methods ---
    #' @description
    #' Mark `n` enrolled subjects at a given time.
    #'
    #' Enrollment applies to subjects who are not yet enrolled `NA`.
    #'
    #' @param n `integer` Number of subjects to enroll.
    #' @param time `numeric` Enrollment time.
    #'
    #' @examples
    #' pop <- Population$new("Test", vector_to_dataframe(rnorm(10)))
    #' pop$set_enrolled(n = 4, time = 2)
    set_enrolled = function(n, time) {
      idx <- which(is.na(self$enrolled))
      self$enrolled[sample(idx, n, replace = FALSE)] <- time
    },

    #' @description
    #' Mark `n` subjects as dropped at a given time.
    #'
    #' Dropout applies only to subjects who:
    #' - have enrolled, and
    #' - have not already dropped.
    #'
    #' @param n `integer` # of subjects to drop
    #' @param time `numeric` Dropout time
    #'
    #' @examples
    #' pop <- Population$new("Test", vector_to_dataframe(rnorm(10)))
    #' pop$set_enrolled(n = 5, time = 1)
    #' pop$set_dropped(n = 2, time = 3)
    set_dropped = function(n, time) {
      idx <- which(is.na(self$dropped) & !is.na(self$enrolled))
      self$dropped[sample(idx, n, replace = FALSE)] <- time
    },

    #' @description
    #' Replace the underlying subject data and reset enrollment/dropout status.
    #'
    #' @param data `data.frame` Subject-level data frame with columns:
    #' - `subject_id` `integer`
    #' - `data` `numeric`
    #' - `population_name` `character`
    #' - may contain more columns
    #'
    #' @examples
    #' pop <- Population$new("ResetDemo", vector_to_dataframe(rnorm(5)))
    #' pop$set_data(
    #'   data.frame(
    #'     subject_id = 1:8,
    #'     endpoint = rnorm(8),
    #'     population_name = "ResetDemo"
    #'   )
    #' )
    set_data = function(data) {
      if(!("arm" %in% names(data))){
        data$arm <- self$name
      }
      self$data <- data
      self$n <- nrow(self$data)
      self$dropped  <- rep(NA, self$n)
      self$enrolled <- rep(NA, self$n)
    }

  ) # end public
) # end class

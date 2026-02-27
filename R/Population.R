#' Population: Manage a patient population
#'
#' @description
#' The `Population` class stores subject-level data and manages enrollment and
#' dropout times for each subject.
#'
#' Use `set_enrolled()` and `set_dropped()` to assign enrollment and dropout
#' times to random subsets of subjects.
#'
#' @seealso [Trial] for multi-arm simulations, [Timer] for managing timepoints.
#'
#' @examples
#' # Basic example: vector input
#' pop <- Population$new(name = "Control", data = vector_to_dataframe(rnorm(10)))
#'
#' pop$n # number of subjects
#' head(pop$data) # generated subject-level data
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
    #' - `id` `integer`
    #' - `arm` `character`
    #' - `readout_time` `numeric`
    #' - `data` `numeric`
    #' - may contain more columns
    data = NULL,

    #' @field enrolled `numeric` Vector of enrollment times for each subject.
    enrolled = NULL,

    #' @field dropped `numeric` Vector of dropout times for each subject.
    dropped = NULL,

    #' @field n `integer` Number of unique subjects in the population.
    n = NULL,

    #' @field n_readouts `integer` Number of readout_times in the population.
    n_readouts = NULL,

    # --- constructor ---
    #' @description
    #' Create a new `Population` instance.
    #'
    #' @param name `character` Unique identifier for the population.
    #' @param data `data.frame` with columns: `id`, `arm`, `readout_time`,
    #'   `data`, and optionally more columns.
    #' @param enrolled `numeric` Optional enrollment times (auto-initialized if `NULL`).
    #' @param dropped `numeric` Optional dropout times (auto-initialized if `NULL`).
    #' @param n `integer` Auto-computed from data (optional).
    #' @param n_readouts `integer` Auto-computed from data (optional).
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
      n = NULL,
      n_readouts = NULL
    ) {
      stopifnot(is.character(name))
      self$name <- name
      if (!("arm" %in% names(data))) {
        data$arm <- name
      }

      # Check required data frame columns
      col_names <- c("id", "arm", "readout_time")
      missing_cols <- setdiff(col_names, names(data))
      if (length(missing_cols) > 0) {
        stop(sprintf("Data frame is missing required columns: %s", paste(missing_cols, sep = ", ")))
      }
      if (length(names(data)) < 4) {
        stop(sprintf("Data frame is missing endpoint data."))
      }

      self$data <- data
      self$n <- length(unique(self$data$id))
      self$n_readouts <- (dim(self$data)[1] / self$n)

      # Initialize enrollment/dropout status if not provided
      ifelse(
        is.null(enrolled),
        self$enrolled <- rep(NA, self$n),
        self$enrolled <- enrolled
      )
      ifelse(
        is.null(dropped),
        self$dropped <- rep(NA, self$n),
        self$dropped <- dropped
      )
    },

    # --- methods ---
    #' @description
    #' Mark subjects as enrolled at a given time.
    #'
    #' Enrollment applies only to unenrolled subjects (`NA`).
    #'
    #' @param n `integer` Number of subjects to enroll.
    #' @param time `numeric` Enrollment time.
    #'
    #' @seealso [Population], [Trial].
    #'
    #' @examples
    #' pop <- Population$new("Test", vector_to_dataframe(rnorm(10)))
    #' pop$set_enrolled(n = 4, time = 2)
    set_enrolled = function(n, time) {
      # input validation
      n <- as.integer(n)
      if (length(n) != 1L || is.na(n) || n < 0L) {
        stop("`n` must be a single non-negative integer.")
      }

      # Don't enroll more subjects than available
      idx <- which(is.na(self$enrolled))
      n_use <- min(n, length(idx))
      if (n_use == 0L) return(invisible(self))

      pick <- idx[sample.int(length(idx), n_use, replace = FALSE)]
      self$enrolled[pick] <- time
    },

    #' @description
    #' Mark subjects as dropped at a given time.
    #'
    #' Dropout applies only to enrolled, not-yet-dropped subjects.
    #'
    #' @param n `integer` Number of subjects to drop.
    #' @param time `numeric` Dropout time.
    #'
    #' @seealso [Population], [Trial].
    #'
    #' @examples
    #' pop <- Population$new("Test", vector_to_dataframe(rnorm(10)))
    #' pop$set_enrolled(n = 5, time = 1)
    #' pop$set_dropped(n = 2, time = 3)
    set_dropped = function(n, time) {
      # input validation
      n <- as.integer(n)
      if (length(n) != 1L || is.na(n) || n < 0L) {
        stop("`n` must be a single non-negative integer.")
      }

      # Don't drop more subjects than eligible (enrolled and not yet dropped)
      idx <- which(is.na(self$dropped) & !is.na(self$enrolled))
      n_use <- min(n, length(idx))
      if (n_use == 0L) return(invisible(self))

      pick <- idx[sample.int(length(idx), n_use, replace = FALSE)]
      self$dropped[pick] <- time
    },

    #' @description
    #' Replace underlying subject data and reset enrollment/dropout status.
    #'
    #' @param data `data.frame` with columns: `id`, `data`, `arm`,
    #'   `readout_time`, and optionally more columns.
    #'
    #' @seealso [Population].
    #'
    #' @examples
    #' pop <- Population$new("ResetDemo", vector_to_dataframe(rnorm(5)))
    #' pop$set_data(
    #'   data.frame(
    #'     id = 1:8,
    #'     data = rnorm(8),
    #'     arm = "ResetDemo",
    #'     readout_time = 0
    #'   )
    #' )
    set_data = function(data) {
      if (!("arm" %in% names(data))) {
        data$arm <- self$name
      }
      self$data <- data
      self$n <- nrow(self$data)
      self$dropped <- rep(NA, self$n)
      self$enrolled <- rep(NA, self$n)
    }
  ) # end public
) # end class

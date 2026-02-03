#' @export
Population <- R6::R6Class(
  classname = "Population",
  public = list(
    # fields
    name = NULL,
    data = NULL,
    dropped = NULL,
    enrolled = NULL,
    n=NULL,

    # constructor
    initialize = function(
    name,
    data = NULL,
    dropped = NULL,
    enrolled = NULL,
    n=NULL,
    ...
    ) {
      stopifnot(is.character(name))
      self$name <- name

      if (is.vector(data)) {
        self$data <- data.frame(
          subject_id = seq_along(data),
          data = data,
          population_name = self$name
        )
      } else {
        self$data <- data
      }

      self$n <- length(unique(self$data$subject_id))
      self$dropped <- rep(NA, self$n)
      self$enrolled <- rep(NA, self$n)
    },

    # methods
    set_dropped = function(n, time) {
      idx <- which(is.na(self$dropped) & !is.na(self$enrolled))
      self$dropped[sample(idx, n, replace = FALSE)] <- time
    },

    set_enrolled = function(n, time) {
      idx <- which(is.na(self$enrolled))
      self$enrolled[sample(idx, n, replace = FALSE)] <- time
    },

    set_data = function(data) {
      if (is.vector(data)) {
        self$data <- data.frame(
          subject_id = seq_along(data),
          data = data,
          population_name = self$name
        )
      } else {
        self$data <- data
      }
      self$n <- nrow(self$data)
      self$dropped  <- rep(NA, self$n)
      self$enrolled <- rep(NA, self$n)
    }

  ) # end public
) # end class

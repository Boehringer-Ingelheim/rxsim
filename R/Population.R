Population <- R6::R6Class(
  classname = "Population",
  public = list(
    # fields
    name = NULL,
    data = NULL,
    dropped = NULL,
    enrolled = NULL,

    # constructor
    initialize = function(
    name,
    data = NULL,
    dropped = NULL,
    enrolled = NULL,
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

      n <- nrow(self$data)
      self$dropped <- rep(NA, n)
      self$enrolled <- rep(NA, n)
    },

    # methods
    set_dropped = function(n, time) {
      potential=self$dropped[is.na(self$dropped)&!is.na(self$enrolled)]
      id <- sample(seq_len(length(potential)), n, replace = FALSE)
      self$dropped[is.na(self$dropped)&!is.na(self$enrolled)] [id]<- time
    },

    set_enrolled = function(n, time) {
      potential=self$enrolled[is.na(self$enrolled)]
      id <- sample(seq_len(length(potential)), n, replace = FALSE)
      self$enrolled[is.na(self$enrolled)] [id]<- time
    },

    get_data = function() {
      subset(self$data,!is.na(self$enrolled))
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
     n <- nrow(self$data)
      self$dropped <- rep(NA, n)
      self$enrolled <- rep(NA, n)
    }

  ) # end public
) # end class

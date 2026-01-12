Population <- R6::R6Class(
  "Population",

  public = list(
    name = NULL,
    data = NULL,
    dropped = NULL,
    enrolled = NULL,

    initialize = function(name, data = NULL, dropped = NULL, enrolled = NULL, ...) {
      stopifnot(is.character(name))
      self$name <- name

      if (is.vector(data)) {
        self$data <- cbind(seq_along(data), data)
      } else {
        self$data <- data
      }

      n <- nrow(self$data)
      self$dropped <- rep(NA, n)
      self$enrolled <- rep(NA, n)
    },



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
      self$data[!is.na(self$enrolled), ]
    },
    set_data = function(data) {
      if (is.vector(data)) {
        self$data <- cbind(seq_along(data), data)
      } else {
        self$data <- data
      }
     n <- nrow(self$data)
      self$dropped <- rep(NA, n)
      self$enrolled <- rep(NA, n)
    }
  )
)

ann <- Population$new(name="Ann", rnorm(50))
enroll<-cbind(1:10,rep(5,10))
drop<-cbind(1:10,rep(1,10))

for(i in 1:nrow(enroll)){
  ann$set_enrolled(time=enroll[i,1],n=enroll[i,2])
  ann$set_dropped(time=drop[i,1],n=drop[i,2])

  }


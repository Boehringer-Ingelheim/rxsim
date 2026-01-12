# --- Fix run_readers so it passes current_time into cond$func ---
Timers <- R6::R6Class(
  classname = "Timers",
  public = list(
    name = NULL,
    timelist = NULL,
    reader_conditions = NULL,

    initialize = function(name, timelist = NULL, reader_conditions = NULL) {
      stopifnot(is.character(name))
      self$name <- name
      self$timelist <- if (is.null(timelist)) list() else timelist
      self$reader_conditions <- if (is.null(reader_conditions)) list() else reader_conditions
    },

    add_timepoint = function(dropper, enroller) {
      stopifnot(is.numeric(dropper), is.numeric(enroller))
      tp <- list(dropper = dropper, enroller = enroller)
      self$timelist <- append(self$timelist, list(tp))
    },

    add_reader_condition = function(time = NULL, n_events = NULL, threshold = NULL, func = NULL) {
      cond <- list(time = time, n_events = n_events, threshold = threshold, func = func)
      self$reader_conditions <- append(self$reader_conditions, list(cond))
    },

    get_n_timepoints = function() length(self$timelist),

    get_timepoint = function(i) {
      if (i <= 0 || i > length(self$timelist)) stop("Index out of range")
      self$timelist[[i]]
    },

    run_readers = function(locked_data, current_time, n_events = NULL) {
      results <- list()
      for (i in seq_along(self$reader_conditions)) {
        cond <- self$reader_conditions[[i]]
        time_ok <- is.null(cond$time) || cond$time == current_time
        events_ok <- is.null(cond$n_events) || (!is.null(n_events) && n_events >= cond$n_events)
        threshold_ok <- is.null(cond$threshold) || (is.numeric(cond$threshold) && n_events >= cond$threshold)

        if (time_ok && events_ok && threshold_ok) {
          if (is.function(cond$func)) {
            # Pass both locked_data and current_time into the function
            results[[paste0("reader_", i)]] <- cond$func(locked_data, current_time)
          } else {
            results[[paste0("reader_", i)]] <- locked_data
          }
        }
      }
      results
    }
  )
)

# --- Usage ---
t <- Timers$new(name = "Experiment1")

t$add_timepoint(1, 5)
t$add_timepoint(2, 10)

# Reader condition with pre-check on current_time
t$add_reader_condition(
  time = 1,
  n_events = 5,
  func = function(data, time) {
    if (time < 6) {
      data <- 0
    }
    mean(data)
  }
)

t$add_reader_condition(
  time = 2,
  n_events = 10,
  func = function(data, time) sum(data)
)

locked_data <- c(1, 2, 3, 4, 5)
results <- t$run_readers(locked_data, current_time = 1, n_events = 10)

print(results)

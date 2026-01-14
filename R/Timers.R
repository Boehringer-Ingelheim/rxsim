
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

    # New: add reader with dplyr-style predicates
    # - '...' are filter-like boolean expressions, e.g., status == "active", visit >= 3
    # - 'func' will be called as func(filtered_data, current_time)
    # - 'name' becomes the result key ("reader_<name>")
    # Legacy params kept but ignored (warn) to ease migration
    add_reader_condition = function(..., func = NULL, name = NULL,
                                    time = NULL, n_events = NULL, threshold = NULL,
                                    .env = parent.frame()) {
      if (!requireNamespace("rlang", quietly = TRUE)) {
        stop("This version requires {rlang}.")
      }
      # Capture filter predicates as quosures (with caller env)
      where_quos <- rlang::enquos(..., .named = FALSE)


      cond <- list(
        where = where_quos,
        func  = func,
        name  = name
      )
      self$reader_conditions <- append(self$reader_conditions, list(cond))
      invisible(self)
    },

    get_n_timepoints = function() length(self$timelist),

    get_timepoint = function(i) {
      if (i <= 0 || i > length(self$timelist)) stop("Index out of range")
      self$timelist[[i]]
    },

    # New run_readers:
    # - Applies each reader's own filter predicates (cond$where) to locked_data
    # - Calls cond$func(filtered_data, current_time)
    # - Optional controls:
    #     .skip_empty: skip calling func if per-reader filtered data is empty
    #     .name_prefix: prefix for result keys (default "reader_")
    run_readers = function(locked_data, current_time,
                           .skip_empty = FALSE,
                           .name_prefix = "reader_") {

      if (!requireNamespace("dplyr", quietly = TRUE) ||
          !requireNamespace("rlang", quietly = TRUE)) {
        stop("run_readers requires {dplyr} and {rlang}.")
      }
      stopifnot(is.data.frame(locked_data))

      results <- list()

      for (i in seq_along(self$reader_conditions)) {
        cond <- self$reader_conditions[[i]]

        key <- if (!is.null(cond$name) && nzchar(cond$name)) {
          paste0(.name_prefix, cond$name)
        } else {
          paste0(.name_prefix, i)
        }

        # Per-reader filtering (dplyr semantics: NA in predicates drops rows)
        df_i <- if (!is.null(cond$where) && length(cond$where) > 0) {
          dplyr::filter(locked_data, !!!cond$where)
        } else {
          locked_data
        }

        if (.skip_empty && nrow(df_i) == 0L) {
          next
        }

        if (is.function(cond$func) && nrow(df_i) != 0L) {
          results[[key]] <- cond$func(df_i, current_time)
        } else {
          results[[key]] <- df_i
          warning(sprintf("Condition '%s' has no valid func; returning filtered data.", key), call. = FALSE)
        }
      }

      results
    }
  )
)

# Example reader funcs
summary_reader <- function(dat, t_now) {
  list(t = t_now, n = nrow(dat),
       mean_value = if ("value" %in% names(dat) && nrow(dat) > 0) mean(dat$value) else NA_real_)
}
ids_reader <- function(dat, t_now) unique(dat$subject_id)

# Sample data
set.seed(1)
df <- data.frame(
  subject_id = rep(1:5, each = 4),
  visit      = rep(1:4, times = 5),
  status     = sample(c("active", "inactive"), 20, replace = TRUE),
  value      = rnorm(20),
  center     = sample(c("EU","US"), 20, replace = TRUE)
)

timers <- Timers$new("TidyTimers")
# Each reader carries its own dplyr-like predicates
timers$add_reader_condition(status == "active", visit >= 3, func = summary_reader, name = "active_v3_sum")
timers$add_reader_condition(center == "EU", value > 0,       func = ids_reader,     name = "eu_ids")

out <- timers$run_readers(df, Sys.time())
str(out, 1)

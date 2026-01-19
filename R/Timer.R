Timer <- R6::R6Class(
  classname = "Timer",
  public = list(
    # fields
    name = NULL,
    timelist = NULL,
    conditions = NULL,

    # constructor
    initialize = function(
    name,
    timelist = NULL,
    conditions = NULL
    ) {
      stopifnot(is.character(name))
      self$name <- name
      self$timelist <- if (is.null(timelist)) list() else timelist
      self$conditions <- if (is.null(conditions)) list() else conditions
    },

    # methods
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
    add_condition = function(
    ...,
    analysis = NULL,
    name = NULL,
    time = NULL,
    n_events = NULL,
    threshold = NULL,
    .env = parent.frame()
    ) {
      # Capture filter predicates as quosures (with caller env)
      where_quos <- rlang::enquos(..., .named = FALSE)

      cond <- list(
        where = where_quos,
        analysis  = analysis,
        name  = name
      )
      self$conditions <- append(self$conditions, list(cond))
      invisible(self)
    },

    get_n_timepoints = function() length(self$timelist),

    get_timepoint = function(i) {
      if (i <= 0 || i > length(self$timelist)) stop("Index out of range")
      self$timelist[[i]]
    },

    # New check_conditions:
    # - Applies each reader's own filter predicates (cond$where) to locked_data
    # - Calls cond$func(filtered_data, current_time)
    # - Optional controls:
    #     .skip_empty: skip calling func if per-reader filtered data is empty
    #     .name_prefix: prefix for result keys (default "reader_")
    check_conditions = function(
    locked_data,
    current_time,
    .skip_empty = FALSE,
    .name_prefix = "reader_"
    ) {
      stopifnot(is.data.frame(locked_data))

      results <- list()

      for (i in seq_along(self$conditions)) {
        cond <- self$conditions[[i]]

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
          warning(sprintf("\n Condition '%s' has no valid func; returning filtered data.", key), call. = FALSE)
        }
      }

      results
    }

  ) # end public
) # end class

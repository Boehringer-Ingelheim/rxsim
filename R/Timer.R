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
    add_timepoint = function(time,arm,dropper, enroller) {
      stopifnot(is.numeric(dropper), is.numeric(enroller))
      tp <- list(time=time,arm=arm,dropper = dropper, enroller = enroller)
      self$timelist <- append(self$timelist, list(tp))
    },

    # New: add reader with dplyr-style predicates
    # - '...' are filter-like boolean expressions, e.g., status == "active", visit >= 3
    # - 'analysis' will be called as analysis(filtered_data, current_time)
    # - 'name' becomes the result key ("reader_<name>")
    # Legacy params kept but ignored (warn) to ease migration
    add_condition = function(
    ...,
    analysis = NULL,
    name = NULL
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

    get_end_timepoint = function() max(sapply(self$timelist,function(x){ x$time})),
    get_n_arms = function() length(unique(sapply(self$timelist,function(x) x$arm))),
    get_unique_times = function() unique(sapply(self$timelist,function(x) x$time)),

    get_timepoint = function(arm, i) {
      # Basic validation
      if (missing(arm)) stop("`arm` is required.")
      if (missing(i))   stop("`i` is required.")

      # Extract columns from list with type safety
      times <- vapply(self$timelist, function(x) x$time, FUN.VALUE = numeric(1))
      arms  <- vapply(self$timelist, function(x) x$arm,  FUN.VALUE = character(1))

      # Find matching indices
      idx <- which(times == i & arms == arm)

      # Handle match outcomes
      if (length(idx) == 0L) {
        return(NULL)
      }
      if (length(idx) > 1L) {
        stop(sprintf("Multiple timepoints found for arm = %s and time = %s.", arm, as.character(i)))
      }

      self$timelist[[idx]]
    },

    # New check_conditions:
    # - Applies each reader's own filter predicates (cond$where) to locked_data
    # - Calls cond$func(filtered_data, current_time)
    check_conditions = function(
    locked_data,
    current_time
    ) {
      stopifnot(is.data.frame(locked_data))

      results <- list()

      for (i in seq_along(self$conditions)) {
        cond <- self$conditions[[i]]

        key <- ifelse(
          !is.null(cond$name) && nzchar(cond$name),
          cond$name,
          i
        )

        # Per-reader filtering (dplyr semantics: NA in predicates drops rows)
        df_i <- if (!is.null(cond$where) && length(cond$where) > 0) {
          dplyr::filter(locked_data, !!!cond$where)
        } else {
          locked_data
        }

        if (nrow(df_i) == 0L) {
          warning(sprintf(" skipping check as filtered is empty \n"), call. = FALSE)
          next
        }

        if (is.function(cond$analysis)) {
          results[[key]] <- cond$analysis(df_i, current_time)
        } else {
          results[[key]] <- df_i
          warning(sprintf(" returning filtered data as is because condition '%s' has no applicable analysis \n", key), call. = FALSE)
        }
      }

      results
    }

  ) # end public
) # end class

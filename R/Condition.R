#' Condition: Stateful trigger and analysis unit
#'
#' @description
#' A `Condition` represents a single trigger rule that can be evaluated against
#' a data frame at a given time. It encapsulates:
#' \itemize{
#'   \item filtering logic using `dplyr::filter()` expressions
#'   \item optional analysis execution
#'   \item trigger bookkeeping (cooldown, trigger count, last trigger time)
#' }
#'
#' `Condition` objects are typically managed by a [`Timer`] instance, but can be
#' constructed and evaluated independently for testing or advanced use cases.
#'
#' @details
#' A condition triggers when:
#' \enumerate{
#'   \item The filtered data contains at least one row
#'   \item The cooldown period has elapsed since the last trigger
#'   \item The maximum number of triggers has not been exceeded
#' }
#'
#' On a successful trigger, the condition either:
#' \itemize{
#'   \item applies an analysis function to the filtered data, or
#'   \item returns the filtered data unchanged (with a warning)
#' }
#'
#' Trigger state is updated automatically after each successful trigger.
#'
#' @section Fields:
#' \describe{
#'   \item{\code{where}}{`list` of quosures capturing filtering expressions
#'     suitable for `dplyr::filter()`.}
#'   \item{\code{analysis}}{`function` or `NULL`. A function with signature
#'     `(data, current_time)` returning analysis results.}
#'   \item{\code{name}}{`character` or `NULL`. Optional unique identifier
#'     for the condition.}
#'   \item{\code{cooldown}}{`numeric`. Minimum time between consecutive triggers.}
#'   \item{\code{max_triggers}}{`integer` or `Inf`. Maximum number of allowed triggers.}
#'   \item{\code{trigger_count}}{`integer`. Number of times the condition has triggered.}
#'   \item{\code{last_trigger_time}}{`numeric`. Time of the most recent trigger.}
#' }
#'
#' @section Methods:
#' \describe{

#' }
#'
#' @seealso
#' \itemize{
#'   \item [`Timer`] for managing multiple conditions and timepoints
#'   \item [`dplyr::filter()`] for condition syntax
#' }
#'
#' @examples
#' library(dplyr)
#'
#' df <- data.frame(
#'   id = 1:6,
#'   arm = c("A","A","B","B","A","B"),
#'   status = c("active","inactive","active","active","inactive","active"),
#'   visit = c(1,2,1,3,3,2)
#' )
#'
#' analysis_fn <- function(dat, current_time) {
#'   out <- aggregate(id ~ arm, dat, length)
#'   out$current_time <- current_time
#'   out
#' }
#'
#' cond <- Condition$new(
#'   where = rlang::enquos(status == "active"),
#'   analysis = analysis_fn,
#'   name = "active_only",
#'   cooldown = 1,
#'   max_triggers = 2
#' )
#'



#' Condition: Stateful trigger and analysis unit
#'
#' @description
#' Represents a single trigger condition with optional analysis,
#' cooldown logic, and trigger count tracking.
#'
Condition <- R6::R6Class(
  classname = "Condition",

  public = list(
    # --- fields ---
    where = NULL,
    analysis = NULL,
    name = NULL,
    cooldown = 0,
    max_triggers = 1L,
    trigger_count = 0L,
    last_trigger_time = NA_real_,

    # --- constructor ---
    initialize = function(
    where,
    analysis = NULL,
    name = NULL,
    cooldown = 0,
    max_triggers = 1L
    ) {
      self$where <- where
      self$analysis <- analysis
      self$name <- name

      cooldown <- as.numeric(cooldown)
      if (length(cooldown) != 1L || cooldown < 0 || is.na(cooldown)) {
        stop("`cooldown` must be a single non-negative number.")
      }

      max_triggers <- as.integer(max_triggers)
      if (length(max_triggers) != 1L || max_triggers < 0 || is.na(max_triggers)) {
        stop("`max_triggers` must be a non-negative integer (use Inf for unlimited).")
      }

      self$cooldown <- cooldown
      self$max_triggers <- max_triggers
    },

    # --- methods ---

    check_conditions = function(
    locked_data,
    current_time
    ) {
      stopifnot(is.data.frame(locked_data))

      results <- list()



        key <- ifelse(
          !is.null(self$name) && nzchar(self$name),
          self$name,
          i
        )

        # Per-reader filtering (dplyr semantics: NA in predicates drops rows)
        df_i <- if (!is.null(self$where) && length(self$where) > 0) {
          dplyr::filter(locked_data, !!!self$where)
        } else {
          locked_data
        }


        if (is.null(self$trigger_count)) self$trigger_count <- 0L
        if (is.null(self$max_triggers))  self$max_triggers <- 1L
        if (is.null(self$last_trigger_time)) self$last_trigger_time <- NA_real_
        if (is.null(self$cooldown)) self$cooldown <- 0

        # match_now <- nrow(df_i) > 0L
        #
        #
        # # If no match, skip
        # if (!match_now) {
        #  # self$conditions[[i]] <- cond
        #   next
        # }
        #
        #
        # # Hard cap on number of triggers
        # if (is.finite(self$max_triggers) && self$trigger_count >= self$max_triggers) {
        #   #self$conditions[[i]] <- cond
        #   next
        # }
        #
        # # Check cooldown
        # if (is.finite(self$last_trigger_time)) {
        #   if ((current_time - self$last_trigger_time) < self$cooldown) {
        #     #self$conditions[[i]] <- cond
        #     next
        #   }
        # }


        if (is.function(self$analysis)) {
          results[[key]] <- self$analysis(df_i, current_time)
        } else {
          results[[key]] <- df_i
          warning(sprintf(" returning filtered data as is because condition '%s' has no applicable analysis \n", key), call. = FALSE)
        }

        # Update trigger info after a successful trigger
        self$trigger_count <- self$trigger_count + 1L
        self$last_trigger_time <- current_time

        # Persist state back into Timer
        #self$conditions[[i]] <- cond

      }

     # results

)
)



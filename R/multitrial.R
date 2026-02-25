
replicate_trial<- function(trial,n=1)
{
trials<-list()
for(i in 1:n){
temp<- Trial$new(
name = paste(trial$name, i),
timer = trial$timer$clone(),
population = lapply(trial$population,function(x){x$clone()}))
trials <- append(trials, temp)
}
return(trials)
}


#' Generate a Population Object from a Named Generator
#'
#' Helper to create a `Population` R6 object by running a provided generator function.
#'
#' @param name Character scalar; the population name (e.g., arm label).
#' @param generator A function with no required arguments that returns a data.frame-like object.
#'
#' @return A `Population` R6 object constructed as `Population$new(name = name, data = generator())`.
#'
#' @examples
#' gen_control <- function() data.frame(id = 1:5, value = rnorm(5))
#' pop <- gen_population(name = "control", generator = gen_control)
#'
#' @export
gen_population <-function(name=NULL,generator=NULL)
{
population<-Population$new(
name = name,
data = generator()
)
return(population)
}


create_multiple_trials_gen_pop<- function(trial_name='name',data_gen_list=NULL, timer=timer, n = 1)
{
trials <- list()
for (i in 1:n) {
# Create a cloned trial with updated name and cloned components
temp <- Trial$new(
name       = paste(trial_name, i),
timer      = timer$clone(),
population = lapply(data_gen_list, function(x) {gen_population(x[[1]],x[[2]]) })
)
trials <- append(trials, list(temp))
}
return(trials)
}
run_multiple_trial<-function(trials=NULL)
{
lapply(trials,function(x){ x$run()})
}

# ======================================================================
# Example usage (can be placed in a vignette or examples section)
# ======================================================================

#--- Configuration / Generators ---
sample_size <- 100
arms        <- c("pbo", "trt")
allocation  <- c(1, 1)
enrollment  <- list(end_time = c(4, 8, 12, 16), rate = c(5, 10, 15, 20))
dropout     <- list(end_time = c(4, 8, 12, 16), rate = c(1, 2, 4, 8))

timepoints <- gen_timepoints(
  sample_size = sample_size,
  arms        = arms,
  allocation  = allocation,
  enrollment  = enrollment,
  dropout     = dropout
)

t <- Timer$new(name = "trial_timer")
add_timepoints(t, timepoints)

t$add_condition(
  time %in% c(6, 12),
  analysis = function(df, time) {
    df_enrolled <- df |>
      dplyr::filter(!is.na(enroll_time))
    stats::t.test(value ~ arm, data = df_enrolled)$p.value
  },
  name = "t_test_final"
)

gen_control <- function() {
  data.frame(
    subject_id = 1:50,
    value = stats::rnorm(50)
  )
}

gen_trt <- function() {
  data.frame(
    subject_id = 1:50,
    value = stats::rnorm(50, mean = 0.2)
  )
}

pops_gen <- list(
  list("pbo", gen_control),
  list("trt", gen_trt)
)

# Create 100 trials and run them
trials <- create_multiple_trials_gen_pop(
  trial_name    = "trial",
  data_gen_list = pops_gen,
  timer         = t,
  n             = 100
)

results <- run_multiple_trial(trials)

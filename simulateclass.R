
#problem 3 grid search if you want to have different trial designs with grid search
#
# trial 1 which is one arm trial one look at final
#
#

scenario <- tidyr::expand_grid(
  sample_size = c(50, 100),
  mu = c(5, 10),
  enroll_end = list(c(4,8)),
  endroll_rate = list(c(10, 10))
  data_functions = list()
)


# trial 2 we have 2 arms 2 looks at interim at final
#

scenario <- tidyr::expand_grid(
  sample_size = c(50, 100),
  mu = c(5, 10),
  mu_pbo = c(6,12),
  enroll_end = list(c(4,8)),
  endroll_rate = list(c(10, 10))
  data_functions = list()
)

# trial 3 we have 6 arms
# 2 looks at interim at final
# interimno go at final bmcpmod

#The objective is one simulate function that can handle all of three




simulate<-function(name, timer, data_ga)
{
  populations<-NULL
  for(i in 1:length(data_gen))
  {
   populations.append(data_gen(n,parms))
  }
  trial <- Trial$new(
    name = name,
    timer = timer,
    population = populations
  )
  plan <- gen_timepoints(
    sample_size,
    c("A", "B"),
    c(1,1),
    data.frame(
      enroll
    ))
  # --- Run ---
  trial$run()

  prettify_results(trial$results)

}


scenario <- tidyr::expand_grid(
  sample_size = c(50, 100),
  mu = c(5, 10),
  mu_pbo = c(6,12),
  enroll_end = list(c(4,8)),
  endroll_rate = list(c(10, 10))
  data_functions = list()
)

simulate <- function(scenario, generator, n_reps){
  plan <- gen_timepoints(
    sample_size,
    c("A", "B"),
    c(1,1),
    data.frame(
      enroll
    ))

  t <- new

  pop <- Population$new(name = "pbo"
                        data = generator(n, mu))

  pop2 <- Population$new(name ="trt"
                         data =  generaotr(n, mu_pbo))

  trial <- new

  trial$run()



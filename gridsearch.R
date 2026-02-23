
#problem 2 grid search
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



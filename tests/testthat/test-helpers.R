# as_population_data()

test_that("as_population_data has the correct column names", {
  pop <- as_population_data(c(1.0, 2.0, 3.0))
  expect_equal(colnames(pop), c("id", "data", "readout_time"))
})

test_that("as_population_data sets id as sequential integers starting at 1", {
  pop <- as_population_data(c(10.0, 20.0, 30.0))
  expect_equal(pop$id, 1:3)
})

test_that("as_population_data preserves input data values", {
  vals <- c(1.5, -0.5, 3.0)
  pop <- as_population_data(vals)
  expect_equal(pop$data, vals)
})

test_that("as_population_data sets readout_time to 0 for all rows", {
  pop <- as_population_data(c(1.0, 2.0, 3.0))
  expect_equal(pop$readout_time, rep(0, 3))
})

test_that("as_population_data works with a single value", {
  pop <- as_population_data(42.0)
  df <- data.frame(id = 1L, data = 42.0, readout_time = 0)
  expect_equal(pop, df)
})

# add_timepoints()

test_that("add_timepoints errors when timer is not a Timer instance", {
  df <- data.frame(time = 1, arm = "A", enroll = 3L, drop = 0L)
  expect_error(add_timepoints("not_a_timer", df), "`timer` must be a Timer instance")
  expect_error(add_timepoints(list(), df), "`timer` must be a Timer instance")
})

test_that("add_timepoints errors when df is not a data.frame", {
  t <- Timer$new(name = "T")
  expect_error(add_timepoints(t, "not_a_df"), "`df` must be a data.frame")
  expect_error(add_timepoints(t, 1:5), "`df` must be a data.frame")
})

test_that("add_timepoints adds the correct number of timepoints", {
  t <- Timer$new(name = "T")
  df <- data.frame(
    time   = c(1, 2, 3),
    arm    = rep("A", 3),
    enroll = c(5L, 3L, 2L),
    drop   = c(0L, 1L, 1L)
  )
  add_timepoints(t, df)
  expect_equal(length(t$timelist), 3L)
})

test_that("add_timepoints returns the timer invisibly", {
  t <- Timer$new(name = "T")
  df <- data.frame(time = 1, arm = "A", enroll = 3L, drop = 0L)
  result <- add_timepoints(t, df)
  expect_true(inherits(result, "Timer"))
})

test_that("add_timepoints stores correct enroll and drop values", {
  t <- Timer$new(name = "T")
  df <- data.frame(time = 5, arm = "B", enroll = 4L, drop = 2L)
  add_timepoints(t, df)
  tp <- t$get_timepoint("B", 5)
  expect_equal(tp$enroll, 4L)
  expect_equal(tp$drop, 2L)
})

# get_col_names()

test_that("get_col_names always includes the four fixed time columns", {
  pop <- Population$new(name = "P", data = as_population_data(rnorm(5)))
  result <- get_col_names(pop)
  expect_true(all(c("time", "enroll_time", "drop_time", "measure_time") %in% result))
})

test_that("get_col_names returns columns from a single population", {
  pop <- Population$new(name = "P", data = as_population_data(rnorm(5)))
  result <- get_col_names(pop)
  expect_true(all(c("id", "data") %in% result))
})

test_that("get_col_names returns combined columns from a list of populations", {
  pop1 <- Population$new(
    name = "P1",
    data = data.frame(id = 1:5, score = rnorm(5), readout_time = 0)
  )
  pop2 <- Population$new(
    name = "P2",
    data = data.frame(id = 1:5, weight = rnorm(5), readout_time = 0)
  )
  result <- get_col_names(list(pop1, pop2))
  expect_true(all(c("score", "weight") %in% result))
})

test_that("get_col_names returns unique column names (no duplicates)", {
  pop1 <- Population$new(name = "P1", data = as_population_data(rnorm(5)))
  pop2 <- Population$new(name = "P2", data = as_population_data(rnorm(5)))
  result <- get_col_names(list(pop1, pop2))
  expect_equal(length(result), length(unique(result)))
})


# collect_results()

make_trial <- function(seed = 1L) {
  pop <- Population$new(name = "A", data = as_population_data(rnorm(10)))
  t   <- Timer$new(name = "T")
  t$add_timepoint(time = 1, arm = "A", enroll = 5L, drop = 0L)
  t$add_timepoint(time = 2, arm = "A", enroll = 3L, drop = 1L)

  cond <- Condition$new(
    where    = calendar_trigger(1),
    analysis = function(df, current_time) data.frame(n = nrow(df)),
    name     = "interim"
  )
  cond2 <- Condition$new(
    where    = calendar_trigger(2),
    analysis = function(df, current_time) data.frame(n = nrow(df)),
    name     = "final"
  )

  Trial$new(
    name       = "T1",
    seed       = seed,
    timer      = t,
    population = list(pop),
    conditions = list(cond, cond2)
  )
}

test_that("collect_results errors when input is not a Trial or list", {
  expect_error(collect_results("not_a_trial"), "`trials` must be a Trial object")
  expect_error(collect_results(42), "`trials` must be a Trial object")
})

test_that("collect_results errors when a list contains non-Trial elements", {
  trial <- make_trial()
  trial$run()
  expect_error(
    collect_results(list(trial, "oops")),
    "`trials` must be a Trial object or a list of Trial objects."
  )
})

test_that("collect_results returns an empty data.frame for a trial with no results", {
  trial <- make_trial()
  result <- collect_results(trial)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
})

test_that("collect_results result has the expected columns", {
  trial <- make_trial()
  trial$run()
  result <- collect_results(trial)
  expect_true(all(c("replicate", "timepoint", "analysis") %in% names(result)))
})

test_that("collect_results sets replicate = 1 for a single Trial", {
  trial <- make_trial()
  trial$run()
  result <- collect_results(trial)
  expect_true(all(result$replicate == 1L))
})

test_that("collect_results filters by analysis name", {
  trial <- make_trial()
  trial$run()
  result <- collect_results(trial, analysis = "interim")
  expect_true(all(result$analysis == "interim"))
})

test_that("collect_results returns all analyses when analysis = NULL", {
  trial <- make_trial()
  trial$run()
  result <- collect_results(trial)
  expect_true(all(result$analysis %in% c("interim", "final")))
})

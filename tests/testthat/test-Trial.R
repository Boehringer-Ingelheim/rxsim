test_that("Trial run: subject_id globally unique across arms with same n_readouts", {
  # Arm A: 3 subjects, 2 readouts each -> 6 rows
  dataA <- data.frame(
    id           = rep(1:3, each = 2),
    data         = c(1, 1, 2, 2, 3, 3),
    readout_time = rep(c(0, 1), 3)
  )
  popA <- Population$new("A", data = dataA)

  # Arm B: 3 subjects, 2 readouts each -> 6 rows
  dataB <- data.frame(
    id           = rep(1:3, each = 2),
    data         = c(4, 4, 5, 5, 6, 6),
    readout_time = rep(c(0, 1), 3)
  )
  popB <- Population$new("B", data = dataB)

  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 3L, dropper = 0L)
  timer$add_timepoint(time = 1, arm = "B", enroller = 3L, dropper = 0L)

  trigger_by_calendar(1, timer, analysis = function(df, current_time) df)

  trial <- Trial$new(
    name       = "two_arm_repeated",
    seed       = 1,
    timer      = timer,
    population = list(popA, popB)
  )

  trial$run()

  snap <- trial$locked_data[["time_1"]]

  # 6 subjects x 2 readouts = 12 rows
  testthat::expect_equal(nrow(snap), 12L)

  # 6 globally unique subject IDs
  testthat::expect_equal(sort(unique(snap$subject_id)), 1:6)

  # Each subject ID appears exactly 2 times (each, not cycled with times)
  testthat::expect_true(all(table(snap$subject_id) == 2L))
})

test_that("Trial initialize: errors when populations have different n_readouts", {
  popA <- Population$new("A", data = data.frame(
    id = 1:3, data = 1:3, readout_time = 0
  ))
  popB <- Population$new("B", data = data.frame(
    id = rep(1:3, each = 2), data = rep(1:3, each = 2), readout_time = rep(c(0, 1), 3)
  ))

  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 3L, dropper = 0L)
  timer$add_timepoint(time = 1, arm = "B", enroller = 3L, dropper = 0L)

  testthat::expect_error(
    Trial$new(name = "bad_trial", timer = timer, population = list(popA, popB)),
    "All populations must have the same n_readouts"
  )
})

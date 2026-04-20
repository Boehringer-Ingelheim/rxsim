test_that("Trial run: subject_id globally unique with mixed n_readouts across arms", {
  # Arm A: 3 subjects, 1 readout each -> 3 rows
  dataA <- data.frame(
    id           = c(1, 2, 3),
    value        = c(1, 2, 3),
    readout_time = 0
  )
  popA <- Population$new("A", data = dataA, n_readouts = 1L)

  # Arm B: 3 subjects, 2 readouts each -> 6 rows
  dataB <- data.frame(
    id           = c(1, 1, 2, 2, 3, 3),
    value        = c(4, 4, 5, 5, 6, 6),
    readout_time = rep(c(0, 1), 3)
  )
  popB <- Population$new("B", data = dataB, n_readouts = 2L)

  timer <- Timer$new("t")
  timer$add_timepoint(time = 1, arm = "A", enroller = 3L, dropper = 0L)
  timer$add_timepoint(time = 1, arm = "B", enroller = 3L, dropper = 0L)

  trigger_by_calendar(1, timer, analysis = function(df, current_time) df)

  trial <- Trial$new(
    name       = "mixed_readouts",
    seed       = 1,
    timer      = timer,
    population = list(popA, popB)
  )

  trial$run()

  snap <- trial$locked_data[["time_1"]]

  expect_equal(nrow(snap), 9L)

  expect_equal(sort(unique(snap$subject_id)), 1:6)

  rows_A <- snap[snap$subject_id %in% 1:3, ]
  expect_true(all(table(rows_A$subject_id) == 1))

  rows_B <- snap[snap$subject_id %in% 4:6, ]
  expect_true(all(table(rows_B$subject_id) == 2))
})

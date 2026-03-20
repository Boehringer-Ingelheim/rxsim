testthat::test_that("Population initialize: basic type checks and arm check", {

  set.seed(1)
  df <- data.frame(
    id = 1:5,
    readout_time = 1,
    value = rnorm(5)
  )

  pop <- Population$new("Arm A", data = df)

  # Field types
  testthat::expect_true(is.character(pop$name))
  testthat::expect_true(is.data.frame(pop$data))
  testthat::expect_true(is.numeric(pop$enrolled))
  testthat::expect_true(is.numeric(pop$dropped))
  testthat::expect_true(is.integer(pop$n))
  testthat::expect_true(is.integer(pop$n_readouts))

  # Values
  testthat::expect_equal(pop$name, "Arm A")
  testthat::expect_true("arm" %in% names(pop$data))
  testthat::expect_true(all(pop$data$arm == "Arm A"))

  testthat::expect_equal(pop$n, length(unique(pop$data$id)))
  testthat::expect_equal(length(pop$enrolled), pop$n)
  testthat::expect_equal(length(pop$dropped), pop$n)
  testthat::expect_true(all(is.na(pop$enrolled)))
  testthat::expect_true(all(is.na(pop$dropped)))
})

testthat::test_that("Population initialize: errors on wrong name type", {


  df <- data.frame(
    id = 1:3,
    readout_time = 1,
    value = rnorm(3),
    arm = "X"
  )

  testthat::expect_error(Population$new(123, data = df))
  testthat::expect_error(Population$new(NULL, data = df))
})

testthat::test_that("Population initialize: errors when required columns are missing", {


  # Missing readout_time
  df_bad <- data.frame(
    id = 1:5,
    value = rnorm(5),
    arm = "Arm A"
  )
  testthat::expect_error(
    Population$new("Arm A", data = df_bad),
    "Data frame is missing required columns",
    fixed = FALSE
  )

  # Missing endpoint data
  df_no_endpoint <- data.frame(
    id = 1:5,
    arm = "Arm A",
    readout_time = 1
  )
  testthat::expect_error(
    Population$new("Arm A", data = df_no_endpoint),
    "Data frame is missing endpoint data",
    fixed = FALSE
  )
})

testthat::test_that("Population initialize: computes n and n_readouts correctly (repeated measures)", {


  set.seed(1)
  df <- rbind(
    data.frame(id = 1:10, readout_time = 1, value = rnorm(10)),
    data.frame(id = 1:10, readout_time = 4, value = rnorm(10))
  )

  pop <- Population$new("Arm A", data = df)

  testthat::expect_equal(pop$n, 10)
  testthat::expect_equal(pop$n_readouts, 2)

  testthat::expect_equal(length(pop$enrolled), 10)
  testthat::expect_equal(length(pop$dropped), 10)
})

testthat::test_that("Population initialize: inputs correctly provided enrolled/dropped vectors", {


  df <- data.frame(
    id = 1:5,
    readout_time = 1,
    value = rnorm(5)
  )

  enrolled <- c(1, NA, 4, NA, NA)
  dropped  <- c(NA, NA, NA, 3, NA)

  pop <- Population$new("Arm A", data = df, enrolled = enrolled, dropped = dropped)

  testthat::expect_true(is.numeric(pop$enrolled)) #Replace in first test
  testthat::expect_true(is.numeric(pop$dropped))

  testthat::expect_equal(pop$enrolled, enrolled)
  testthat::expect_equal(pop$dropped, dropped)
})

# testthat::test_that("set_enrolled: validates n type and preserves numeric time values", {
#
#
#   set.seed(10)
#   df <- data.frame(
#     id = 1:10,
#     readout_time = 1,
#     value = rnorm(10)
#   )
#   pop <- Population$new("Arm A", data = df)
#
#   # bad n
#   testthat::expect_error(pop$set_enrolled(n = -1, time = 1), "`n` must be a single non-negative integer")
#   testthat::expect_error(pop$set_enrolled(n = NA, time = 1), "`n` must be a single non-negative integer")
#   testthat::expect_error(pop$set_enrolled(n = c(1, 2), time = 1), "`n` must be a single non-negative integer")
#
#   # time confirm time stored as provided
#   pop$set_enrolled(n = 3, time = 4.5)
#   testthat::expect_true(any(pop$enrolled == 4.5, na.rm = TRUE))
#   testthat::expect_true(is.numeric(pop$enrolled))
# })

testthat::test_that("set_enrolled: Don't enroll more subjects than available and dont replace already enrolled ", {


  set.seed(1)
  df <- data.frame(
    id = 1:10,
    readout_time = 1,
    value = rnorm(10)
  )
  pop <- Population$new("Arm A", data = df)

  pop$set_enrolled(n = 4, time = 1)
  enrolled_before <- pop$enrolled

  pop$set_enrolled(n = 999, time = 2)

  # should cap to n
  testthat::expect_equal(sum(!is.na(pop$enrolled)), pop$n)

  # original enrolled at time 1 remain 1
  testthat::expect_true(all(pop$enrolled[!is.na(enrolled_before)] == 1))
})

testthat::test_that("set_enrolled: milestone counts by time match requests", {


  set.seed(1)
  df <- data.frame(
    id = 1:12,
    readout_time = 1,
    value = rnorm(12)
  )
  pop <- Population$new("Arm A", data = df)

  pop$set_enrolled(n = 2, time = 1)
  pop$set_enrolled(n = 3, time = 2)
  pop$set_enrolled(n = 4, time = 2)

  tt <- table(pop$enrolled, useNA = "ifany")

  testthat::expect_equal(as.integer(tt[["1"]]), 2)
  testthat::expect_equal(as.integer(tt[["2"]]), 7)
  testthat::expect_equal(sum(!is.na(pop$enrolled)), 9)
  testthat::expect_equal(sum(is.na(pop$enrolled)), pop$n - 9)
})

testthat::test_that("set_enrolled / set_dropped: n = 0 does nothing", {


  set.seed(1)
  df <- data.frame(
    id = 1:8,
    readout_time = 1,
    value = rnorm(8)
  )
  pop <- Population$new("Arm A", data = df)

  enrolled_before <- pop$enrolled
  dropped_before  <- pop$dropped

  pop$set_enrolled(n = 0, time = 1)
  pop$set_dropped(n = 0, time = 1)

  testthat::expect_equal(pop$enrolled, enrolled_before)
  testthat::expect_equal(pop$dropped, dropped_before)
})

testthat::test_that("set_enrolled and set_dropped validate n and apply eligibility rules", {

  set.seed(1)

  df <- data.frame(
    id = 1:10,
    readout_time = 1,
    value = rnorm(10)
  )
  pop <- Population$new("Arm A", data = df)


  # set_enrolled: input validation

  testthat::expect_error(
    pop$set_enrolled(n = -1, time = 1),
    "`n` must be a single non-negative integer",
    fixed = TRUE
  )
  testthat::expect_error(
    pop$set_enrolled(n = NA, time = 1),
    "`n` must be a single non-negative integer",
    fixed = TRUE
  )
  testthat::expect_error(
    pop$set_enrolled(n = c(1, 2), time = 1),
    "`n` must be a single non-negative integer",
    fixed = TRUE
  )

  # set_enrolled: time stored as provided + type
  pop$set_enrolled(n = 3, time = 4.5)
  testthat::expect_true(any(pop$enrolled == 4.5, na.rm = TRUE))
  testthat::expect_true(is.numeric(pop$enrolled))

  # set_dropped: eligibility + input validation

  dropped_before <- pop$dropped
  pop$set_dropped(n = 3, time = 2)

  # Enroll more subjects, then drop
  pop$set_enrolled(n = 5, time = 1)
  pop$set_dropped(n = 2, time = 3)

  testthat::expect_equal(sum(!is.na(pop$dropped)), 5)
  testthat::expect_true(all(pop$dropped[!is.na(pop$dropped)] %in% c(2,3)))

  dropped_idx <- which(!is.na(pop$dropped))
  testthat::expect_true(all(!is.na(pop$enrolled[dropped_idx])))

  testthat::expect_true(is.numeric(pop$dropped))
})

testthat::test_that("set_dropped: milestone counts by time match requests and subset-of-enrolled holds", {


  set.seed(1)
  df <- data.frame(
    id = 1:12,
    readout_time = 1,
    value = rnorm(12)
  )
  pop <- Population$new("Arm A", data = df)

  pop$set_enrolled(n = 10, time = 1)
  pop$set_dropped(n = 3, time = 3)
  pop$set_dropped(n = 2, time = 4)

  dt <- table(pop$dropped, useNA = "ifany")
  testthat::expect_equal(as.integer(dt[["3"]]), 3)
  testthat::expect_equal(as.integer(dt[["4"]]), 2)
  testthat::expect_equal(sum(!is.na(pop$dropped)), 5)

  dropped_idx <- which(!is.na(pop$dropped))
  testthat::expect_true(all(!is.na(pop$enrolled[dropped_idx])))
  testthat::expect_true(sum(!is.na(pop$dropped)) <= sum(!is.na(pop$enrolled)))
})

testthat::test_that("set_dropped: caps to eligible slots and does not re-drop", {


  set.seed(1)
  df <- data.frame(
    id = 1:6,
    readout_time = 1,
    value = rnorm(6)
  )
  pop <- Population$new("Arm A", data = df)

  pop$set_enrolled(n = 3, time = 1)
  pop$set_dropped(n = 999, time = 2)

  testthat::expect_equal(sum(!is.na(pop$dropped)), 3)

  dropped_before <- sum(!is.na(pop$dropped))
  pop$set_dropped(n = 2, time = 3)
  testthat::expect_equal(sum(!is.na(pop$dropped)), dropped_before)
})

testthat::test_that("set_dropped: never assigns dropout to non-enrolled slots", {


  set.seed(1)
  df <- data.frame(
    id = 1:12,
    readout_time = 1,
    value = rnorm(12)
  )
  pop <- Population$new("Arm A", data = df)

  pop$set_enrolled(n = 5, time = 1)
  pop$set_dropped(n = 5, time = 2)

  testthat::expect_true(all(is.na(pop$dropped) | !is.na(pop$enrolled)))
})

testthat::test_that("set_data: resets enrolled/dropped and updates n", {


  set.seed(14)
  df <- rbind(
    data.frame(id = 1:10, readout_time = 1, value = rnorm(10)),
    data.frame(id = 1:10, readout_time = 4, value = rnorm(10))
  )
  pop <- Population$new("Arm A", data = df)

  pop$set_enrolled(n = 5, time = 1)
  pop$set_dropped(n = 2, time = 2)

  df2 <- data.frame(
    id = 1:8,
    readout_time = 1,
    value = rnorm(8)
  )
  pop$set_data(df2)


  testthat::expect_equal(pop$n, nrow(df2))
  testthat::expect_equal(length(pop$enrolled), pop$n)
  testthat::expect_equal(length(pop$dropped), pop$n)
  testthat::expect_true(all(is.na(pop$enrolled)))
  testthat::expect_true(all(is.na(pop$dropped)))
})

testthat::test_that("Set_data: Test that repeated measures parameters are inputted and calculated correctly", {


  df <- data.frame(
    id = 1:10,
    readout_time = 1,
    value = rnorm(10)
  )
  pop <- Population$new("Arm A", data = df)

  df2 <- rbind(
    data.frame(id = 1:10, readout_time = 1, value = rnorm(10)),
    data.frame(id = 1:10, readout_time = 4, value = rnorm(10))
  )
  pop$set_data(df2)

  testthat::expect_equal(pop$n, length(unique(df2$id)))
  testthat::expect_equal(pop$n_readouts, nrow(df2) / pop$n)
  testthat::expect_equal(length(pop$enrolled), pop$n)
  testthat::expect_equal(length(pop$dropped), pop$n)

})

# Statistical ground-truth tests for the stochastic assignment paths.
# Each test compares realised output against an analytic baseline (or a
# cross-path reference), not a snapshot. Seeds are fixed, so the tests are
# deterministic; tolerances are sized several sigma wide so they validate the
# distribution without being flaky.

# Pull per-subject follow-up time (drop_time - enroll_time) from a full reveal.
.followups <- function(adaptive, reps, ss, enroll_rate, drop_rate) {
  do.call(c, lapply(seq_len(reps), function(s) {
    set.seed(s)
    sch <- stochastic_schedule(ss, "A", 1,
             enrollment = function(n) rexp(n, enroll_rate),
             dropout    = function(n) rexp(n, drop_rate))
    timer <- Timer$new("t"); add_timepoints(timer, sch)
    pop <- Population$new("A", as_population_data(rnorm(ss)))
    reveal <- condition_calendar_time(max(timer$get_unique_times()),
                                      analysis = function(df, t) df)
    tr <- Trial$new("e", seed = s, timer = timer, population = list(pop),
                    conditions = list(reveal), adaptive = adaptive)
    suppressWarnings(tr$run())
    if (length(tr$locked_data) == 0L) return(numeric(0))
    sn <- tr$locked_data[[length(tr$locked_data)]]
    sn <- sn[!duplicated(sn$subject_id), c("enroll_time", "drop_time")]
    sn <- sn[!is.na(sn$drop_time), , drop = FALSE]
    sn$drop_time - sn$enroll_time
  }))
}

test_that("enroll inter-event times are Exp(rate): mean = 1/rate", {
  # stochastic_schedule places enrollments at cumsum(enrollment(n)); the gaps
  # between sorted enroll times recover the iid draws, so E[gap] = 1/rate.
  set.seed(101)
  rate <- 0.5
  gaps <- unlist(lapply(1:40, function(i)
    diff(sort(stochastic_schedule(150L, "A", 1,
              enrollment = function(n) rexp(n, rate), dropout = NULL)$time))))
  testthat::expect_equal(mean(gaps), 1 / rate, tolerance = 0.1)  # ~5+ sigma
})

test_that("drop arm allocation follows multinomial(ratio)", {
  # Drop events get an arm via sample(arms, prob = ratio); pooled proportions
  # must match the allocation ratio.
  set.seed(102)
  arms <- c("A", "B", "C"); ratio <- c(3, 2, 1) / 6
  cnt <- c(0, 0, 0)
  for (i in 1:40) {
    sch <- stochastic_schedule(150L, arms, c(3, 2, 1),
             enrollment = function(n) rexp(n, 1), dropout = function(n) rexp(n, 1))
    d <- sch[sch$drop == 1L, ]
    cnt <- cnt + as.integer(table(factor(d$arm, levels = arms)))
  }
  testthat::expect_equal(as.numeric(cnt / sum(cnt)), ratio, tolerance = 0.03)
})

test_that("enroll arm counts equal the exact allocation target", {
  # shuffled_arms is a permutation of rep(arms, target), so per-arm enrollment
  # totals are exact (not just approximate).
  set.seed(103)
  arms <- c("A", "B", "C")
  sch <- stochastic_schedule(120L, arms, c(3, 2, 1),
           enrollment = function(n) rexp(n, 1), dropout = NULL)
  observed <- vapply(arms, function(a) sum(sch$enroll[sch$arm == a]), integer(1))
  testthat::expect_equal(unname(observed), c(60L, 40L, 20L))
})

test_that("fixed and adaptive paths produce the same follow-up distribution", {
  # Both paths assign drops at random among eligible subjects; the per-subject
  # follow-up (drop - enroll) distributions must be statistically identical.
  f <- .followups(adaptive = FALSE, reps = 25, ss = 80L, enroll_rate = 0.05, drop_rate = 0.03)
  a <- .followups(adaptive = TRUE,  reps = 25, ss = 80L, enroll_rate = 0.05, drop_rate = 0.03)
  ks <- suppressWarnings(stats::ks.test(f, a))
  testthat::expect_gt(ks$p.value, 0.01)                       # not distinguishable
  testthat::expect_equal(median(f), median(a), tolerance = 0.1)
})

test_that("fixed-path drop assignment is uncorrelated when eligibility doesn't bind", {
  # enroll ~ Exp(1) finishes long before drop ~ Exp(0.01) begins, so every
  # subject is drop-eligible and random assignment makes enroll_time and
  # drop_time independent: cor ~ 0, far below the earliest-enrolled artefact
  # ceiling (n-1)/(n+5) ~ 0.94.
  ss <- 100L
  cors <- unlist(lapply(1:60, function(s) {
    set.seed(s)
    sch <- stochastic_schedule(ss, "A", 1,
             enrollment = function(n) rexp(n, 1), dropout = function(n) rexp(n, 0.01))
    timer <- Timer$new("t"); add_timepoints(timer, sch)
    pop <- Population$new("A", as_population_data(rnorm(ss)))
    reveal <- condition_calendar_time(max(timer$get_unique_times()),
                                      analysis = function(df, t) df)
    tr <- Trial$new("e", seed = s, timer = timer, population = list(pop),
                    conditions = list(reveal))
    suppressWarnings(tr$run())
    sn <- tr$locked_data[[length(tr$locked_data)]]
    sn <- sn[!duplicated(sn$subject_id), ]
    sn <- sn[!is.na(sn$drop_time), ]
    if (nrow(sn) < 3L) return(NA_real_)
    cor(sn$enroll_time, sn$drop_time)
  }))
  ceiling_cor <- (ss - 1) / (ss + 5)
  testthat::expect_lt(abs(mean(cors, na.rm = TRUE)), 0.2)
  testthat::expect_lt(mean(cors, na.rm = TRUE), 0.5 * ceiling_cor)
})

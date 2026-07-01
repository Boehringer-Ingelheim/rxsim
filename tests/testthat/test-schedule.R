testthat::test_that("deterministic_schedule: handles small sample sizes across many arms", {
  set.seed(1)
  sch <- deterministic_schedule(
    sample_size = 2,
    arms = c("A", "B", "C"),
    allocation = c(1, 1, 1),
    enrollment = list(end_time = 1, rate = 3),
    dropout = list(end_time = 1, rate = 0)
  )

  testthat::expect_true(is.data.frame(sch))
  testthat::expect_equal(sum(sch$enroll), 2L)
  testthat::expect_true(all(sch$enroll >= 0L))
  testthat::expect_true(all(sch$drop >= 0L))
})

testthat::test_that("deterministic_schedule: keeps per-arm enrollment at target", {
  arms <- c("A", "B", "C")
  allocation <- c(3, 2, 1)
  sample_size <- 12

  sch <- deterministic_schedule(
    sample_size = sample_size,
    arms = arms,
    allocation = allocation,
    enrollment = list(end_time = c(2, 4), rate = c(4, 3)),
    dropout = list(end_time = c(2, 4), rate = c(0, 1))
  )

  expected <- as.integer(allocation / sum(allocation) * sample_size)
  names(expected) <- arms
  observed <- vapply(arms, function(a) sum(sch$enroll[sch$arm == a]), integer(1))
  testthat::expect_equal(observed, expected)
  testthat::expect_equal(sum(observed), sample_size)
})

testthat::test_that("stochastic_schedule: enrollment totals match sample_size", {
  set.seed(99)
  sch <- stochastic_schedule(
    sample_size = 25,
    arms = c("A", "B"),
    allocation = c(1, 2),
    enrollment = function(n) rexp(n, 1),
    dropout = function(n) rexp(n, 0.5)
  )

  testthat::expect_equal(sum(sch$enroll), 25L)
  testthat::expect_equal(sum(sch$drop), 25L)
})

### Statistical ground-truth tests for stochastic assignment ###

test_that("enroll inter-event times are Exp(rate): mean = 1/rate", {
  # stochastic_schedule places enrollments at cumsum(rexp(n, rate)); cumsum is
  # monotone so diff(sort(time)) recovers the iid Exp(rate) gaps. E[gap] = 1/rate,
  # sd(gap) = 1/rate. 40 schedules x 149 gaps = 5960 iid draws.
  # CLT: mean ~ N(2, (1/rate)/sqrt(5960) = 0.02591). tol 0.1 => 3.86 sigma,
  # P(a random-seed run exceeds tol) = 2*pnorm(-3.86) = 1.13e-4. Seed fixes it.
  set.seed(101)
  rate <- 0.5
  gaps <- unlist(lapply(1:40, function(i)
    diff(sort(stochastic_schedule(150L, "A", 1,
              enrollment = function(n) rexp(n, rate), dropout = NULL)$time))))
  obs_mean <- mean(gaps)
  cat(sprintf("enroll-gap mean: obs %.4f target %.4f (n=%d)\n", obs_mean, 1 / rate, length(gaps)))
  testthat::expect_equal(obs_mean, 1 / rate, tolerance = 0.1)
})

test_that("drop arm allocation follows multinomial(ratio)", {
  # Drop arm ~ iid Categorical(ratio) via sample(arms, prob = ratio). N = 40*150 =
  # 6000 draws. Binomial-proportion sd = sqrt(p(1-p)/N); tightest arm A (p=0.5):
  # sd 0.00645, tol 0.03 => 4.65 sigma, P(fail) ~ 3.4e-6. Seed fixes it.
  set.seed(102)
  arms <- c("A", "B", "C"); ratio <- c(3, 2, 1) / 6
  cnt <- c(0, 0, 0)
  for (i in 1:40) {
    sch <- stochastic_schedule(150L, arms, c(3, 2, 1),
             enrollment = function(n) rexp(n, 1), dropout = function(n) rexp(n, 1))
    d <- sch[sch$drop == 1L, ]
    cnt <- cnt + as.integer(table(factor(d$arm, levels = arms)))
  }
  obs <- as.numeric(cnt / sum(cnt))
  cat(sprintf("drop-arm props: obs [%.3f %.3f %.3f] target [%.3f %.3f %.3f]\n",
              obs[1], obs[2], obs[3], ratio[1], ratio[2], ratio[3]))
  testthat::expect_equal(obs, as.numeric(ratio), tolerance = 0.03)
})

test_that("enroll arm counts equal the exact allocation target", {
  # shuffled_arms is a permutation of rep(arms, target), so per-arm enrollment
  # totals are EXACT (deterministic; P(fail) = 0, seed irrelevant to counts).
  set.seed(103)
  arms <- c("A", "B", "C"); allocation <- c(3, 2, 1); sample_size <- 120L
  target <- as.integer(sample_size * allocation / sum(allocation))  # 60, 40, 20
  sch <- stochastic_schedule(sample_size, arms, allocation,
           enrollment = function(n) rexp(n, 1), dropout = NULL)
  observed <- vapply(arms, function(a) sum(sch$enroll[sch$arm == a]), integer(1))
  cat(sprintf("enroll arm counts: obs [%d %d %d] target [%d %d %d]\n",
              observed[1], observed[2], observed[3], target[1], target[2], target[3]))
  testthat::expect_equal(unname(observed), target)
})

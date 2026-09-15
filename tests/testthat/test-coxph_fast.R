test_that("coxph_fast returns named vector of length 5", {
  set.seed(1)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_length(res, 5L)
  expect_named(res, c("coef", "exp(coef)", "se(coef)", "lower .95", "upper .95"))
})

test_that("coxph_fast exp(coef) is positive", {
  set.seed(2)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_gt(unname(res["exp(coef)"]), 0)
})

test_that("coxph_fast CI lower <= exp(coef) <= upper", {
  set.seed(3)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_lte(unname(res["lower .95"]), unname(res["exp(coef)"]) + 1e-10)
  expect_gte(unname(res["upper .95"]), unname(res["exp(coef)"]) - 1e-10)
})

test_that("coxph_fast coef = log(exp(coef))", {
  set.seed(4)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_equal(unname(res["coef"]), log(unname(res["exp(coef)"])),
               tolerance = 1e-12)
})

test_that("coxph_fast presorted=TRUE and presorted=FALSE agree", {
  set.seed(5)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)
  ord   <- order(time)

  res_pre <- coxph_fast(time[ord], event[ord], group[ord],
                        control = 1, presorted = TRUE)
  res_uns <- coxph_fast(time, event, group,
                        control = 1, presorted = FALSE)
  expect_equal(unname(res_pre), unname(res_uns), tolerance = 1e-12)
})

test_that("coxph_fast conf.level argument changes CI label names", {
  set.seed(6)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res90 <- coxph_fast(time, event, group, control = 1, conf.level = 0.90)
  expect_named(res90, c("coef", "exp(coef)", "se(coef)", "lower .90", "upper .90"))
})

test_that("coxph_fast handles factor group correctly", {
  set.seed(7)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- factor(rep(c("control", "treatment"), each = n / 2))

  res <- coxph_fast(time, event, group, control = "control")
  expect_gt(unname(res["exp(coef)"]), 0)
  expect_false(any(is.na(res)))
})

test_that("coxph_fast agrees with coxph (no ties, tolerance 1e-4)", {
  skip_if_not_installed("survival")
  set.seed(8)
  n     <- 200
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res <- coxph_fast(time, event, group, control = 1)

  fit <- survival::coxph(
    survival::Surv(time, event) ~ I(group != 1),
    ties    = "breslow",
    control = survival::coxph.control(
      eps        = 1e-12,
      toler.chol = .Machine$double.eps^0.875
    )
  )
  hr_cox <- unname(exp(stats::coef(fit)[1L]))

  expect_equal(unname(res["exp(coef)"]), hr_cox, tolerance = 1e-4)
})

test_that("coxph_fast returns all-NA when no events", {
  time  <- c(1, 2, 3, 4)
  event <- c(0, 0, 0, 0)
  group <- c(1, 1, 2, 2)

  res <- coxph_fast(time, event, group, control = 1)
  expect_true(all(is.na(res)))
})

test_that("coxph_fast stops when input lengths differ", {
  expect_error(
    coxph_fast(c(1, 2, 3), c(1, 0), c(1, 1, 2), control = 1),
    "same length"
  )
})

test_that("coxph_fast stratified agrees with coxph + strata() (tolerance 1e-4)", {
  skip_if_not_installed("survival")
  set.seed(4821)
  n <- 300
  strat  <- sample(1:3, n, replace = TRUE)
  group  <- sample(c("A", "B"), n, replace = TRUE)
  strata_haz <- c(1, 1.8, 0.6)[strat]
  time_true  <- rexp(n, rate = strata_haz * exp(0.5 * (group == "B")))
  cens       <- rexp(n, rate = 0.3)
  event      <- as.integer(time_true <= cens)
  time       <- pmin(time_true, cens)

  res <- coxph_fast(time, event, group, control = "A", strata = strat)

  fit <- survival::coxph(
    survival::Surv(time, event) ~ I(group != "A") + survival::strata(strat),
    ties    = "breslow",
    control = survival::coxph.control(
      eps        = 1e-12,
      toler.chol = .Machine$double.eps^0.875
    )
  )
  hr_cox <- unname(exp(stats::coef(fit)[1L]))

  expect_equal(unname(res["exp(coef)"]), hr_cox, tolerance = 1e-4)
})

test_that("coxph_fast single stratum reduces to unstratified", {
  set.seed(9)
  n     <- 100
  time  <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event <- rep(1L, n)
  group <- rep(c(1L, 2L), each = n / 2)

  res_strat <- coxph_fast(time, event, group, control = 1,
                          strata = rep(1L, n))
  res_uns   <- coxph_fast(time, event, group, control = 1)

  expect_equal(unname(res_strat), unname(res_uns), tolerance = 1e-12)
})

test_that("coxph_fast stratified presorted=TRUE and presorted=FALSE agree", {
  set.seed(10)
  n      <- 150
  strata <- sample(1:3, n, replace = TRUE)
  time   <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.15))
  event  <- rep(1L, n)
  group  <- rep(c(1L, 2L), each = n / 2)

  res_uns <- coxph_fast(time, event, group, control = 1, strata = strata)

  ord <- order(strata, time)
  res_pre <- coxph_fast(time[ord], event[ord], group[ord], control = 1,
                        strata = strata[ord], presorted = TRUE)

  expect_equal(unname(res_pre), unname(res_uns), tolerance = 1e-12)
})

test_that("coxph_fast handles a stratum with only one group", {
  # Stratum 2 contains only group 1 subjects; should not error and should
  # contribute zero events/expectation from that stratum, matching the
  # convention used by the stratified log-rank implementation.
  set.seed(11)
  n1     <- 80
  n2     <- 20
  time   <- c(rexp(n1, 0.1), rexp(n1, 0.15), rexp(n2, 0.2))
  event  <- rep(1L, n1 * 2 + n2)
  group  <- c(rep(1L, n1), rep(2L, n1), rep(1L, n2))
  strata <- c(rep(1L, n1 * 2), rep(2L, n2))

  res <- coxph_fast(time, event, group, control = 1, strata = strata)
  expect_false(any(is.na(res)))
})

test_that("coxph_fast stops when strata length differs from time", {
  expect_error(
    coxph_fast(c(1, 2, 3), c(1, 0, 1), c(1, 1, 2), control = 1,
              strata = c(1, 1)),
    "same length"
  )
})

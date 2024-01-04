test_that("LRT works on edge cases", {
  x <- c(10, 0, 0, 0, 0)

  expect_true(lrt_men_g4(x = x, g1 = 0, g2 = 0, pp = TRUE, dr = FALSE)$p_value != 0)
  expect_true(lrt_men_g4(x = x, g1 = 0, g2 = 0, pp = FALSE, dr = FALSE)$p_value != 0)
})


test_that("qq plot is unif in some cases", {
  g1 <- 2
  g2 <- 1
  alpha <- 1/6
  xi1 <- 1/3
  xi2 <- 1/3
  pp <- FALSE
  dr <- FALSE
  n <- 1000
  iter <- 1000
  pvec <- rep(NA_real_, iter)
  stat <- rep(NA_real_, iter)
  df <- rep(NA_real_, iter)
  gf <- matrix(ncol = 5, nrow = iter)
  for (i in seq_len(iter)) {
    x <- simf1g(n = n, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
    lout <- lrt_men_g4(x = x, g1 = g1, g2 = g2, pp = pp, dr = dr, alpha = alpha, xi1 = xi1, xi2 = xi2)
    pvec[[i]] <- lout$p_value
    stat[[i]] <- lout$statistic
    df[[i]] <- lout$df
    gf[i, ] <- offspring_gf_2(alpha = lout$alpha, xi1 = lout$xi1, xi2 = lout$xi2, p1 = g1, p2 = g2)
  }

  table(df)

  hwep::qqpvalue(pvals = pvec)
  qqplot(x = ppoints(iter), y = pvec)
  abline(a = 0, b = 1, col = 2, lty = 2)

  pvec2 <- stats::pchisq(q = stat, df = 2, lower.tail = FALSE)
  qqplot(x = ppoints(iter), y = pvec2)
  abline(a = 0, b = 1, col = 2, lty = 2)
})

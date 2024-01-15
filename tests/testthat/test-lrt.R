test_that("LRT works on edge cases", {
  x <- c(10, 0, 0, 0, 0)

  expect_true(lrt_men_g4(x = x, g1 = 0, g2 = 0, pp = TRUE, dr = FALSE)$p_value != 0)
  expect_true(lrt_men_g4(x = x, g1 = 0, g2 = 0, pp = FALSE, dr = FALSE)$p_value != 0)
})


test_that("qq plot is unif in some cases", {
  skip("too slow")


  g1 <- 1
  g2 <- 0
  alpha <- 1/12
  xi1 <- 1/3
  xi2 <- 1/3
  pp <- TRUE
  dr <- TRUE
  n <- 1000
  iter <- 100
  pvec <- rep(NA_real_, iter)
  stat <- rep(NA_real_, iter)
  df <- rep(NA_real_, iter)
  aest <- rep(NA_real_, iter)
  xi1est <- rep(NA_real_, iter)
  xi2est <- rep(NA_real_, iter)
  gf <- matrix(ncol = 5, nrow = iter)
  for (i in seq_len(iter)) {
    if (i %% 10 == 0) {
      cat(i, " of ", iter, "\n")
    }
    x <- simf1g(n = n, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
    lout <- lrt_men_g4(x = x, g1 = g1, g2 = g2, pp = pp, dr = dr, alpha = alpha, xi1 = xi1, xi2 = xi2)
    pvec[[i]] <- lout$p_value
    stat[[i]] <- lout$statistic
    df[[i]] <- lout$df
    #gf[i, ] <- offspring_gf_2(alpha = lout$alpha, xi1 = lout$xi1, xi2 = lout$xi2, p1 = g1, p2 = g2)
    aest[[i]] <- lout$alpha
    xi1est[[i]] <- lout$xi1
    xi2est[[i]] <- lout$xi2
  }

  table(df)

  ## hwep::qqpvalue(pvals = pvec)
  qqplot(x = ppoints(iter), y = pvec, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)

  ## competitors
  pvec2 <- stats::pchisq(q = stat, df = 2, lower.tail = FALSE)
  qqplot(x = ppoints(iter), y = pvec2, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)

  pvec3 <- stats::pchisq(q = stat, df = ifelse(aest < alpha - 1e-7, 4, 3), lower.tail = FALSE)
  qqplot(x = ppoints(iter), y = pvec3, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)
})


test_that("GL qq plot is unif in some cases", {
  skip("too slow")


  g1 <- 4
  g2 <- 2
  alpha <- 0
  xi1 <- 1/3
  xi2 <- 1/3
  pp <- TRUE
  dr <- TRUE
  rd <- 10
  n <- 100
  iter <- 100
  pvec <- rep(NA_real_, iter)
  stat <- rep(NA_real_, iter)
  df <- rep(NA_real_, iter)
  aest <- rep(NA_real_, iter)
  xi1est <- rep(NA_real_, iter)
  xi2est <- rep(NA_real_, iter)
  gf <- matrix(ncol = 5, nrow = iter)
  for (i in seq_len(iter)) {
    cat(i, " of ", iter, "\n")
    gl <- simf1gl(n = n, g1 = g1, g2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2, rd = rd)
    lout <- lrt_men_gl4(gl = gl, g1 = g1, g2 = g2, pp = pp, dr = dr, alpha = alpha, xi1 = xi1, xi2 = xi2)
    pvec[[i]] <- lout$p_value
    stat[[i]] <- lout$statistic
    df[[i]] <- lout$df
    #gf[i, ] <- offspring_gf_2(alpha = lout$alpha, xi1 = lout$xi1, xi2 = lout$xi2, p1 = g1, p2 = g2)
    aest[[i]] <- lout$alpha
    xi1est[[i]] <- lout$xi1
    xi2est[[i]] <- lout$xi2
  }

  table(df)

  ## hwep::qqpvalue(pvals = pvec)
  qqplot(x = ppoints(iter), y = pvec, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)

  ## competitors
  pvec2 <- stats::pchisq(q = stat, df = 3, lower.tail = FALSE)
  qqplot(x = ppoints(iter), y = pvec2, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)

  pvec3 <- stats::pchisq(q = stat, df = ifelse(aest > 1/6 - 1e-7, 2, 1), lower.tail = FALSE)
  qqplot(x = ppoints(iter), y = pvec3, xlim = c(0, 1), ylim = c(0, 1))
  abline(a = 0, b = 1, col = 2, lty = 2)
})

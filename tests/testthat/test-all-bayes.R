test_that("glknown works", {
  skip("takes too long")
  pdf <- expand.grid(p1 = 0:4, p2 = 0:4)
  for (i in seq_len(nrow(pdf))) {
    x <- round(100 * offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = pdf$p1[[i]], p2 = pdf$p2[[i]]))
    trash <- capture.output(
    suppressWarnings(
        {
        t1 <- bayes_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = TRUE, chains = 1, iter = 100)
        t2 <- bayes_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = TRUE, chains = 1, iter = 100)
        t3 <- bayes_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = FALSE, chains = 1, iter = 100)
        t4 <- bayes_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = FALSE)
        }
      )
    )
  }
})

test_that("glpknown works", {
  skip("takes too long")
  pdf <- expand.grid(p1 = 0:4, p2 = 0:4)
  for (i in seq_len(nrow(pdf))) {
    set.seed(i)
    gl <- simf1gl(n = 10, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]])
    trash <- capture.output(
    suppressWarnings(
        {
        t1 <- bayes_men_gl4(gl = gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = TRUE, chains = 1, iter = 100)
        t2 <- bayes_men_gl4(gl = gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = TRUE, chains = 1, iter = 100)
        t3 <- bayes_men_gl4(gl = gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = FALSE, chains = 1, iter = 100)
        t4 <- bayes_men_gl4(gl = gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = FALSE, chains = 1, iter = 100)
        }
      )
    )
  }
})

test_that("lrt_men_g4 works", {
  expect_no_error({
    pdf <- expand.grid(p1 = 0:4, p2 = 0:4)
    for (i in seq_len(nrow(pdf))) {
      set.seed(i)
      x <- round(100 * offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = pdf$p1[[i]], p2 = pdf$p2[[i]]))
      t1 <- lrt_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = TRUE)
      t2 <- lrt_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = FALSE)
      t3 <- lrt_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = TRUE)
      t4 <- lrt_men_g4(x = x, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = FALSE)
    }
  })
})

test_that("lrt_men_gl4 works", {
  skip("takes too long")
  expect_no_error({
    pdf <- expand.grid(p1 = 0:4, p2 = 0:4)
    for (i in seq_len(nrow(pdf))) {
      set.seed(i)
      gl <- simf1gl(n = 10, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]])
      t1 <- lrt_men_gl4(gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = TRUE)
      t2 <- lrt_men_gl4(gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = TRUE, dr = FALSE)
      t3 <- lrt_men_gl4(gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = TRUE)
      t4 <- lrt_men_gl4(gl, g1 = pdf$p1[[i]], g2 = pdf$p2[[i]], pp = FALSE, dr = FALSE)
      t5 <- lrt_men_gl4(gl, pp = TRUE, dr = TRUE)
      t6 <- lrt_men_gl4(gl, pp = TRUE, dr = FALSE)
      t7 <- lrt_men_gl4(gl, pp = FALSE, dr = TRUE)
      t8 <- lrt_men_gl4(gl, pp = FALSE, dr = FALSE)
    }
  })
})

test_that("hard data alt", {
  skip("takes too long")
  load(file = "./hard_alt.RData")
  lout <- lrt_men_g4(x = x, g1 = g1, g2 = g2)
  bout <- bayes_men_g4(x = x, g1 = g1, g2 = g2)
  lout$p_value
  bout$lbf

  lout <- lrt_men_gl4(gl = gl, g1 = g1, g2 = g2)
  bout <- bayes_men_gl4(gl = gl, g1 = g1, g2 = g2, chains = 1)
  lout$p_value
  bout$lbf

  g1 <- 0
  g2 <- 0
  x <- simf1g(n = 20, g1 = g1, g2 = g2, alpha = 1/6, xi1 = 1/3, xi2 = 1/3)
  lout <- lrt_men_g4(x = x, g1 = g1, g2 = g2)
  bout <- bayes_men_g4(x = x, g1 = g1, g2 = g2, chains = 1)
  lout$p_value
  bout$lbf

  gl <- simf1gl(n = 100, rd = 10, g1 = g1, g2 = g2, alpha = 0, xi1 = 1/3, xi2 = 1/3)
  lout <- lrt_men_gl4(gl = gl, g1 = g1, g2 = g2)
  bout <- bayes_men_gl4(gl = gl, g1 = g1, g2 = g2, chains = 1)
  lout$p_value
  bout$lbf
})

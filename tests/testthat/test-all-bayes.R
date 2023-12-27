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

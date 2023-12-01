#####################
## Likelihood methods when genotypes are known
#####################

#' Likelihood under three parameter model when genotypes are known
#'
#' @param x A vector of length 5. \code{x[i]} is the count of individuals with
#'     genotype \code{i-1}.
#' @param alpha The probability of quadrivalent formation.
#' @param beta The probability of double reduction given quadrivalent formation.
#' @param gamma1 The probability of AA_aa pairing for parent 1.
#' @param gamma2 The probability of AA_aa pairing for parent 2.
#' @param g1 Parent 1's genotype.
#' @param g2 Parent 2's genotype.
#' @param log_p A logical. Should we return the log likelihood or not?
#'
#' @return The (log) likelihood.
#'
#' @examples
#' x <- c(1, 4, 5, 3, 1)
#' tau <- 0.5
#' beta <- 0.1
#' gamma1 <- 0.5
#' gamma2 <- 0.3
#' g1 <- 1
#' g2 <- 2
#' like_gknown(
#'   x = x,
#'   tau = tau,
#'   beta = beta,
#'   gamma1 = gamma1,
#'   gamma2 = gamma2,
#'   g1 = g1,
#'   g2 = g2)
#'
#' @author David Gerard
#'
#' @export
like_gknown <- function(x, tau, beta, gamma1, gamma2, g1, g2, log_p = TRUE) {
  stopifnot(length(x) == 5)
  stopifnot(tau >= 0, tau <= 1,
            beta >= 0, beta <= 1,
            gamma1 >= 0, gamma1 <= 1,
            x >= 0,
            g1 >= 0, g1 <= 4,
            g2 >= 0, g2 <= 4)

  gf <- offspring_gf_3(
    tau = tau,
    beta = beta,
    gamma1 = gamma1,
    gamma2 = gamma2,
    p1 = g1,
    p2 = g2)

  return(stats::dmultinom(x = x, prob = gf, log = log_p))
}

#' Likelihood ratio test assuming no double reduction and no preferential pairing
#'
#' @inheritParams like_gknown
#'
#' @return A list of length three with the following elements
#' \describe{
#'   \item{\code{statistic}}{The likelihood ratio test statistic.}
#'   \item{\code{df}}{The degrees of freedom.}
#'   \item{\code{p_value}}{The p-value.}
#' }
#'
#' @author David Gerard
lrt_ndr_npp_g4 <- function(x, g1, g2) {

  l1 <- stats::dmultinom(x = x, prob = x / sum(x), log = TRUE)
  l0 <- like_gknown(
    x = x,
    tau = 1,
    beta = 0,
    gamma1 = 1/3,
    gamma2 = 1/3,
    g1 = g1,
    g2 = g2,
    log_p = TRUE)

  llr <- -2 * (l0 - l1)
  df <- 4
  p <- stats::pchisq(q = llr, df = df, lower.tail = FALSE, log.p = FALSE)

  return(list(statistic = llr, df = df, p_value = p))
}

#' Calculate degrees of freedom of test
#'
#' @param alpha estimated double reduction rate
#' @param xi1 estimated preferential pairing parameter of parent 1
#' @param xi2 estimated preferential pairing parameter of parent 2
#' @param g1 genotype of parent 1.
#' @param g2 genotype of parent 2
#' @param drbound the maximum possible value of alpha
#'
#' @examples
#' find_df(alpha = 1/6, xi1 = 1/3, xi2 = 1/3, g1 = 2, g2 = 2)
#' find_df(alpha = 1/12, xi1 = 5/33, xi2 = 23/33, g1 = 2, g2 = 2)
#' find_df(alpha = 1/12, xi1 = 5/33, xi2 = 1/3, g1 = 2, g2 = 2)
#' find_df(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, g1 = 2, g2 = 2)
#' find_df(alpha = 0, xi1 = 1/3, xi2 = 1/3, g1 = 2, g2 = 2)
#' find_df(alpha = 0, xi1 = 1/3, xi2 = 1/3, g1 = 2, g2 = 1)
#' find_df(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, g1 = 2, g2 = 1)
#' find_df(alpha = 1/12, xi1 = 1/3, xi2 = 1/3, g1 = 1, g2 = 1)
#' find_df(alpha = 0, xi1 = 1/3, xi2 = 1/3, g1 = 1, g2 = 1)
#'
#'
#' @author David Gerard
#'
#' @export
find_df <- function(alpha, xi1, xi2, g1, g2, drbound = 1/6) {
  TOL <- sqrt(.Machine$double.eps)

  lxi <- (1/3) * (alpha / (1 - alpha)) * ((1 - drbound) / drbound)
  uxi <- 1 - (2/3) * (alpha / (1 - alpha)) * ((1 - drbound) / drbound)

  ## consider alpha on boundary
  if (alpha >= drbound - TOL) {
    return(c(df = 0))
  } else if (alpha > TOL) {
    df <- 1
  } else {
    df <- 0
  }

  ## Consider xi1 on boundary
  if (g1 == 2) {
    df <- df + 1
    if (abs(xi1 - lxi) <= TOL) {
      df <- df - 1
    } else if (abs(uxi - xi1) <= TOL) {
      df <- df - 1
    }
  }

  ## Consider xi2 on boundary
  if (g2 == 2) {
    df <- df + 1
    if (abs(xi2 - lxi) <= TOL) {
      df <- df - 1
    } else if (abs(uxi - xi2) <= TOL) {
      df <- df - 1
    }
  }

  return(c(df = df))
}



lrt_dr_pp_g4 <- function(x, g1, g2, drbound = 1/6) {

  ## max likelihood under alt
  l1 <- stats::dmultinom(x = x, prob = x / sum(x), log = TRUE)

  ## Find MLE under null

  ## Get df and calculate p-value


}





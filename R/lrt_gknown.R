#####################
## Likelihood methods when genotypes are known
#####################

#' Likelihood under three parameter model when genotypes are known
#'
#' @param x A vector of length 5. \code{x[i]} is the count of individuals with
#'     genotype \code{i-1}.
#' @param tau The probability of quadrivalent formation.
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
#' This is when all genotypes are known.
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
  p_value <- stats::pchisq(q = llr, df = df, lower.tail = FALSE, log.p = FALSE)

  return(list(statistic = llr, p_value = p_value, df = df))
}

## Methods when dr and pp are not known ---------------------------------------

#' Objective when estimating both double reduction and preferential pairing
#'
#' @param par either just alpha (when no parent genotype is 2)
#'   or (tau, beta, gamma1) or (tau, beta, gamma2) or
#'   (tau, beta, gamma1, gamma2).
#' @param x offspring genotype counts
#' @param g1 first parent genotype
#' @param g2 second parent genotype
#'
#' @author David Gerard
#'
#' @noRd
obj_dr_pp <- function(par, x, g1, g2) {
  if (g1 != 2 && g2 != 2) {
    stopifnot(length(par) == 1)
    obj <- like_gknown(
      x = x,
      tau = 1,
      beta = par[[1]],
      gamma1 = 1/3,
      gamma2 = 1/3,
      g1 = g1,
      g2 = g2,
      log_p = TRUE)
  } else if (g1 == 2 && g2 != 2){
    stopifnot(length(par) == 3)
    obj <- like_gknown(
      x = x,
      tau = par[[1]],
      beta = par[[2]],
      gamma1 = par[[3]],
      gamma2 = 1/3,
      g1 = g1,
      g2 = g2,
      log_p = TRUE)
  } else if (g1 != 2 && g2 == 2){
    stopifnot(length(par) == 3)
    obj <- like_gknown(
      x = x,
      tau = par[[1]],
      beta = par[[2]],
      gamma1 = 1/3,
      gamma2 = par[[3]],
      g1 = g1,
      g2 = g2,
      log_p = TRUE)
  } else {
    stopifnot(length(par) == 4)
    obj <- like_gknown(
      x = x,
      tau = par[[1]],
      beta = par[[2]],
      gamma1 = par[[3]],
      gamma2 = par[[4]],
      g1 = g1,
      g2 = g2,
      log_p = TRUE)
  }
  return(obj)
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
#' @noRd
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

#' Initial parameter values
#'
#' @param g1 parent 1 genotype
#' @param g2 parent 2 genotype.
#' @param drbound The maximum double reduction rate
#' @param type How should we initialize? "random" or "half"?
#'
#' @return A list of length 3. The first element is the initial values,
#'     the second is the lower bound of parameters, the third element
#'     is the upper bound of parameters.
#'
#' @author David Gerard
#'
#' @noRd
lrt_init <- function(g1, g2, drbound = 1/6, type = c("random", "half")) {
  fudge <- 10^-6
  type <- match.arg(type)
  if (type == "random") {
    mult <- stats::runif(4)
  } else if (type == "half") {
    mult <- c(0.5, 0.5, 1/3, 1/3)
  }

  if (g1 != 2 && g2 != 2) {
    out <- list(
      par = c(alpha = drbound * mult[[1]]),
      lower = fudge,
      upper = drbound)
  } else if (g1 == 2 & g2 != 2) {
    out <- list(
      par = c(tau = mult[[1]], beta = drbound * mult[[2]], gamma1 = mult[[3]]),
      lower = rep(fudge, length.out = 3),
      upper = c(1 - fudge, drbound, 1 - fudge))
  } else if (g1 != 2 & g2 == 2) {
    out <- list(
      par = c(tau = mult[[1]], beta = drbound * mult[[2]], gamma2 = mult[[3]]),
      lower = rep(fudge, length.out = 3),
      upper = c(1 - fudge, drbound, 1 - fudge))
  } else {
    out <- list(
      par = c(tau = mult[[1]], beta = drbound * mult[[2]], gamma1 = mult[[3]], gamma2 = mult[[4]]),
      lower = rep(fudge, length.out = 4),
      upper = c(1 - fudge, drbound, 1 - fudge, 1 - fudge))
  }
  return(out)
}


#' LRT when both double reduction and preferential pairing are not known.
#'
#' All genotypes are known.
#'
#' @section Impossible genotypes:
#' Some offspring genotype combinations are impossible given the parental
#' gentoypes. If these impossible genotypes combinations show up, we return a
#' p-value of 0, a log-likelihood ratio statistic of Infinity, and missing
#' values for all other return items. The impossible genotypes are:
#' \describe{
#'   \item{\code{g1 = 0 && g2 = 0}}{Only offspring genotypes of 0 are possible.}
#'   \item{\code{g1 = 4 && g2 = 4}}{Only offspring genotypes of 4 are possible.}
#'   \item{\code{g1 = 0 && g2 = 4 || g1 == 4 && g2 == 0}}{Only offspring genotypes of 2 are possible.}
#'   \item{\code{g1 = 0 && g2 %in% c(1, 2, 3) || g1 = %in% c(1, 2, 3) && g2 == 0}}{Only offspring genotypes of 0, 1, and 2 are possible.}
#'   \item{\code{g1 = 4 && g2 %in% c(1, 2, 3) || g1 = %in% c(1, 2, 3) && g2 == 4}}{Only offspring genotypes of 2, 3, and 4 are possible.}
#' }
#'
#'
#' @inheritParams like_gknown
#' @param drbound the upper bound on the double reduction rate.
#' @param ntry The number of times to run the gradient ascent.
#'
#' @author David Gerard
#'
#' @examples
#' alpha <- 1/12
#' xi1 <- 1/3
#' xi2 <- 1/3
#' n <- 1000
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 0, p2 = 0)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 0, g2 = 0)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 0)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 0)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 1)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 1)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 2)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 3)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 3)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 1, p2 = 4)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 1, g2 = 4)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 0)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 0)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 1)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 1)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 2)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 3)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 3)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 2, p2 = 4)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 2, g2 = 4)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 0)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 0)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 1)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 1)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 2)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 2)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 3)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 3)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 3, p2 = 4)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 3, g2 = 4)
#'
#' gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = 4, p2 = 4)
#' x <- offspring_geno(gf = gf, n = n)
#' lrt_dr_pp_g4(x = x, g1 = 4, g2 = 4)
#'
#' @export
lrt_dr_pp_g4 <- function(x, g1, g2, drbound = 1/6, ntry = 5) {

  ## Deal with impossible values ----
  if (g1 == 0 && g2 == 0) {
    ret <- list(
      statistic = ifelse(all(x[2:5] == 0), Inf, 0),
      p_value = ifelse(all(x[2:5] == 0), 1, 0),
      df = NA_real_,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_)
    return(ret)
  } else if (g1 == 0 && g2 == 4 | g1 == 4 && g2 == 0) {
    ret <- list(
      statistic = ifelse(all(x[c(1, 2, 4, 5)] == 0), Inf, 0),
      p_value = ifelse(all(x[c(1, 2, 4, 5)] == 0), 1, 0),
      df = NA_real_,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_)
    return(ret)
  } else if (g1 == 4 && g2 == 4) {
    ret <- list(
      statistic = ifelse(all(x[1:4] == 0), Inf, 0),
      p_value = ifelse(all(x[1:4] == 0), 1, 0),
      df = NA_real_,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_)
    return(ret)
  } else if (g1 %in% c(1, 2, 3) && g2 == 0 | g1 == 0 && g2 %in% c(1, 2, 3)) {
    if (!all(x[4:5] == 0)) {
      ret <- list(
        statistic = Inf,
        p_value = 0,
        df = NA_real_,
        alpha = NA_real_,
        xi1 = NA_real_,
        xi2 = NA_real_)
      return(ret)
    }
  } else if (g1 %in% c(1, 2, 3) && g2 == 4 | g1 == 4 && g2 %in% c(1, 2, 3)) {
    if (!all(x[1:2] == 0)) {
      ret <- list(
        statistic = Inf,
        p_value = 0,
        df = NA_real_,
        alpha = NA_real_,
        xi1 = NA_real_,
        xi2 = NA_real_)
      return(ret)
    }
  }

  ## max likelihood under alt
  l1 <- stats::dmultinom(x = x, prob = x / sum(x), log = TRUE)

  ## Find MLE under null
  l0 <- -Inf
  for (i in seq_len(ntry)) {
    params <- lrt_init(g1 = g1, g2 = g2, drbound = drbound, type = "random")
    oout <- stats::optim(
      par = params$par,
      fn = obj_dr_pp,
      method = "L-BFGS-B",
      lower = params$lower,
      upper = params$upper,
      control = list(fnscale = -1),
      x = x,
      g1 = g1,
      g2 = g2)

    if (oout$value > l0) {
      l0 <- oout$value
    }
  }

  ## Get df and calculate p-value
  if (g1 != 2 & g2 != 2) {
    alpha <- oout$par[[1]]
    xi1 <- NA_real_
    xi2 <- NA_real_
  } else if (g1 == 2 & g2 != 2) {
    tau <- oout$par[[1]]
    beta <- oout$par[[2]]
    gamma1 <- oout$par[[3]]
    two <- three_to_two(tau = tau, beta = beta, gamma = gamma1)
    alpha <- two[[1]]
    xi1 <- two[[2]]
    xi2 <- NA_real_
  } else if (g1 != 2 & g2 == 2) {
    tau <- oout$par[[1]]
    beta <- oout$par[[2]]
    gamma2 <- oout$par[[3]]
    two <- three_to_two(tau = tau, beta = beta, gamma = gamma2)
    alpha <- two[[1]]
    xi1 <- NA_real_
    xi2 <- two[[2]]
  } else if (g1 == 2 & g2 == 2) {
    tau <- oout$par[[1]]
    beta <- oout$par[[2]]
    gamma1 <- oout$par[[3]]
    gamma2 <- oout$par[[4]]
    two <- three_to_two(tau = tau, beta = beta, gamma = gamma1)
    alpha <- two[[1]]
    xi1 <- two[[2]]
    two <- three_to_two(tau = tau, beta = beta, gamma = gamma2)
    stopifnot(alpha == two[[1]])
    xi2 <- two[[2]]
  }
  df_null <- find_df(alpha = alpha, xi1 = xi1, xi2 = xi2, g1 = g1, g2 = g2, drbound = drbound)
  df <- 4 - df_null

  llr <- -2 * (l0 - l1)

  p_value <- stats::pchisq(q = llr, df = df, lower.tail = FALSE)

  ret <-  list(
    statistic = llr,
    p_value = p_value,
    df = df,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2)

  return(ret)
}





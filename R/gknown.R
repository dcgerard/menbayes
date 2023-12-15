###################
## Methods when genotypes are known
###################

#' Checks to see if x is impossible given parental gentoypes
#'
#' @param x The offspring genotype counts
#' @param g1 Parent 1's genotype
#' @param g2 Parent 2's genotype
#' @param df Is double reduction possible?
#'
#' @return TRUE if impossible, FALSE if possible
#'
#' @author David Gerard
#'
#' @examples
#' is_impossible(x = c(1, 2, 3, 0, 0), g1 = 0, g2 = 3, dr = TRUE)
#'
#'
#' @noRd
is_impossible <- function(x, g1, g2, dr = TRUE) {
  if (g1 == 0 && g2 == 0) {
    return(!all(x[2:5] == 0))
  } else if (g1 == 4 & g2 == 4) {
    return(!all(x[1:4] == 0))
  } else if (g1 == 0 && g2 == 4 || g1 == 4 && g2 == 0) {
    return(!all(x[c(1, 2, 4, 5)] == 0))
  } else if (dr) {
    if (g1 == 0 && g2 %in% c(1, 2, 3) || g1 %in% c(1, 2, 3) && g2 == 0) {
      return(!all(x[4:5] == 0))
    } else if (g1 == 4 && g2 %in% c(1, 2, 3) || g1 %in% c(1, 2, 3) && g2 == 4) {
      return(!all(x[1:2] == 0))
    }
  } else if (!dr) {
    if (g1 == 0 && g2 == 1 || g1 == 1 && g2 == 0) {
      return(!all(x[3:5] == 0))
    } else if (g1 == 0 && g2 == 2 || g1 == 2 && g2 == 0 || g1 == 1 && g2 == 1) {
      return(!all(x[4:5] == 0))
    } else if (g1 == 0 && g2 == 3 || g1 == 3 && g2 == 0) {
      return(!all(x[c(1, 4, 5)] == 0))
    } else if (g1 == 1 && g2 == 4 || g1 == 4 && g2 == 1) {
      return(!all(x[c(1, 2, 5)] == 0))
    } else if (g1 == 2 && g2 == 4 || g1 == 4 && g2 == 2 || g1 == 3 && g2 == 3) {
      return(!all(x[1:2] == 0))
    } else if (g1 == 3 && g2 == 4 || g1 == 4 && g2 == 3) {
      return(!all(x[1:3] == 0))
    } else if (g1 == 1 && g2 == 2 || g1 == 2 && g2 == 1) {
      return(!all(x[5] == 0))
    } else if (g1 == 1 && g2 == 3 || g1 == 3 && g2 == 1) {
      return(!all(x[c(1, 5)] == 0))
    } else if (g1 == 2 && g2 == 3 || g1 == 3 && g2 == 2) {
      return(!all(x[1] == 0))
    }
  }

  return(FALSE)
}

#' Bayesian test for segregation distortion in tetraploids when genotypes are known.
#'
#' @section Impossible genotypes:
#' Some offspring genotype combinations are impossible given the parental
#' gentoypes. If these impossible genotypes combinations show up, we return a
#' log-BF of \code{-Inf} and missing values for all other return items.
#' The impossible genotypes are (when double reduction is allowed):
#' \describe{
#'   \item{\code{g1 == 0 && g2 == 0}}{Only offspring genotypes of 0 are possible.}
#'   \item{\code{g1 == 4 && g2 == 4}}{Only offspring genotypes of 4 are possible.}
#'   \item{\code{g1 == 0 && g2 == 4 || g1 == 4 && g2 == 0}}{Only offspring genotypes of 2 are possible.}
#'   \item{\code{g1 == 0 && g2 \%in\% c(1, 2, 3) || g1 \%in\% c(1, 2, 3) && g2 == 0}}{Only offspring genotypes of 0, 1, and 2 are possible.}
#'   \item{\code{g1 == 4 && g2 \%in\% c(1, 2, 3) || g1 \%in\% c(1, 2, 3) && g2 == 4}}{Only offspring genotypes of 2, 3, and 4 are possible.}
#' }
#' When double reduction is 0, the list of impossible genotypes is more
#' complicated.
#'
#' @section Unidentified parameters:
#' When \code{g1 = 2} or \code{g2 = 2} (or both), the model is not identified
#' and those estimates (\code{alpha}, \code{xi1}, and \code{xi2}) are
#' meaningless. Do NOT interpret them. The log-BF is fine, though.
#'
#' @inheritParams lrt_men_g4
#' @param ... Additional arguments to pass to Stan.
#'
#' @return A list with the following elements
#' \describe{
#'   \item{\code{lbf}}{The log Bayes Factor.}
#'   \item{\code{m0}}{Marginal likelihood under null.}
#'   \item{\code{ma}}{Marginal likelihood under alternative.}
#'   \item{\code{alpha}}{The posterior mean double reduction rate.}
#'   \item{\code{xi1}}{The posterior mean preferential pairing parameter of parent 1.}
#'   \item{\code{xi2}}{The posterior mean preferential pairing parameter of parent 2.}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(1)
#' gf <- offspring_gf_2(alpha = 1/12, xi1 = 0.1, xi2 = 0.5, p1 = 2, p2 = 2)
#' x <- offspring_geno(gf = gf, n = 20)
#' bayes_men_g4(x = x, g1 = 2, g2 = 2, chains = 1, iter = 1000)
#'
#' @export
bayes_men_g4 <- function(
    x,
    g1,
    g2,
    drbound = 1/6,
    shape1 = 5/9,
    shape2 = 10/9,
    pp = TRUE,
    dr = TRUE,
    ...) {
  ## check input
  stopifnot(length(x) == 5,
            x >= 0,
            g1 >= 0,
            g1 <= 4,
            g2 >= 0,
            g2 <= 4,
            drbound > 1e-6,
            drbound <= 1,
            is.logical(pp),
            is.logical(dr),
            length(g1) == 1,
            length(g2) == 1,
            length(drbound) == 1,
            length(pp) == 1,
            length(dr) == 1)

  ## log marginal likelihood under alternative
  ma <- marg_alt_g(x = x)

  ## Check if data are impossible under null
  if (is_impossible(x = x, g1 = g1, g2 = g2, dr = dr)) {
    ret <- list(
      lbf = -Inf,
      m0 = -Inf,
      ma = ma,
      alpha = NA_real_,
      xi1 = NA_real_,
      xi2 = NA_real_
    )
    return(ret)
  }

  ## log marginal likelihood under null
  if (pp && dr) {
    stout <- marg_f1_dr_pp_g4(x = x, g1 = g1, g2 = g2, drbound = drbound, shape1 = shape1, shape2 = shape2, output = "all", ...)
    m0 <- stout[[1]]
    alpha <- mean(as.data.frame(stout[[2]])$alpha)
    xi1 <- mean(as.data.frame(stout[[2]])$xi1)
    xi2 <- mean(as.data.frame(stout[[2]])$xi2)
  } else if (pp && !dr) {
    stout <- marg_f1_ndr_pp_g4(x = x, g1 = g1, g2 = g2, shape1 = shape1, shape2 = shape2, output = "all", ...)
    m0 <- stout[[1]]
    alpha <- 0
    xi1 <- mean(as.data.frame(stout[[2]])$gamma1)
    xi2 <- mean(as.data.frame(stout[[2]])$gamma2)
  } else if (!pp && dr) {
    stout <- marg_f1_dr_npp_g4(x = x, g1 = g1, g2 = g2, drbound = drbound, output = "all", ...)
    m0 <- stout[[1]]
    alpha <- mean(as.data.frame(stout[[2]])$alpha)
    xi1 <- 1/3
    xi2 <- 1/3
  } else if (!pp && !dr) {
    m0 <- marg_f1_ndr_npp_g4(x = x, g1 = g1, g2 = g2)
    alpha <- 0
    xi1 <- 1/3
    xi2 <- 1/3
  }
  if (g1 != 2) {
    xi1 <- NA_real_
  }
  if (g2 != 2) {
    xi2 <- NA_real_
  }

  ret <- list(
    lbf = m0 - ma,
    m0 = m0,
    ma = ma,
    alpha = alpha,
    xi1 = xi1,
    xi2 = xi2
  )

  return(ret)
}

#' Stan version of marg_alt_g(). Not to be used.
#'
#' @noRd
marg_alt_g_stan <- function(x, beta = rep(1, length(x)), lg = TRUE, ...) {
  ploidy <- length(x) - 1
  stopifnot(length(x) == length(beta))
  stan_dat <- list(K = ploidy, x = x, beta = beta)
  stan_out <- rstan::sampling(object = stanmodels$alt_g,
                              data = stan_dat,
                              verbose = FALSE,
                              show_messages = FALSE,
                              ...)
  bridge_out <- bridgesampling::bridge_sampler(stan_out, verbose = FALSE, silent = TRUE)

  if (lg) {
    mx <- bridge_out$logml
  } else {
    mx <- exp(bridge_out$logml)
  }
  return(mx)
}

#' Density of dirichlet-multinomial
#'
#' @param x the counts
#' @param alpha the dirichlet parameters
#' @param lg should we log marginal likelihood?
#'
#' @noRd
ddirmult <- function (x, alpha, lg = FALSE)
{
    stopifnot(length(x) == length(alpha))
    stopifnot(all(alpha > 0))
    asum <- sum(alpha)
    n <- sum(x)
    ll <- lgamma(asum) +
      lgamma(n + 1) -
      lgamma(n + asum) +
      sum(lgamma(x + alpha)) -
      sum(lgamma(alpha)) -
      sum(lgamma(x + 1))
    if (!lg) {
        ll <- exp(ll)
    }
    return(ll)
}

#' Marginal likelihood under alternative when genotypes are known
#'
#' @param x The genotype counts
#' @param beta The prior hyperparameters
#' @param lg A logical. Should we log the marginal likelihood or not?
#' @param ... Additional arguments to pass to \code{\link[rstan]{sampling}()}.
#'
#' @return The marginal likelihood.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L)
#' marg_alt_g(x = x, chains = 1)
#'
#' @noRd
marg_alt_g <- function(x, beta = rep(1, length(x)), lg = TRUE, ...) {
  ploidy <- length(x) - 1
  stopifnot(length(x) == length(beta))
  mx <- ddirmult(x = x, alpha = beta, lg = lg)
  return(mx)
}

#' Marginal likelihood, no double reduction, no preferential pairing, genotypes known.
#'
#' @inheritParams marg_f1_dr_pp_g4
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L) ## simulated via null
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_ndr_npp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' x <- c(25L, 24L, 20L, 15L, 16L) ## simulated via alt
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_ndr_npp_g4(x = x, g1 = g1, g2 = g2)
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' @noRd
marg_f1_ndr_npp_g4 <- function(x,
                               g1,
                               g2,
                               lg = TRUE) {
  gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = g1, p2 = g2)
  return(stats::dmultinom(x = x, prob = gf, log = lg))
}

#' Marginal likelihood, double reduction, no preferential pairing, genotypes known.
#'
#' Here, parental genotypes are known.
#'
#' @inheritParams marg_f1_dr_pp_g4
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L) ## simulated via null
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_dr_npp_g4(x = x, g1 = g1, g2 = g2) # -8.8
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' x <- c(25L, 24L, 20L, 15L, 16L) ## simulated via alt
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_dr_npp_g4(x = x, g1 = g1, g2 = g2) # -48.6
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' @noRd
marg_f1_dr_npp_g4 <- function(x,
                              g1,
                              g2,
                              drbound = 1/6,
                              mixprop = 0.001,
                              lg = TRUE,
                              output = c("marg", "all"),
                              ...) {
  stopifnot(length(x) == 5,
            length(g1) == 1,
            length(g2) == 1)
  stopifnot(g1 >= 0, g1 <= 4, g2 >= 0, g2 <= 4)
  stan_dat <- list(x = x,
                   drbound = drbound,
                   g1 = g1,
                   g2 = g2,
                   mixprop = mixprop)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_npp_g4,
                              data = stan_dat,
                              verbose = FALSE,
                              show_messages = FALSE,
                              ...)
  bridge_out <- bridgesampling::bridge_sampler(stan_out, verbose = FALSE, silent = TRUE)

  if (lg) {
    mx <- bridge_out$logml
  } else {
    mx <- exp(bridge_out$logml)
  }

  samps <- as.data.frame(stan_out)
  all <- list(mx, samps)

  output <- match.arg(output)
  if (output == "marg") {
    return(mx)
  } else {
    return(all)
  }

}

#' Marginal likelihood, no double reduction, preferential pairing, genotypes known.
#'
#' @inheritParams marg_f1_dr_pp_g4
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L) ## simulated via null
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_ndr_pp_g4(x = x, g1 = g1, g2 = g2) # -9.3
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' x <- c(25L, 24L, 20L, 15L, 16L) ## simulated via alt
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_ndr_pp_g4(x = x, g1 = g1, g2 = g2) # -40.8
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' @noRd
marg_f1_ndr_pp_g4 <- function(x,
                              g1,
                              g2,
                              shape1 = 5/9,
                              shape2 = 10/9,
                              mixprop = 0.001,
                              lg = TRUE,
                              output = c("marg", "all"),
                              ...) {
  stopifnot(length(x) == 5,
            length(g1) == 1,
            length(g2) == 1)
  stopifnot(g1 >= 0, g1 <= 4, g2 >= 0, g2 <= 4)
  stan_dat <- list(x = x,
                   g1 = g1,
                   g2 = g2,
                   mixprop = mixprop,
                   shape1 = shape1,
                   shape2 = shape2)
  stan_out <- rstan::sampling(object = stanmodels$marg_ndr_pp_g4,
                              data = stan_dat,
                              verbose = FALSE,
                              show_messages = FALSE,
                              ...)
  bridge_out <- bridgesampling::bridge_sampler(stan_out, verbose = FALSE, silent = TRUE)

  if (lg) {
    mx <- bridge_out$logml
  } else {
    mx <- exp(bridge_out$logml)
  }

  samps <- as.data.frame(stan_out)
  all <- list(mx, samps)

  output <- match.arg(output)
  if (output == "marg") {
    return(mx)
  } else {
    return(all)
  }
}

#' Marginal likelihood in tetraploid F1 population
#'
#' Here, parental genotypes are known.
#'
#' @param x The genotype counts of the offspring. \code{x[i]} is the
#'     number of offspring with genotype \code{i-1}.
#' @param g1 The first parent's genotype.
#' @param g2 The second parent's genotype.
#' @param mixprop The mixing proportion with the uniform for mixing purposes.
#' @param lg A logical. Should we log the marginal likelihood (\code{TRUE})
#'     or not (\code{FALSE})?
#' @param output Return either the marginal likelihood with (\code{"marg"})
#'    or the full STAN output with (\code{"all"}).
#' @param ... Additional arguments to pass to \code{\link[rstan]{sampling}()}.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(3L, 21L, 52L, 20L, 4L) ## simulated via null
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_dr_pp_g4(x = x, g1 = g1, g2 = g2) # -9.4
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' x <- c(25L, 24L, 20L, 15L, 16L) ## simulated via alt
#' g1 <- 2
#' g2 <- 2
#' mnull <- marg_f1_dr_pp_g4(x = x, g1 = g1, g2 = g2) # -32.8
#' malt <- -15.34 ## output of marg_alt_g(x = x)
#' mnull - malt ## log-BF
#'
#' @noRd
marg_f1_dr_pp_g4 <- function(x,
                             g1,
                             g2,
                             drbound = 1/6,
                             shape1 = 5/9,
                             shape2 = 10/9,
                             mixprop = 0.001,
                             lg = TRUE,
                             output = c("marg", "all"),
                             ...) {
  stopifnot(length(x) == 5,
            length(g1) == 1,
            length(g2) == 1)
  stopifnot(g1 >= 0, g1 <= 4, g2 >= 0, g2 <= 4)
  stan_dat <- list(x = x,
                   drbound = drbound,
                   g1 = g1,
                   g2 = g2,
                   mixprop = mixprop,
                   shape1 = shape1,
                   shape2 = shape2)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_pp_g4,
                              data = stan_dat,
                              verbose = FALSE,
                              show_messages = FALSE,
                              ...)
  bridge_out <- bridgesampling::bridge_sampler(stan_out, verbose = FALSE, silent = TRUE)

  if (lg) {
    mx <- bridge_out$logml
  } else {
    mx <- exp(bridge_out$logml)
  }

  samps <- as.data.frame(stan_out)
  all <- list(mx, samps)

  output <- match.arg(output)
  if (output == "marg") {
    return(mx)
  } else {
    return(all)
  }
}

#' Chi Square test when genotypes are known
#'
#' This chi-squared test is run under the assumption of no double reduction
#' and no preferential pairing.
#'
#' @param x Vector of observed genotype counts
#' @param g1 Parent 1's genotype
#' @param g2 Parent 2's genotype
#'
#' @return The Chi Square statistic and p-value
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(1, 2, 4, 3, 0)
#' g1 <- 2
#' g2 <- 2
#' chisq_ndr_npp_g4(x, g1, g2)
#'
#' x <- c(10, 25, 10, 0, 0)
#' g1 <- 1
#' g2 <- 1
#' chisq_ndr_npp_g4(x, g1, g2)
#'
#' @export
chisq_ndr_npp_g4 <- function(x, g1, g2){
  TOL <- sqrt(.Machine$double.eps)
  gf <- menbayes::offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = g1, p2 = g2)
  which_zero <- gf < TOL
  gf[which_zero] <- 0

  if (sum(x[which_zero]) > 0.5) { ## if any incompatibility, p-value is 0
    ret <- list(statistic = Inf,
                p_value = 0,
                df = NA_real_)
    return(ret)
  } else {
    x <- x[!which_zero]
    gf <- gf[!which_zero]
  }

  chout <- stats::chisq.test(x = x, p = gf)
  ret <- list(statistic = chout$statistic[[1]],
              p_value = chout$p.value[[1]],
              df = chout$parameter[[1]])
  return(ret)

}

#############################
## Methods using genotype log-likelihoods
#############################

#' Marginal likelihood under alternative using genotype log-likelihoods
#'
#' @inheritParams marg_f1_dr_pp_gl4
#' @param beta The vector of hyperparameters.
#' @param tol Change prior based on number of 0's
#'
#' @return The marginal likelihood.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' x <- c(1, 2, 7, 3, 3)
#' gl <- hwep::simgl(nvec = x, rdepth = 100, od = 0, bias = 1, seq = 0)
#' marg_alt_gl(gl = gl, chains = 1)
#'
#' @noRd
marg_alt_gl <- function(gl, beta = rep(1, 5), tol = 1e-2, lg = TRUE, ...) {

  ## hack to account for lots of zeros
  xest <- table(factor(apply(X = gl, MARGIN = 1, FUN = which.max) - 1, 1:ncol(gl) - 1))
  cond1 <- xest > nrow(gl) * tol
  cond2 <- xest > nrow(gl) * tol * 2
  beta[!cond2] <- 1/3
  gl <- gl[, cond1, drop = FALSE]
  beta <- beta[cond1]

  ploidy <- ncol(gl) - 1
  nind <- nrow(gl)
  stopifnot(ncol(gl) == length(beta))
  stan_dat <- list(N = nind, K = ploidy, gl = gl, beta = beta)
  stan_out <- rstan::sampling(object = stanmodels$alt_gl,
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

#' Marginal likelihood, no double reduction, no preferential pairing, genotype likelihoods.
#'
#' @inheritParams marg_f1_dr_pp_gl4
#' @param alpha The known fixed rate of double reduction.
#' @param xi1 The known rate of preferential pairing for parent 1. A value of 1/3
#'     corresponds to no preferential pairing.
#' @param xi2 The known rate of preferential pairing for parent 2. A value of 1/3
#'     corresponds to no preferential pairing.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(gf = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_npp_gl4(
#'     gl = gl,
#'     p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf)
#'     )
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_npp_gl4(
#'     gl = gl,
#'     p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf)
#'     )
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @noRd
marg_f1_ndr_npp_gl4 <- function(gl,
                                p1_gl = rep(log(0.2), 5),
                                p2_gl = rep(log(0.2), 5),
                                alpha = 0,
                                xi1 = 1/3,
                                xi2 = 1/3,
                                lg = TRUE,
                                output = c("marg", "all"), ...) {
  stopifnot(ncol(gl) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5)
  p1_gl <- p1_gl - updog::log_sum_exp(p1_gl)
  p2_gl <- p2_gl - updog::log_sum_exp(p2_gl)
  glmat <- matrix(data = NA_real_, nrow = 5, ncol = 5)
  for (i in 0:4) {
    for (j in 0:4) {
      gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = i, p2 = j)
      lgf <- log(gf)
      glmat[i + 1, j + 1] <- sum(apply(X = t(gl) + lgf, MARGIN = 2, FUN = updog::log_sum_exp)) +
        p1_gl[[i + 1]] +
        p2_gl[[j + 1]]
    }
  }
  mx <- updog::log_sum_exp(glmat)
  if (!lg) {
    mx <- exp(mx)
  }
  return(mx)
}

#' Marginal likelihood, double reduction, no preferential pairing, genotype likelihoods.
#'
#' @inheritParams marg_f1_dr_pp_gl4
#' @param xi1 The known rate of preferential pairing for parent 1. A value of 1/3
#'     corresponds to no preferential pairing.
#' @param xi2 The known rate of preferential pairing for parent 2. A value of 1/3
#'     corresponds to no preferential pairing.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' \dontrun{
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(gf = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_npp_gl4(
#'     gl = gl, p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     chains = 1
#'     )
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_npp_gl4(
#'     gl = gl,
#'     p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     chains = 1
#'     )
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#' }
#'
#' @noRd
marg_f1_dr_npp_gl4 <- function(gl,
                               p1_gl = rep(log(0.2), 5),
                               p2_gl = rep(log(0.2), 5),
                               drbound = 1/6,
                               xi1 = 1/3,
                               xi2 = 1/3,
                               ts1 = 1,
                               ts2 = 1,
                               mixprop = 0.001,
                               lg = TRUE,
                               output = c("marg", "all"),
                               ...) {
  stopifnot(ncol(gl) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5,
            length(mixprop) == 1)
  stopifnot(mixprop > 0, mixprop <= 1)
  p1_gl <- p1_gl - updog::log_sum_exp(p1_gl)
  p2_gl <- p2_gl - updog::log_sum_exp(p2_gl)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   drbound = drbound,
                   p1_gl = p1_gl,
                   p2_gl = p2_gl,
                   mixprop = mixprop,
                   ts1 = ts1,
                   ts2 = ts2,
                   xi1 = xi1,
                   xi2 = xi2)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_npp_gl4,
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

#' Marginal likelihood, no double reduction, preferential pairing, genotype likelihoods.
#'
#' @inheritParams marg_f1_dr_pp_gl4
#' @param alpha The known fixed rate of double reduction.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' \dontrun{
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(gf = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_pp_gl4(
#'     gl = gl,
#'     p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     chains = 1
#'     )
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_pp_gl4(
#'     gl = gl,
#'     p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     chains = 1
#'     )
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#' }
#'
#' @noRd
marg_f1_ndr_pp_gl4 <- function(gl,
                               p1_gl = rep(log(0.2), 5),
                               p2_gl = rep(log(0.2), 5),
                               alpha = 0,
                               shape1 = 5/9,
                               shape2 = 10/9,
                               mixprop = 0.001,
                               lg = TRUE,
                               output = c("marg", "all"),
                               ...) {
  stopifnot(ncol(gl) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5,
            length(mixprop) == 1)
  stopifnot(mixprop > 0, mixprop <= 1)
  p1_gl <- p1_gl - updog::log_sum_exp(p1_gl)
  p2_gl <- p2_gl - updog::log_sum_exp(p2_gl)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   p1_gl = p1_gl,
                   p2_gl = p2_gl,
                   mixprop = mixprop,
                   shape1 = shape1,
                   shape2 = shape2,
                   alpha = alpha)
  stan_out <- rstan::sampling(object = stanmodels$marg_ndr_pp_gl4,
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

#' Marginal likelihood, double reduction, preferential pairing, genotype likelihoods.
#'
#' @param gl A matrix of genotype log-likelihoods. The rows index the
#'     individuals and the columns index the possible genotypes. So
#'     \code{gl[i, k]} is the genotype log-likelihood for individual i and
#'     genotype k-1.
#' @param p1_gl The vector of genotype likelihoods of parent 1.
#' @param p2_gl The vector of genotype likelihoods of parent 2.
#' @param drbound The maximum rate of double reduction.
#' @param mixprop The mixing proportion with the uniform for mixing purposes.
#' @param lg A logical. Should we log the marginal likelihood (\code{TRUE})
#'     or not (\code{FALSE})?
#' @param output Should we return just the marginal likelihood (\code{"marg"})
#'     or both the marginal likelihood and all of the output of Stan
#'     (\code{"all"})?
#' @param ... Additional parameters sent to \code{\link[rstan]{sampling}()}.
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' \dontrun{
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(gf = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_pp_gl4(
#'     gl = gl,
#'     p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     chains = 1
#'     )
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_pp_gl4(
#'     gl = gl,
#'     p1_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     p2_gl = c(-Inf, -Inf, 0, -Inf, -Inf),
#'     chains = 1
#'     )
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#' }
#'
#' @noRd
marg_f1_dr_pp_gl4 <- function(gl,
                              p1_gl = rep(log(0.2), 5),
                              p2_gl = rep(log(0.2), 5),
                              drbound = 1/6,
                              shape1 = 5/9,
                              shape2 = 10/9,
                              ts1 = 1,
                              ts2 = 1,
                              mixprop = 0.001,
                              lg = TRUE,
                              output = c("marg", "all"),
                              ...) {
  stopifnot(ncol(gl) == 5,
            length(p1_gl) == 5,
            length(p2_gl) == 5,
            length(mixprop) == 1)
  stopifnot(mixprop > 0, mixprop <= 1)
  p1_gl <- p1_gl - updog::log_sum_exp(p1_gl)
  p2_gl <- p2_gl - updog::log_sum_exp(p2_gl)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   drbound = drbound,
                   p1_gl = p1_gl,
                   p2_gl = p2_gl,
                   mixprop = mixprop,
                   shape1 = shape1,
                   shape2 = shape2,
                   ts1 = ts1,
                   ts2 = ts2)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_pp_gl4,
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



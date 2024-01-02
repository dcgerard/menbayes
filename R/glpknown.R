###############
## Marginal likelihoods when parental genotypes are known, but offspring
## only have genotype likelihoods
###############

#' Bayesian test for segregation distortion in tetraploids using genotype likelihoods
#'
#' @section Unidentified parameters:
#' When \code{g1 = 2} or \code{g2 = 2} (or both), the model is not identified
#' and those estimates (\code{alpha}, \code{xi1}, and \code{xi2}) are
#' meaningless. Do NOT interpret them. The log-BF is fine, though.
#'
#' @inheritParams bayes_men_g4
#' @param gl The genotype log-likelihoods. The rows index the individuals
#'    and the columns index the genotypes.
#' @param g1 Either parent 1's genotype, or parent 1's genotype log-likelihoods.
#' @param g2 Either parent 2's genotype, or parent 2's genotype log-likelihoods.
#'
#' @examples
#' \dontrun{
#' ## null sims
#' set.seed(1)
#' g1 <- 1
#' g2 <- 2
#' gl <- simf1gl(n = 25, g1 = g1, g2 = g2)
#' bayes_men_gl4(gl = gl, g1 = g1, g2 = g2, chains = 1, iter = 1000)
#'
#' ## genotypes of parents not known
#' g1_gl <- log(c(0.2, 0.5, 0.2, 0.09, 0.01))
#' g2_gl <- log(c(0.09, 0.2, 0.4, 0.3, 0.01))
#' bayes_men_gl4(gl = gl, g1 = g1_gl, g2 = g2_gl, chains = 1, iter = 1000)
#'
#' ## alt sims
#' gl <- hwep::simgl(nvec = rep(5, 5), rdepth = 10)
#' bayes_men_gl4(gl = gl, g1 = g1, g2 = g2, chains = 1, iter = 1000)
#'
#' bayes_men_gl4(gl = gl, g1 = g1_gl, g2 = g2_gl, chains = 1, iter = 1000)
#' }
#'
#' @author David Gerard
#'
#' @export
bayes_men_gl4 <- function(
    gl,
    g1,
    g2,
    drbound = 1/6,
    shape1 = 5/9,
    shape2 = 10/9,
    ts1 = 1,
    ts2 = 1,
    pp = TRUE,
    dr = TRUE,
    alpha = 0,
    xi1 = 1/3,
    xi2 = 1/3,
    ...) {
  ## check input ----
  stopifnot(
    ncol(gl) == 5,
    drbound > 1e-6,
    drbound <= 1,
    is.logical(pp),
    is.logical(dr),
    length(drbound) == 1,
    length(pp) == 1,
    length(dr) == 1,
    alpha >= 0,
    alpha <= 1,
    xi1 >= 0,
    xi1 <= 1,
    xi2 >= 0,
    xi2 <= 1,
    length(alpha) == 1,
    length(xi1) == 1,
    length(xi2) == 1)
  if (length(g1) == 1 && length(g2) == 1) {
    pknown <- TRUE
    g1 <- round(g1)
    g2 <- round(g2)
    stopifnot(g1 >= 0, g1 <= 4, g2 >= 0, g2 <= 4)
  } else if (length(g1) == 5 && length(g2) == 5) {
    pknown <- FALSE
    ## normalize to sum to 1
    g1 <- g1 - log_sum_exp(g1)
    g2 <- g2 - log_sum_exp(g2)
  } else {
    stop("g1 and g2 should either both be length 1 (known genotypes)\nor both be length 5 (genotype log-likelihoods).")
  }

  ## log marginal likelihood under alternative ----
  ma <- marg_alt_gl(gl = gl, ...)

  ## log marginal likelihood under null ----
  if (pp && dr && pknown) {
    stout <- marg_f1_dr_pp_glpknown4(
      gl = gl,
      p1 = g1,
      p2 = g2,
      drbound = drbound,
      shape1 = shape1,
      shape2 = shape2,
      ts1 = ts1,
      ts2 = ts2,
      lg = TRUE,
      output = "all",
      ...)
    m0 <- stout[[1]]
    alpha <- mean(as.data.frame(stout[[2]])$alpha)
    xi1 <- mean(as.data.frame(stout[[2]])$xi1)
    xi2 <- mean(as.data.frame(stout[[2]])$xi2)
  } else if (pp && !dr && pknown) {
    stout <- marg_f1_ndr_pp_glpknown4(
      gl = gl,
      p1 = g1,
      p2 = g2,
      shape1 = shape1,
      shape2 = shape2,
      output = "all",
      alpha = alpha,
      ...)
    m0 <- stout[[1]]
    alpha <- alpha
    xi1 <- mean(as.data.frame(stout[[2]])$gamma1)
    xi2 <- mean(as.data.frame(stout[[2]])$gamma2)
  } else if (!pp && dr && pknown) {
    stout <- marg_f1_dr_npp_glpknown4(
      gl = gl,
      p1 = g1,
      p2 = g2,
      drbound = drbound,
      ts1 = ts1,
      ts2 = ts2,
      output = "all",
      xi1 = xi1,
      xi2 = xi2,
      ...)
    m0 <- stout[[1]]
    alpha <- mean(as.data.frame(stout[[2]])$alpha)
    xi1 <- xi1
    xi2 <- xi2
  } else if (!pp && !dr && pknown) {
    m0 <- marg_f1_ndr_npp_glpknown4(gl = gl, p1 = g1, p2 = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
    alpha <- alpha
    xi1 <- xi1
    xi2 <- xi2
  } else if (pp && dr && !pknown) {
    stout <- marg_f1_dr_pp_gl4(
      gl = gl,
      p1_gl = g1,
      p2_gl = g2,
      drbound = drbound,
      shape1 = shape1,
      shape2 = shape2,
      ts1 = ts1,
      ts2 = ts2,
      lg = TRUE,
      output = "all",
      ...)
    m0 <- stout[[1]]
    alpha <- mean(as.data.frame(stout[[2]])$alpha)
    xi1 <- mean(as.data.frame(stout[[2]])$xi1)
    xi2 <- mean(as.data.frame(stout[[2]])$xi2)
  } else if (pp && !dr && !pknown) {
    stout <- marg_f1_ndr_pp_gl4(
      gl = gl,
      p1_gl = g1,
      p2_gl = g2,
      shape1 = shape1,
      shape2 = shape2,
      output = "all",
      alpha = alpha,
      ...)
    m0 <- stout[[1]]
    alpha <- alpha
    xi1 <- mean(as.data.frame(stout[[2]])$gamma1)
    xi2 <- mean(as.data.frame(stout[[2]])$gamma2)
  } else if (!pp && dr && !pknown) {
    stout <- marg_f1_dr_npp_gl4(
      gl = gl,
      p1_gl = g1,
      p2_gl = g2,
      drbound = drbound,
      ts1 = ts1,
      ts2 = ts2,
      output = "all",
      xi1 = xi1,
      xi2 = xi2,
      ...)
    m0 <- stout[[1]]
    alpha <- mean(as.data.frame(stout[[2]])$alpha)
    xi1 <- xi1
    xi2 <- xi2
  } else if (!pp && !dr && !pknown) {
    m0 <- marg_f1_ndr_npp_gl4(gl = gl, p1_gl = g1, p2_gl = g2, alpha = alpha, xi1 = xi1, xi2 = xi2)
    alpha <- alpha
    xi1 <- xi1
    xi2 <- xi2
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

#' Marginal likelihood, no double reduction, no preferential pairing, parent genotypes known, offspring genotypes unknown
#'
#' @inheritParams marg_f1_dr_pp_glpknown4
#' @param alpha The known fixed rate of double reduction.
#' @param xi1 The known rate of preferential pairing for parent 1. A value of 1/3
#'     corresponds to no preferential pairing.
#' @param xi2 The known rate of preferential pairing for parent 2. A value of 1/3
#'     corresponds to no preferential pairing.
#'
#' @examples
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(gf = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_npp_glpknown4(gl = gl, p1 = 2, p2 = 2)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_npp_glpknown4(gl = gl, p1 = 2, p2 = 2)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @author Mira Thakkar and David Gerard
#'
#' @noRd
marg_f1_ndr_npp_glpknown4 <- function(gl,
                                      p1,
                                      p2,
                                      alpha = 0,
                                      xi1 = 1/3,
                                      xi2 = 1/3,
                                      lg = TRUE) {
  gf <- offspring_gf_2(alpha = alpha, xi1 = xi1, xi2 = xi2, p1 = p1, p2 = p2)
  lgf <- log(gf)
  mx <- sum(apply(X = t(gl) + lgf, MARGIN = 2, FUN = updog::log_sum_exp))
  if(!lg) {
    mx <- exp(mx)
  }
  return(mx)
}

#' Marginal likelihood, double reduction, no preferential pairing, parent genotypes known, offspring genotypes unknown
#'
#' @inheritParams marg_f1_dr_pp_glpknown4
#' @param xi1 The known rate of preferential pairing for parent 1. A value of 1/3
#'     corresponds to no preferential pairing.
#' @param xi2 The known rate of preferential pairing for parent 2. A value of 1/3
#'     corresponds to no preferential pairing.
#'
#' @examples
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(gf = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_npp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_npp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @author Mira Thakkar and David Gerard
#'
#' @noRd
marg_f1_dr_npp_glpknown4 <- function(gl,
                                     p1,
                                     p2,
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
            length(p1) == 1,
            length(p2) == 1)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   drbound = drbound,
                   g1 = p1,
                   g2 = p2,
                   mixprop = mixprop,
                   ts1 = ts1,
                   ts2 = ts2,
                   xi1 = xi1,
                   xi2 = xi2)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_npp_glpknown4,
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


#' Marginal likelihood, no double reduction, preferential pairing, parent genotypes known, offspring genotypes unknown
#'
#' @inheritParams marg_f1_dr_pp_glpknown4
#' @param alpha The known fixed rate of double reduction.
#'
#' @examples
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(gf = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_pp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_ndr_pp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @author Mira Thakkar and David Gerard
#'
#' @noRd
marg_f1_ndr_pp_glpknown4 <- function(gl,
                                     p1,
                                     p2,
                                     alpha = 0,
                                     shape1 = 5/9,
                                     shape2 = 10/9,
                                     mixprop = 0.001,
                                     lg = TRUE,
                                     output = c("marg", "all"),
                                     ...) {
  stopifnot(ncol(gl) == 5,
            length(p1) == 1,
            length(p2) == 1)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   g1 = p1,
                   g2 = p2,
                   mixprop = mixprop,
                   shape1 = shape1,
                   shape2 = shape2,
                   alpha = alpha)
  stan_out <- rstan::sampling(object = stanmodels$marg_ndr_pp_glpknown4,
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

#' Marginal likelihood, double reduction, preferential pairing, parent genotypes known, offspring genotypes unknown
#'
#' @param gl A matrix of offspring genotype log-likelihoods. The rows index the
#'     individuals and the columns index the possible genotypes. So
#'     \code{gl[i, k]} is the offspring genotype log-likelihood for individual i and
#'     genotype k-1.
#' @param p1 Genotype of parent 1.
#' @param p2 Genotype of parent 2.
#' @param mixprop The mixing proportion with the uniform for mixing purposes.
#' @param lg A logical. Should we log the marginal likelihood (\code{TRUE}) or
#'     not (\code{FALSE})?
#' @param output Should we return just the marginal likelihood (\code{"marg"})
#'     or both the marginal likelihood and all of the output of Stan
#'     (\code{"all"})?
#' @param drbound The maximum rate of double reduction.
#' @param ... Additional parameters sent to \code{\link[rstan]{sampling}()}.
#'
#' @examples
#' ## null sims
#' set.seed(1)
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = 2, p2 = 2)
#' gcount <- offspring_geno(gf = gf, n = 20)
#' gvec <- gcount_to_gvec(gcount)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_pp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -48.81 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' ## alt sims
#' set.seed(1)
#' gvec <- rep(0:4, each = 4)
#' fout <- po_gl(genovec = gvec, p1_geno = 2, p2_geno = 2, ploidy = 4)
#' gl <- fout$genologlike
#' mnull <- marg_f1_dr_pp_glpknown4(gl = gl, p1 = 2, p2 = 2, chains = 1)
#' malt <- -46.82 ## result of running marg_alt_gl(gl)
#' mnull - malt ## log-BF
#'
#' @author Mira Thakkar and David Gerard
#'
#' @noRd
marg_f1_dr_pp_glpknown4 <- function(gl,
                                    p1,
                                    p2,
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
            length(p1) == 1,
            length(p2) == 1)
  stan_dat <- list(gl = gl,
                   N = nrow(gl),
                   drbound = drbound,
                   g1 = p1,
                   g2 = p2,
                   mixprop = mixprop,
                   shape1 = shape1,
                   shape2 = shape2,
                   ts1 = ts1,
                   ts2 = ts2)
  stan_out <- rstan::sampling(object = stanmodels$marg_dr_pp_glpknown4,
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


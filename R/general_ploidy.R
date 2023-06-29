
#' Marginal likelihood for F1 populations of any ploidy.
#'
#' @inheritParams marg_f1_dr_pp_g4
#' @param beta The concentration parameters for the dirichlet of the gamere
#'     frequencies.
#'
#' @author David Gerard
#'
#' @examples
#' x <- c(4, 8, 12, 16, 12, 8, 4) ## generated under null
#' mnull <- marg_f1_g(x = x)
#' malt <- marg_alt_g(x = x)
#' mnull - malt
#'
#' x <- rep(9, 7) ## generated under alt
#' mnull <- marg_f1_g(x = x)
#' malt <- marg_alt_g(x = x)
#' mnull - malt
#'
marg_f1_g <- function(x,
                      beta = rep(1, (length(x) + 1) / 2),
                      mixprop = 0.001,
                      lg = TRUE,
                      ...) {
  ploidy <- length(x) - 1
  stopifnot(ploidy / 2 + 1 == length(beta))
  stan_dat <- list(ploidy = ploidy,
                   phalf = ploidy / 2 + 1,
                   x = x,
                   mixprop = mixprop,
                   beta = beta)
  stan_out <- rstan::sampling(object = stanmodels$marg_f1_g,
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

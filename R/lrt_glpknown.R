## Offspring genotype likelihoods, parent genotypes

#' Likelihood under three parameter model when using offspring genotypes
#' likelihoods but parent genotypes are known.
#'
#' @inheritParams like_gknown
#' @param gl The matrix of genotype likelihoods of the offspring. Rows index
#'     The individuals, columns index the genotypes.
#'
#' @author David Gerard
#'
#' @examples
#' p1 <- 1
#' p2 <- 0
#' gf <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = p1, p2 = p2)
#' x <- offspring_geno(gf = gf, n = 10)
#' genovec <- gcount_to_gvec(gcount = x)
#' uout <- po_gl(genovec = genovec, p1_geno = p1, p2_geno = p2, ploidy = 4)
#' gl <- uout$genologlike
#' like_glpknown(
#'   gl = gl,
#'   tau = 1/2,
#'   beta = 1/12,
#'   gamma1 = 1/3,
#'   gamma2 = 1/3,
#'   g1 = p1,
#'   g2 = p2,
#'   log_p = TRUE)
#'
#' @export
like_glpknown <- function(gl, tau, beta, gamma1, gamma2, g1, g2, log_p = TRUE) {
  stopifnot(ncol(gl) == 5)
  stopifnot(tau >= 0, tau <= 1,
            beta >= 0, beta <= 1,
            gamma1 >= 0, gamma1 <= 1,
            g1 >= 0, g1 <= 4,
            g2 >= 0, g2 <= 4)

  gf <- offspring_gf_3(
    tau = tau,
    beta = beta,
    gamma1 = gamma1,
    gamma2 = gamma2,
    p1 = g1,
    p2 = g2)

  return(llike_li(B = gl, lpivec = log(gf)))
}


#' Likelihood ratio test assuming no double reduction and no preferential pairing
#'
#' This is when offspring genotypes are not known
#'
#' @inheritParams like_glpknown
#'
#' @return A list of length three with the following elements
#' \describe{
#'   \item{\code{statistic}}{The likelihood ratio test statistic.}
#'   \item{\code{df}}{The degrees of freedom.}
#'   \item{\code{p_value}}{The p-value.}
#' }
#'
#' @author David Gerard
lrt_ndr_npp_glpknown4 <- function(gl, g1, g2) {
  ## MLE under alternative
  log_qhat1 <- c(em_li(B = gl))
  l1 <- llike_li(B = gl, lpivec = log_qhat1)

  ## null likelihood
  qhat2 <- offspring_gf_2(alpha = 0, xi1 = 1/3, xi2 = 1/3, p1 = g1, p2 = g2)
  l0 <- llike_li(B = gl, lpivec = log(qhat2))

  ## calculate test statistic and run test
  llr <- -2 * (l0 - l1)
  df <- 4
  p_value <- stats::pchisq(q = llr, df = df, lower.tail = FALSE, log.p = FALSE)

  return(list(llr = llr, p_value = p_value, df = df))
}

## LRT when dr and pp are estimated --------------------------------------

#' LRT when both double reduction and preferential pairing are not known.
#'
#' Offspring use genotype likelihoods, parent genotypes are known.
#'
#' @inheritParams like_glpknown
#' @param drbound Maximum value of double reduction.
#'
#' @author David Gerard
lrt_dr_pp_glpknown4 <- function(gl, g1, g2, drbound = 1/6) {

}

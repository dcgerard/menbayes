#' Tetraploid genotype frequencies of gametes when one parent's genotype is known
#'
#' @param alpha The double reduction rate
#' @param xi The preferential pairing parameter
#' @param ell The parental genotype
#'
#' @return The gamete genotype frequencies
#'
#' @author Mira Thakkar and David Gerard
#'
#' @examples
#' alpha <- 1/6
#' xi <- 1/3
#' pvec_tet(alpha = alpha, xi = xi, ell = 0)
#' pvec_tet(alpha = alpha, xi = xi, ell = 1)
#' pvec_tet(alpha = alpha, xi = xi, ell = 2)
#' pvec_tet(alpha = alpha, xi = xi, ell = 3)
#' pvec_tet(alpha = alpha, xi = xi, ell = 4)
#'
#' @export
pvec_tet <- function(alpha, xi, ell) {
  if (ell > 4 | ell < 0 | is.na(ell)){
    stop("Invalid input")
  } else if(ell == 0){
    pv <- c(
      1,
      0,
      0
    )
  } else if(ell == 1){
    pv <- c(
      0.5 + 0.25 * alpha,
      0.5 - 0.5 * alpha,
      0.25 * alpha
    )
  } else if(ell == 2){
    pv <- c(
      0.5 * alpha + 0.25 * (1 - alpha) * (1 - xi),
      0.5 * (1 - alpha) * (1 + xi),
      0.5 * alpha + 0.25 * (1 - alpha) * (1 - xi)
    )
  } else if(ell == 3){
    pv <- c(
      0.25 * alpha,
      0.5 - 0.5 * alpha,
      0.5 + 0.25 * alpha
    )
  } else if(ell == 4){
    pv <- c(
      0,
      0,
      1
    )
  }
  return(pv)
}


#' Hexaploid genotype frequencies of gametes when one parent's genotype is known
#'
#' This is experimental and should not be used right now.
#'
#' @inheritParams pvec_tet
#'
#' @author David Gerard
#'
#' @examples
#' alpha <- 3/10
#' xi <- 0.5
#' pvec_hex(alpha = alpha, xi = xi, ell = 0)
#' pvec_hex(alpha = alpha, xi = xi, ell = 1)
#' pvec_hex(alpha = alpha, xi = xi, ell = 2)
#' pvec_hex(alpha = alpha, xi = xi, ell = 3)
#' pvec_hex(alpha = alpha, xi = xi, ell = 4)
#' pvec_hex(alpha = alpha, xi = xi, ell = 5)
#' pvec_hex(alpha = alpha, xi = xi, ell = 6)
#'
#'
#' @noRd
pvec_hex <- function(alpha, xi, ell) {
  if (ell > 6 | ell < 0 | is.na(ell)){
    stop("Invalid input")
  } else if (ell == 0) {
    pv <- c(1, 0, 0, 0)
  } else if (ell == 1) {
    pv <- c(
      2 * alpha / 3 + (1 - alpha) / 2,
      alpha / 6 + (1 - alpha) / 2,
      alpha / 6,
      0
    )
  } else if (ell == 2) {
    pv <- c(
      2 * alpha / 5 + (1 - xi) * (1 - alpha) / 4,
      4 * alpha / 15 + (xi + (1 - xi) / 2) * (1 - alpha),
      4 * alpha / 15 + (1 - xi) * (1 - alpha) / 4,
      alpha / 15
    )
  } else if (ell == 3) {
    pv <- c(
      alpha / 5 + (1 - xi) * (1 - alpha) / 8,
      3 * alpha / 10 + (xi / 2 + 3 * (1 - xi) / 8) * (1 - alpha),
      3 * alpha / 10 + (xi / 2 + 3 * (1 - xi) / 8) * (1 - alpha),
      alpha / 5 + (1 - xi) * (1 - alpha) / 8
    )
  } else if (ell == 4) {
    pv <- c(
      alpha / 15,
      4 * alpha / 15 + (1 - xi) * (1 - alpha) / 4,
      4 * alpha / 15 + (xi + (1 - xi) / 2) * (1 - alpha),
      2 * alpha / 5 + (1 - xi) * (1 - alpha) / 4
    )
  } else if (ell == 5) {
    pv <- c(
      0,
      alpha / 6,
      alpha / 6 + (1 - alpha) / 2,
      2 * alpha / 3 + (1 - alpha) / 2
    )
  } else if (ell == 6) {
    pv <- c(0, 0, 0, 1)
  }
  return(pv)
}

#' Function that takes as input the double reduction rate, the preferential
#' pairing rate, and parent genotypes to return zygote genotype frequencies.
#'
#' @param alpha The double reduction rate
#' @param xi The preferential pairing parameter
#' @param p1 The first parent's genotype
#' @param p2 The second parent's genotype
#'
#' @return Zygote genotype frequencies
#'
#' @author Mira Thakkar
#'
#' @examples
#' alpha <- 1/6
#' xi <- 1/3
#' p1 <- 2
#' p2 <- 3
#' offspring_gf(alpha = alpha, xi = xi, p1 = p1, p2 = p2)
#'
#' @export
offspring_gf <- function(alpha, xi, p1, p2){

  pvec1 <- pvec_tet(alpha = alpha, xi = xi, ell = p1)
  pvec2 <- pvec_tet(alpha = alpha, xi = xi, ell = p2)

  qvec <- stats::convolve(pvec1, rev(pvec2), type = "open")

  stopifnot(qvec > -1e-06)
  stopifnot(abs(sum(qvec) - 1) < 1e-06)
  qvec[qvec < 0] <- 0
  qvec <- qvec / sum(qvec)

  return(qvec)
}


#' A function called offspring_geno() Which takes as input the offspring
#' genotype frequencies and a sample size and returns simulated genotypes.
#'
#' @param x Vector of offspring genotype frequencies
#' @param n Sample size
#'
#' @return Simulated genotypes
#'
#' @author Mira Thakkar
#'
#' @examples
#' x <- offspring_gf(alpha = 1/6, xi = 1/3, p1 = 2, p2 = 3)
#' offspring_geno(x = x, n = 10)
#'
#' @export
offspring_geno <- function(x, n){
  sim_gen <- c(stats::rmultinom(n = 1, size = n, prob = x))
  return(sim_gen)
}

#' Converts genotype counts to genotype vectors.
#'
#' @param gcount The vector of genotype counts.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' gcount <- c(1, 2, 3, 0, 5)
#' gcount_to_gvec(gcount = gcount)
gcount_to_gvec <- function(gcount) {
  unlist(mapply(FUN = rep, x = seq_along(gcount) - 1, each = gcount))
}

#' Function which takes as input (i) the parent genotypes,
#' (ii) the offspring genotype freq, (iii) sequencing error rate, (iv) read
#' depth, (v) bias, (vi) overdispersion and returns genotype likelihoods.
#'
#' @param genovec Offspring genotypes
#' @param p1_geno Parent 1 genotype
#' @param p2_geno Parent 2 genotype
#' @param ploidy Ploidy
#' @param seq Sequencing error rate
#' @param rd Read depth
#' @param bias Bias
#' @param od Overdispersion
#'
#' @return Genotype likelihoods
#'
#' @author Mira Thakkar
#'
#' @export
po_gl <- function(genovec, p1_geno, p2_geno, ploidy, seq = 0.01, rd = 10, bias = 1, od = 0.01) {
  n <- length(genovec)
  sizevec <- rep(rd, length.out = n)
  refvec <- updog::rflexdog(sizevec = sizevec, geno = genovec, ploidy = ploidy, seq = seq, bias = bias, od = od)
  p1ref <- updog::rflexdog(sizevec = rd, geno = p1_geno, ploidy = ploidy, seq = seq, bias = bias, od = od)
  p2ref <- updog::rflexdog(sizevec = rd, geno = p2_geno, ploidy = ploidy, seq = seq, bias = bias, od = od)

  fout <- updog::flexdog_full(refvec = refvec,
                              sizevec = sizevec,
                              ploidy = ploidy,
                              model = "f1pp",
                              seq = seq,
                              bias = bias,
                              od = od,
                              update_bias = FALSE,
                              update_seq = FALSE,
                              update_od = FALSE,
                              p1ref = p1ref,
                              p1size = rd,
                              p2ref = p2ref,
                              p2size = rd)

  return(fout)
}

## polymapR's version of checking for segregation distortion


#' Run segregation distortion tests as implemented in the polymapR package.
#'
#' The polymapR tests for segregation distortion by iterating through all
#' possible forms of disomic or polysomic inheritance from either parent,
#' tests for concordance of the offspring genotypes using a chi-squared
#' test, and returns the largest p-value. It sometimes chooses a different
#' p-value based on other heuristics. When \code{type = "menbayes"},
#' we only look at patterns of the given parent genotypes, choosing
#' the largest p-value. When \code{type = "polymapR"}, we return what
#' they use via their heuristics.
#'
#' @param g Either a vector of genotype counts, or a matrix of genotype
#'     likelihoods where the rows index the individuals and the columns
#'     index the genotypes.
#' @param type Either my implementation which approximates that of
#'     polymapR (\code{"menbayes"}) or the exact implementation
#'     through polymapR (\code{"polymapR"}). Note that
#'     polymapR needs to be installed for \code{type = "polymapR"}.
#'
#' @author David Gerard
#'
#' @export
polymapr_test <- function(g, type = c("menbayes", "polymapR")) {
  if (!requireNamespace("polymapR", quietly = TRUE) && type == "polymapR") {
    stop(
      paste0(
        "polymapr_test:\n",
        "polymapR needs to be installed to use type = 'polymapR'.\n",
        "You can install it with install.packages('polymapR')"
        )
      )
  }

}


#' polymapR test when genotypes are known.
#'
#' @param x genotype count vector
#' @param g1 parent 1's gentoype
#' @param g2 parent 2's genotype
#'
#' @author David Gerard
#'
#' @examples
#' x <- c(3, 22, 50, 22, 3)
#' polymapR_package_g(x = x, g1 = 2, g2 = 2)
#'
#' @noRd
polymapR_package_g <- function(x, g1, g2) {
  stopifnot(requireNamespace("polymapR", quietly = TRUE))
  ploidy <- length(x) - 1
  stopifnot(
    g1 >= 0,
    g2 >= 0,
    g1 <= ploidy,
    g2 <= ploidy
  )

  df <- matrix(c(g1, g2, gcount_to_gvec(gcount = x)), nrow = 1)
  fnames <- paste0("F", seq_len(ncol(df) - 2))
  colnames(df) <- c("P1", "P2", fnames)

  df <- rbind(df, 1) ## they forgot a drop = FALSE somewhere

  rownames(df) <- c("M1", "M2")

  cout <- polymapR::checkF1(
    input_type = "discrete",
    dosage_matrix = df,
    parent1 = "P1",
    parent2 = "P2",
    F1 = fnames,
    polysomic = TRUE,
    disomic = TRUE,
    mixed = TRUE,
    ploidy = ploidy)

  ret <- list(
    p_value = cout$checked_F1$Pvalue_bestParentfit[[1]],
    bestfit = as.character(cout$checked_F1$bestParentfit[[1]]),
    frq_invalid = cout$checked_F1$frqInvalid_bestParentfit[[1]]
    )

  return(ret)
}

#' test for segregation distortion, approximating polymapR procedure
#'
#' @inheritParams polymapR
#' @param seg_invalidrate If there is only one class possible, the p-value
#'     is the binomial probability of the invalid number against a true
#'     invalid rate of this, defaults to 0.03.
#'
#' @author David Gerard
#'
#' @noRd
polymapr_approx_g <- function(x, g1, g2, seg_invalidrate = 0.03) {
  ploidy <- 4
  stopifnot(
    length(x) == 5,
    x >= 0,
    g1 >= 0,
    g2 >= 0,
    g1 <= ploidy,
    g2 <= ploidy
  )

  pval <- 0
  bi <- NULL
  frq_invalid <- NULL
  for (i in seq_len(nrow(segtypes))) {
    is_p <- any((segtypes$pardosage[[i]][, 1] == g1) &
                  (segtypes$pardosage[[i]][, 2] == g2))
    if (is_p) {
      fq <- segtypes$freq[[i]]
      not_0 <- fq > sqrt(.Machine$double.eps)
      if (sum(not_0) == 1) {
        chout <- list(p.value = stats::pbinom(q = sum(x[!not_0]), size = sum(x), prob = seg_invalidrate, lower.tail = FALSE))
      } else {
        suppressWarnings(
          chout <- stats::chisq.test(x = x[not_0], p = fq[not_0])
        )
      }
      if (pval < chout$p.value) {
        pval <- chout$p.value
        bi <- i
        frq_invalid <- sum(x[!not_0])
      }
    }
  }

  ret <- list(
    p_value = pval,
    best_fit = segtypes$mod[[i]],
    frq_invalid = frq_invalid
  )
  return(ret)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lrt_gknown.R
\name{lrt_men_g4}
\alias{lrt_men_g4}
\title{Likelihood ratio test for segregation distortion with known genotypes}
\usage{
lrt_men_g4(
  x,
  g1,
  g2,
  drbound = 1/6,
  pp = TRUE,
  dr = TRUE,
  alpha = 0,
  xi1 = 1/3,
  xi2 = 1/3
)
}
\arguments{
\item{x}{A vector of genotype counts. \code{x[i]} is the number of
offspring with genotype \code{i-1}.}

\item{g1}{The genotype of parent 1.}

\item{g2}{The genotype of parent 2.}

\item{drbound}{The maximum rate of double reduction. A default of 1/6
is provided, which is the rate under the complete equational
segregation model of meiosis.}

\item{pp}{A logical. Should we account for preferential pairing
(\code{TRUE}) or not (\code{FALSE})?}

\item{dr}{A logical. Should we account for double reduction
(\code{TRUE}) or not (\code{FALSE})?}

\item{alpha}{If \code{dr = FALSE}, this is the known rate of double
reduction.}

\item{xi1}{If \code{pp = FALSE}, this is the known preferential pairing
parameter of parent 1.}

\item{xi2}{If \code{pp = FALSE}, this is the known preferential pairing
parameter of parent 2.}
}
\value{
A list with the following elements
\describe{
\item{\code{statistic}}{The log-likelihood ratio test statistic.}
\item{\code{df}}{The degrees of freedom.}
\item{\code{p_value}}{The p-value.}
\item{\code{alpha}}{The estimated double reduction rate.}
\item{\code{xi1}}{The estimated preferential pairing parameter of parent 1.}
\item{\code{xi2}}{The estimated preferential pairing parameter of parent 2.}
}
}
\description{
This will run a likelihood ratio test using the genotypes of an F1 population
of tetraploids for the null of Mendelian segregation (accounting for double
reduction and preferential pairing) against the alternative of
segregation distortion. This is when the genotypes are assumed known.
}
\section{Impossible genotypes}{

Some offspring genotype combinations are impossible given the parental
gentoypes. If these impossible genotypes combinations show up, we return a
p-value of 0, a log-likelihood ratio statistic of Infinity, and missing
values for all other return items. The impossible genotypes are:
\describe{
\item{\code{g1 = 0 && g2 = 0}}{Only offspring genotypes of 0 are possible.}
\item{\code{g1 = 4 && g2 = 4}}{Only offspring genotypes of 4 are possible.}
\item{\code{g1 = 0 && g2 = 4 || g1 == 4 && g2 == 0}}{Only offspring genotypes of 2 are possible.}
\item{\code{g1 = 0 && g2 \%in\% c(1, 2, 3) || g1 = \%in\% c(1, 2, 3) && g2 == 0}}{Only offspring genotypes of 0, 1, and 2 are possible.}
\item{\code{g1 = 4 && g2 \%in\% c(1, 2, 3) || g1 = \%in\% c(1, 2, 3) && g2 == 4}}{Only offspring genotypes of 2, 3, and 4 are possible.}
}
}

\section{Unidentified parameters}{

When \code{g1 = 2} or \code{g2 = 2} (or both), the model is not identified
and those estimates (\code{alpha}, \code{xi1}, and \code{xi2}) are
meaningless. Do NOT interpret them.

The estimate of \code{alpha} (double reduction rate) IS identified as
long as at least one parent is simplex, and no parent is duplex.
However, the estimates of the double reduction rate have extremely high
variance.
}

\examples{
set.seed(100)
gf <- offspring_gf_2(alpha = 1/12, xi1 = 0.2, xi2 = 0.6, p1 = 1, p2 = 0)
x <- offspring_geno(gf = gf, n = 100)
lrt_men_g4(x = x, g1 = 1, g2 = 0)

}
\author{
David Gerard
}

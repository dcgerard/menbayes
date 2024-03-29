% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rscripts.R
\name{pvec_tet_2}
\alias{pvec_tet_2}
\title{Tetraploid gamete frequencies of gametes when one parent's genotype is known}
\usage{
pvec_tet_2(alpha, xi, ell)
}
\arguments{
\item{alpha}{The double reduction rate}

\item{xi}{The preferential pairing parameter}

\item{ell}{The parental genotype}
}
\value{
The gamete genotype frequencies
}
\description{
This is under the two parameter model.
}
\examples{
alpha <- 1/6
xi <- 1/3
pvec_tet_2(alpha = alpha, xi = xi, ell = 0)
pvec_tet_2(alpha = alpha, xi = xi, ell = 1)
pvec_tet_2(alpha = alpha, xi = xi, ell = 2)
pvec_tet_2(alpha = alpha, xi = xi, ell = 3)
pvec_tet_2(alpha = alpha, xi = xi, ell = 4)

}
\author{
Mira Thakkar and David Gerard
}

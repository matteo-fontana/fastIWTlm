\name{fastIWT2}
\alias{fastIWT2}
\title{Fast Two-Sample Interval-Wise Testing (C++/Rcpp engine)}
\description{
  RcppEigen + OpenMP implementation of the two-population Interval-Wise
  Test. Computes pointwise and interval-wise permutation p-values together
  with the interval-wise adjusted p-value function.
}
\usage{
fastIWT2(Data1, Data2, Mu, B = 1000, alt = "two.sided", maxrow = 0,
         paired = FALSE, recycle = FALSE, MPthreads = 1)
}
\arguments{
  \item{Data1, Data2}{numeric matrices (statistical units in rows,
    abscissa points in columns) for the two samples.}
  \item{Mu}{numeric vector of length \code{ncol(Data1)}: the mean under the
    null hypothesis, subtracted from \code{Data1}.}
  \item{B}{number of permutations.}
  \item{alt}{alternative hypothesis: one of \code{"two.sided"},
    \code{"greater"}, \code{"less"}.}
  \item{maxrow}{truncation of the interval combination (0 = full
    combination).}
  \item{paired}{logical; use a paired (sign-flip) permutation scheme.}
  \item{recycle}{logical; use circular (recycled) intervals.}
  \item{MPthreads}{number of OpenMP threads.}
}
\value{
  A list with \code{T0_stat} (observed statistics), \code{unadj_pval}
  (pointwise p-values), \code{inter} (interval p-value matrix) and
  \code{adj_pval} (interval-wise adjusted p-values).
}
\seealso{\code{\link{IWT2}} for the pure-R implementation.}
\examples{
\dontrun{
load(system.file("extdata", "test-data.rdata", package = "fastIWTlm"))
res <- fastIWT2(data1, data2, mu, B = 1000, MPthreads = 2)
}
}

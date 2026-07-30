\name{IWT2}
\alias{IWT2}
\title{Two-Sample Interval-Wise Testing (pure R)}
\description{
  Pure-R implementation of the two-population Interval-Wise Test. Provided
  as a readable reference; for large problems use the compiled
  \code{\link{fastIWT2}} engine.
}
\usage{
IWT2(data1, data2, mu = 0, B = 1000, maxrow = 1, paired = FALSE,
     recycle = FALSE, alternative = "two.sided")
}
\arguments{
  \item{data1, data2}{numeric matrices (units in rows, abscissa points in
    columns) for the two samples.}
  \item{mu}{the mean under the null hypothesis, subtracted from
    \code{data1} (scalar or length-\code{p} vector).}
  \item{B}{number of permutations.}
  \item{maxrow}{first row of the interval-combination matrix to compute.}
  \item{paired}{logical; use a paired (sign-flip) permutation scheme.}
  \item{recycle}{logical; use circular (recycled) intervals.}
  \item{alternative}{one of \code{"two.sided"}, \code{"greater"},
    \code{"less"}.}
}
\value{
  A list with \code{T0} (observed statistics), \code{point} (pointwise
  p-values), \code{inter} (interval p-value matrix) and \code{corre}
  (corrected p-value matrix).
}
\seealso{\code{\link{fastIWT2}}, \code{\link{IWTimage}}}

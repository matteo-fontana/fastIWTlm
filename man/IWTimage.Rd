\name{IWTimage}
\alias{IWTimage}
\title{Heatmap of Interval-Wise Testing Results}
\description{
  Base-graphics heatmap of the interval-wise p-value matrix, together with
  the adjusted p-value function and the functional data, for objects of
  class \code{"IWT1"}, \code{"IWT2"} or \code{"IWTaov"}.
}
\usage{
IWTimage(IWT_result, alpha = 0.05, abscissa_range = c(0, 1),
         nlevel = 20, plot_unadjusted = FALSE)
}
\arguments{
  \item{IWT_result}{a result object of class \code{"IWT1"}, \code{"IWT2"}
    or \code{"IWTaov"}.}
  \item{alpha}{significance level used to highlight significant intervals.}
  \item{abscissa_range}{numeric length-2 vector, the range of the abscissa.}
  \item{nlevel}{number of colour levels in the heatmap.}
  \item{plot_unadjusted}{logical; also draw the unadjusted p-value function.}
}
\value{
  Called for its side effect: it produces a multi-panel plot. No value is
  returned.
}
\seealso{\code{\link{IWT2}}, \code{\link{fastIWT2}}}

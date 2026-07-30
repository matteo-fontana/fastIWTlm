\name{autoplot.IWTlm}
\alias{autoplot.IWTlm}
\alias{plot.IWTlm}
\title{ggplot2 Visualisation of IWTlm Objects}
\description{
  Produces \pkg{ggplot2} graphics for an \code{"IWTlm"} object returned by
  \code{\link{fastIWTlm}}: the functional data with the F-test significance
  bands, the functional regression coefficients with t-test significance
  bands, and the adjusted/unadjusted p-value functions. \code{plot.IWTlm}
  is a convenience alias dispatching to \code{autoplot.IWTlm}.
}
\usage{
\method{autoplot}{IWTlm}(IWTlmModel, prompt = TRUE, save_output = FALSE,
        layout_vertical = FALSE, xrange = c(0, 1), alpha1 = 0.05,
        alpha2 = 0.01, plot_adjpval = FALSE, plot_unadjpval = FALSE,
        ylim = NULL, ylimcoeff = NULL, ylab = "Functional Data",
        xlab = "abscissa_smooth", main = NULL, lwd = 1, type = "l",
        cextotal = NULL, ...)

\method{plot}{IWTlm}(x, ...)
}
\arguments{
  \item{IWTlmModel, x}{an object of class \code{"IWTlm"}.}
  \item{prompt}{logical; ask before drawing each successive plot.}
  \item{save_output}{logical; if \code{TRUE}, return the list of
    \code{ggplot} objects (for use with \code{\link{arrange_IWTlm}})
    instead of only drawing them.}
  \item{layout_vertical}{reserved layout flag.}
  \item{xrange}{numeric length-2 vector, the plotted abscissa range.}
  \item{alpha1, alpha2}{the two significance levels used for the
    highlighting bands.}
  \item{plot_adjpval, plot_unadjpval}{logical flags controlling which
    p-value functions are drawn.}
  \item{ylim, ylimcoeff}{y-axis ranges for the data and coefficient plots.}
  \item{ylab, xlab, main}{axis and title labels.}
  \item{lwd, type}{line width and type for the functional-data plot.}
  \item{cextotal}{overall character-expansion factor.}
  \item{\dots}{further arguments (passed through by \code{plot.IWTlm}).}
}
\value{
  Called for its side effect (drawing). If \code{save_output = TRUE},
  invisibly returns a list of the four \code{ggplot} objects:
  \code{data_plot}, \code{coeff_plot}, \code{pval_t_plot},
  \code{pval_plot_F}.
}
\seealso{\code{\link{fastIWTlm}}, \code{\link{arrange_IWTlm}}}

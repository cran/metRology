\name{gplot}
\alias{gplot}
\title{
Grouped plots of type "h"
}
\description{
gplot is primarily used by plot.mandel.kh to produce the underlying grouped data plot.
}
\usage{
	gplot(x, main = NULL, xlab = NULL, ylab = deparse(substitute(x)), 
		ylim = NULL, las = 1, axes = TRUE, cex.axis = 1, 
		frame.plot = axes, lwd = 1, lty = 1, col = par("col"), 
		separators = TRUE, col.sep = "lightgrey", lwd.sep = 1, 
		lty.sep = 1, zero.line = TRUE, 
		lwd.zero = 1, col.zero = 1, lty.zero = 1, 
		spacing=NA, ...)
}

\arguments{
  \item{x}{
A matrix or data frame to be plotted.
}
  \item{main}{Main title for the plot.}
  \item{xlab, ylab}{Labels for x and y axes.}
  \item{ylim}{	the y limits of the plot. 
}
  \item{las}{
	the style of the axis labels; see \code{par} for details.
}
  \item{axes}{
	a logical value indicating whether axes should be drawn
          on the plot.
}
  \item{cex.axis}{
	The magnification to be used for axis annotation
          relative to the current setting of 'cex'.
}
  \item{frame.plot}{
	Logical; If \code{TRUE} a box is drawn around the plot.
}
  \item{lwd, lty, col}{
	Graphical parameters used for the plotted vertical lines corresponding 
	to each value in x.
}
  \item{separators}{
	Logical; if \code{TRUE}, separator lines are drawn between groups of values.
}
  \item{col.sep, lwd.sep, lty.sep}{
	Graphical parameters used for the separator lines.
}
  \item{zero.line}{
	logical; if \code{TRUE} a horizontal line is drawn at zero.
}
  \item{lwd.zero, col.zero, lty.zero}{
	Graphical parameters used for the zero line.
}
  \item{\dots}{
	Other graphical parameters passed to \code{plot}.
}
  \item{spacing}{
	Spacing for data within each group, as a fraction of inter-group spacing. 
	Defaults to 0.3 or less.
}
}
\details{
\code{gplot} produces a plot of type="h", with values in x grouped by row and with 
optional vertical separators between groups. The plotting order (left to right) is 
in order of \code{stack(as.data.frame(t(x)))}; each group corresoponds to a row in x.

Because \code{gplot} is primarily a supporting function for \code{plot.mandel.kh},
it assumes a suitable object will be provided and does minimal checking
to ensure an appropriate object class. Error messages may not be 
very informative.

}
\value{
A numeric vector of mid-points of the groups along the x-axis.
}
\references{
Accuracy (trueness and precision) of measurement methods and results -- Part 2: 
Basic method for the determination of repeatability and reproducibility of a 
standard measurement method. ISO, Geneva (1994).
}
\author{S Ellison \email{s.ellison@lgcgroup.com}
}
\seealso{
\code{\link{plot.mandel.kh}}
}
\examples{
   data(RMstudy)
   h <- with(RMstudy, mandel.h(RMstudy[2:9], g=Lab))
   gplot(h, las=2) 
   	#Note the absence of indicator lines, title etc. 
   	#compared to plot(h)
}
\keyword{ hplot }

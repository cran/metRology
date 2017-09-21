\name{plot.d.ellipse}
\alias{plot.d.ellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot data ellipses
}
\description{
Plots a number of data ellipses specified by }
\usage{
	\method{plot}{d.ellipse}(x, col.ellipse = 1, lty.ellipse = 1, lwd.ellipse = 1, 
		fill = NA, density = NULL, angle = 45, add = FALSE, npoints = 100, 
		xlim = NA, ylim = NA, 
		prinax = FALSE, col.prinax = 1, lty.prinax = 1, lwd.prinax = 1, 
		xlab=NULL, ylab=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	An object of class \code{d.ellipse}
}
  \item{col.ellipse, lty.ellipse, lwd.ellipse}{
	Colour, line type and line width for the ellipse(s). Can be vectors, allowing different
	colour, line type etc. Recycled as necessary to length \code{length(x)}. 
}
  \item{fill, density, angle}{
	Fill colour, line density and line angle for each ellipse in \code{x}. See  \code{\link{polygon}} for details. 
	Can be vectors. Recycled as necessary to length \code{length(x)}. 
}
  \item{add}{
	If \code{TRUE}, ellipses are added to an existing plot. If \code{FALSE} a new plot is created. 
}
  \item{npoints}{
	Controls the number of points used to form each ellipse. See  \code{\link{data.ellipse}} for details.
}
  \item{xlim, ylim}{
	Plot limits. Ignored if \code{add == FALSE}
}
  \item{prinax}{
	If \code{TRUE} (the default), the principal axes are drawn on the plot. 
}
  \item{col.prinax, lty.prinax, lwd.prinax}{
	Colour, line type and line width for principal axes. 
}
  \item{xlab, ylab}{
	Axis labels passed to \code{plot} if \code{add == FALSE}. Defaults to dimension names
	in \code{x} or, if those are \code{NULL}, to "X" and "Y".
}
  \item{\dots}{
	Additional arguments, passed to \code{plot} if \code{add == TRUE}. 
}
}
\details{
	A series of ellipses specified in \code{x} is plotted. 
	The function is primarily used for adding ellipses to a Youden plot.
}
\value{
The function is called for its side effect, which is the drawing of ellipses.
}

\author{
	S L R Ellison
}
\seealso{
	 \code{\link{data.ellipse}}, \code{\link{youden.plot}}
}
\examples{

data(chromium)
cov.Cr <- cov.dellipse(chromium)
dellipse.Cr <- data.ellipse(cov.Cr, plot=FALSE)
plot(dellipse.Cr)

}

\keyword{ dplot }


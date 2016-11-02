\name{youden.plot}
\alias{youden.plot}
\alias{youden.plot.default}
\alias{yplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Youden plots
}
\description{
	A Youden plot is a bivariate scatter plot, named for its use by W. M Youden in 
	interlaboratory studies. This implementation includes data ellipses based on Pearson, 
	Spearman, Kendall or several robust covariance measures.
}
\usage{
youden.plot(x, ...)

yplot(x, ...)

\method{youden.plot}{default}(x, y = NULL, type = c("points", "labels", "both", "outliers"), 
	labels, probs = c(0.95, 0.99), x0, y0,  pch = par("pch"), cex = par("cex"), 
	col = par("col"), bg = par("bg"), main, xlab, ylab,  
	xlim = c("data", "ellipse", "all"), ylim = c("data", "ellipse", "all"), 
	col.axes = 2, lwd.axes = 1, lty.axes = 1, cex.lab = 0.7, pos = 3, 
	out.method = c("F", "chisq", "n"), n.out, p.out = 0.99, 
	add = FALSE,  ...)
}
\arguments{
  \item{x}{
	An R \code{numeric} object. Can be a vector (in which case \code{y} must be specified 
	and of the same length) or a two-column numeric matrix.
}
  \item{y}{
	A numeric vector of the same length as \code{x}. It is an error to provide \code{y} 
	in addition to a two-column matrix for \code{x}.
}
  \item{type}{
	The type of plot produced. See Details.
}
  \item{labels}{
	Character vector of text labels for data points. Defaults to \code{row.names(x)} if
	\code{x} is a matrix or data frame with row names, \code{names(x)} if \code{x} is a 
	named vector, and to \code{1:length(x)} or \code{1:nrow(x)} 
	as appropriate if \code{x} does not have names.
}
  \item{probs}{
	Numeric vector of probabilities for data ellipses.
}
  \item{x0, y0}{
	If specified, data ellipses will be centred on \code{(x0, y0)} instead of using the 
	location calculated from the data.  
}
  \item{pch, cex, col, bg}{
	passed to \code{points}; see \code{\link{points}} for details.
}
  \item{main, xlab, ylab}{
	Plot titles. If missing, titles are based on the names of the objects plotted. 
}
  \item{xlim, ylim}{
	Specifications for horizintal and vertical plot limits. Each can be either a length 2 
	numeric vector (as usual) or a character value matching one of \code{"data"}, 
	\code{"ellipse"} or \code{"all"}. If \code{"data"}, the relevant limits are set to 
	include the range of the data. If \code{"ellipse"}, the relevant limits are set to 
	include the whole of the outermost ellipse. If \code{"all"}, limits are set to include 
	both the data and the outermost ellipse.
}
  \item{col.axes, lwd.axes, lty.axes}{
	Colour, line width and line type for vertical and horizontal location markers
	drawn through the ellipse centre.
}
  \item{cex.lab}{
	Size for data point labels; see \code{\link{text}} for details.
}
  \item{pos}{
	a position specifier for data point labels; see \code{\link{text}} for details.
}
  \item{out.method}{
	Character specifying outlier marking. See Details.
}
  \item{n.out}{
	Number of outliers marked if \code{out.method=="n"}.
}
  \item{p.out}{
	Confidence level at which points are marked as outliers if \code{out.method} is
	one of \code{"F"} or \code{"chisq"}.
}
  \item{add}{
	If \code{TRUE}, ellipses and points are added to an existing plot.
  }
  \item{\dots}{ Named arguments passed to other functions. In particular:
  	\itemize{
	    \item{Arguments \code{cov.method}, \code{cov.control}, \code{scalefn}, and \code{locfn} 
	          will be passed to \code{cov.dellipse} to control the (optionally robust) 
	          location, scale and covariance estimates used for the ellipses.} 

	    \item{Arguments \code{sub}, \code{ann}, \code{axes}, \code{frame.plot}, and \code{asp} 
	          are passed to \code{plot} to add a sub-title, control axes and annotation and, 
	          importantly, control the plot aspect ratio.}
	          
	    \item{Arguments \code{col.ellipse}, \code{lty.ellipse}, \code{lwd.ellipse}, \code{fill}, 
                  \code{density}, \code{angle}, \code{npoints}, \code{prinax}, \code{col.prinax}, 
                  \code{lty.prinax}, and \code{lwd.prinax} are passed to \code{\link{plot.d.ellipse}}
                  to control the appearance of ellipses and ellipse principal axes.}
	}
   }
}
\details{
  \code{type} controls the type of plot produced. Allowed types and their effect are:
  \describe{
    \item{\code{points}}{Points only are drawn.}
    \item{\code{labels}}{Point labels only are drawn}
    \item{\code{both}}{Points are drawn with labels}
    \item{\code{outliers}}{Points are drawn and outlying points are labelled (see below)}
  }
  
  Ellipses are constructed based on a location and covariance matrix constructed from the data 
  by the method specified by \code{cov.method}. \code{probs} specifies the approximate coverage.
  See \code{\link{data.ellipse}} for details of covariance methods and ellipse specification.
  
  The outlier identification method, if any, is specified by \code{out.method} and controlled by 
  one of \code{n.out} or \code{p.out}. If \code{out.method} is \code{"F"} or \code{"chisq"}, 
  points with Mahalanobis distance greater than an upper critical value with probability \code{p.out} 
  are considered to be outliers. The critical values used are
  \describe{
    \item{\code{"F"}}{Mahalanobis distance greater than 
                      \code{ 2 * (n-1) * qf(p.out, 2, n-1) / (n-2)}}
    \item{\code{"chisq"}}{Mahalanobis distance greater than \code{qchisq(p.out, 2)}}
			which( md > 2 * (n-1) * qf(p.out, 2, n-1) / (n-2) ) #F dist
  }
  The Mahalanobis distance is calculated based on the covariance matrix used to consstruct plot ellipses.
  
  If \code{out.method} is \code{"n"}, the outermost \code{n.out} points (judged by Mahalanobis distance)
  are marked as outliers. Specifying \code{out.method="n"} and \code{n.out=0} suppresses outlier 
  identification.
  If outliers are marked, a list of marked outliers is included in the returned list 
  (see Value, below).

  \code{yplot} is an alias for \code{youden.plot}
}
\value{
Invisibly returns the plotted data ellipses as an object of class \code{\link{d.ellipse}}. 
}
\references{
Youden, W.J. and Steiner, E.H.  (1975) \emph{Statistical Manual of the AOAC}. AOAC 
International, Washington, US. ISBN 0-935584-15-3

ISO 13528:2005, Statistical methods for use in proficiency testing by interlaboratory
comparisons, International Organization for Standardization, Geneva (2005)

}
\author{
S L R Ellison (s.ellison@lgcgroup.com)
}

\seealso{
	\code{\link{d.ellipse}}
}
\examples{
data(chromium)
data(potassium)

( yy <- youden.plot(chromium, main="Chromium") )

#With outlier ID (F based)
youden.plot(chromium, main="Chromium", xlim='a', ylim='a', type='o', p.out=0.95)
	#Note use of xlim="a" etc. to ensure both ellipses and data are included.
	
#Top 5 most distant outliers (5 is also the default)
youden.plot(chromium, main="Chromium", xlim='a', ylim='a', type='o', out.method="n", n.out=5)

#With ellipse principal axes
#(useful to specify asp=1 or the axes will not always appear orthogonal)
youden.plot(chromium, main="Chromium", 
	xlim='a', ylim='a', type='o', p.out=0.99, prinax=TRUE, lty.prinax=2, asp=1.0)

youden.plot(potassium, main="Potassium")

#A different pairs plot ...
panel.youden <- function(x, y, ...) youden.plot(x, y, add=TRUE, type="o", cex=1, pos=1, p.out=0.95)
pairs(chromium, upper.panel=panel.youden)


}
\keyword{ hplot }

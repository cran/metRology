\name{barplot.mandel.kh}
\alias{barplot.mandel.kh}
\title{
Barplot of Mandel's h or k statistics
}
\description{
\code{barplot.mandel.kh} produces a bar plot of Mandel's statistics, suitably
grouped and with appropriate indicator lines for unusual values.
}
\usage{
	\method{barplot}{mandel.kh}(height, probs = c(0.95, 0.99), main,
		xlab = attr(height, "grouped.by"), 
		ylab = attr(height, "mandel.type"), separators = TRUE, 
		zero.line = TRUE, ylim, p.adjust = "none", 
		frame.plot = TRUE, ..., 
		col.ind = 1, lty.ind = c(2, 1), lwd.ind = 1, 
		col.sep = "lightgrey", lwd.sep = 1, lty.sep = 1, 
		lwd.zero = 1, col.zero = 1, lty.zero = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{height}{An object of class \code{'mandel.kh'} (the name is for consistency 
  		with \code{barplot}).
}
  \item{probs}{
	Indicator lines are drawn for these probabilities. Note that
	\code{probs} is interpreted as specifying two-tailed probabilities
	for Mandel's h and one-sided (upper tail) probabilities for
	Mandel's k.
}
  \item{main}{
 	a main title for the plot. If missing, the default is
 	\code{ paste(deparse(substitute(x)), " - Mandel's", 
                     attr(x, "mandel.type"), 
                     if(attr(x, "mandel.method") == "robust") "(Robust variant)")}
  }
  \item{xlab}{
  	a label for the x axis; defaults to the \code{grouped.by} attribute for
  	\code{x}.
}
  \item{ylab}{
  	a label for the x axis; defaults to the \code{mandel.type} attribute for
  	\code{x}.
}
  \item{separators}{
	Logical; if \code{TRUE}, separator lines are drawn between groups of values.
}
  \item{zero.line}{
	logical; if \code{TRUE} a horizontal line is drawn at zero.
}
  \item{ylim}{
	the y limits of the plot. For Mandel's k, the default lower limit is zero.
}
  \item{p.adjust}{
	Correction method for probabilities. If not \code{"none"}, passed to \code{p.adjust}
	prior to calculating indicator lines. Usually, indicator lines are drawn without 
	correction (that is, with \code{p.adjust="none"}); specifying a p-value correction 
	effectively turns the Mandel's 	statistics into single outlier tests.
}
  \item{frame.plot}{
	Logical; If \code{TRUE} a box is drawn around the plot.
}
  \item{\dots}{
	Other (usually graphical) parameters passed to \code{barplot}. Note that some 
	parameters appear after \dots to prevent spurious argument matching inside
	\code{barplot.default}.
}
  \item{col.ind, lty.ind, lwd.ind}{
	Graphical parameters used for the indicator lines, recyckled to \code{length(probs)}. 
	For \code{attr(x, "mandel.type")=="h"} the graphical parameters are applied to negative
	as well as positive indicator lines, applied outwards from zero.
}
  \item{col.sep, lwd.sep, lty.sep}{
	Graphical parameters used for the separator lines.
}
  \item{lwd.zero, col.zero, lty.zero}{
	Graphical parameters used for the zero line.
}
}
\details{
Mandel's statistics are traditionally plotted for inter-laboratory study data, 
grouped by laboratory, to give a rapid graphical view of laboratory bias and 
relative precision. This plot produces a grouped, side-by-side bar plot.

For classical Mandel statistics, indicator lines are drawn based on \code{qmandelh}
or  \code{qmandelk} as appropriate. For robust variants, indicator lines use
\code{qnorm} for the \eqn{h} statistic and \code{qf(probs, n, Inf)} for 
the \eqn{k} statistic. Note that this corresponds to taking the robust estimates of 
location and scale as true values, so will be somewhat anticonservative.

}
\value{
barplot.mandel.kh returns a numeric vector of mid-points of the groups 
along the x-axis.
}
\references{
Accuracy (trueness and precision) of measurement methods and results -- Part 2: 
Basic method for the determination of repeatability and reproducibility of a 
standard measurement method. ISO, Geneva (1994).
}
\author{
S Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
\code{\link{mandel.h}}, \code{\link{mandel.k}}, \code{\link{mandel.kh}},
\code{\link{pmandelh}}, \code{\link{pmandelk}}  for probabilities, quantiles etc.

See \code{\link{plot.mandel.kh}} for the 'classic' Mandel plot.
}
\examples{
   data(RMstudy)

   h <- with(RMstudy, mandel.h(RMstudy[2:9], g=Lab))
   barplot(h, las=2) # Lab 4 shows consistent low bias; 
                  # Lab 23 several extreme values.

   #Use colours to identify particular measurands:
   barplot(h, las=2, col=1:8)
   legend("bottomleft", legend=names(h), fill=1:8, cex=0.7, bg="white")
   
   #Example of Mandel's k:
   k <- with(RMstudy, mandel.k(RMstudy[2:9], g=Lab))
   barplot(k, las=2) # Lab 8 looks unusually variable; 
                  # Lab 14 unusually precise
}

\keyword{ hplot }

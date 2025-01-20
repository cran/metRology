% file RKV/man/cplot.Rd
% (c) 2009 S Ellison
\name{ cplot }
\Rdversion{1.1}
\alias{ cplot }

\title{ Consistency plot for Key Comparisons }
\description{
  Produces a consistency plot for typical metrology comparison data.
}
\usage{
  cplot(x,u,labels=names(x), p.adjust.method="holm", ordered=TRUE,
        breaks=c(0,0.001,0.01, 0.05, 0.1,1), 
        col=terrain.colors(length(breaks)-1), log.p=FALSE,
        main=paste("Consistency map -", deparse(substitute(x))), 
        subtitle=NULL, key=FALSE, 
        key.width=2.54, key.height=0.6,\dots)
}
\arguments{
  \item{x}{Vector of reported values}
  \item{u}{Vector of length \code{length(x)} of standard uncertainties}
  \item{labels}{Vector of of length \code{length(x)} labels for x-axis marks.}
  \item{p.adjust.method}{p-value adjustment method; passed to \code{p.adjust}.}
  \item{ordered}{If \code{TRUE} (the default) observations are arranged in ascending order of \code{x}.}
  \item{breaks}{Vector of breaks; passed to \code{image}.}
  \item{col}{Vector of colours of length \code{length(breaks)-1}. Passed to \code{image}.}
  \item{log.p}{If \code{TRUE}, the plot shows \code{-log_10(p)}.}
  \item{main, subtitle}{Main and subtitle for plot.}
  \item{key}{If \code{TRUE} a key is added to the plot.}
  \item{key.width, key.height}{Width and height of key, if plotted. See details for specification.}
  \item{\dots}{Graphical parameters passed to \code{image}. If present, \code{cex.axis} is passed 
               to \code{axis} for main figure axes only.}
}
\details{
	
	Calculates the (square, symmetric matrix of) optionally adjusted p-values for a 
	two-tailed z-test of \code{|x[i]-x[j]|/sqrt(u[i]^2+u[j]^2)} against zero
	and plots the p-values as an image.
	
	\code{p.adjust} is called prior to plotting to correct for multiple comparisons. 
	To suppress adjustment, set \code{p.adjust.method="none"}. 
	
	key.height is a fraction of the figure region height. key.width is the width 
	of the key area in cm, unless under 1, in which case it is interpreted as a 
	fraction of the plot region width.
	
	If \code{log.p} is \code{TRUE} and \code{subtitle} \code{NULL}, 
	a subtitle indicating the use of log.p is added to the plot, 

}
\value{
	Invisibly returns a matrix of pairwise test p-values or, if \code{log.p==TRUE},
	matrix of \code{-log_10(p)}.	
}
\seealso{
  \code{\link{p.adjust}}, \code{\link[graphics]{image}}.
}
\author{S Ellison \email{s.ellison@lgcgroup.com}}
\examples{

  data(Pb)
  cplot(Pb$value, Pb$u, key=TRUE)

}

\keyword{ hplot }

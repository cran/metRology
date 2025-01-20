\name{xs.plot}
\alias{xs.plot}
\alias{xs.plot.default}
\alias{XSplot}

\title{
Mean/Standard deviation plot with confidence region. 
}
\description{
Produces a plot of standard deviations versus means with a confidence region
based on either robust or classical estimates of location and scale.
}
\usage{

xs.plot(x, \dots)

\method{xs.plot}{default}(x,g,s, degfree, labels.arg=NA, mu, sigma, 
        probs=c(0.5, 0.95, 0.99), basis=c("robust","classical"), 
        method=c("chisq","density"), main=paste("X-S plot -", basis, "basis"), 
        xlab=deparse(substitute(x)), ylab, 
        contours=TRUE, col.contours="lightgrey", lty.contours=par("lty"), 
        lwd.contours=par("lwd"),
        label.contours=contours, format.clab="p=\%3.2f",
        pos.clab="bottomright", col.clab=col.contours, cex.clab=0.7,
        cex.label=0.7, pos=3, adj=NULL, 
        pch=par("pch"), col=par("col"), bg=par("bg"), cex=par("cex"), 
        add=FALSE, \dots)

XSplot(x, \dots)
}

\arguments{
  \item{x}{Numeric vector of values to be plotted.
}
  \item{g}{grouping factor of length \code{length(x)}.
}
  \item{s}{numeric vector of standard deviations of length \code{length(x)}
  	associated with x.
}
  \item{degfree}{A single value for degrees of freedom associated with \emph{all} 
  	the standard deviations \code{s}. Ignored if \code{s} is not supplied.
}
  \item{labels.arg}{An optional vector of point labels, 
  	coerced to character on use.
}
  \item{mu}{A single location used to centre the confidence region. The default is 
  	specified by \code{"basis"}; see Details.
}
  \item{sigma}{A measure of dispersion against which deviations x-mu can be compared.
}
  \item{probs}{A vector of probabilities for confidence region contours. 
}
  \item{basis}{Controls the nature of the location and scale estimators used to 
  	produce the confidence contours drawn on the plot. See Details for 
  	specification.
}
  \item{method}{The method used to calculate the confidence region. See Details.
}
  \item{main}{Main title for the plot.
}
  \item{xlab, ylab}{x- and y-axis labels, 
}
  \item{contours}{logical, specifying whether confidence contours should be drawn. 
}
  \item{col.contours, lty.contours, lwd.contours}{Colour, line type and line width 
  	for contour lines.
}
  \item{label.contours}{Logical, controlling whether contour lines are 
  	labelled with approximate probabilities.
}
  \item{format.clab}{Format string for contour labels, passed to \code{sprintf}.
}
  \item{pos.clab}{Specification for location of contour labels.  A vector can be 
  		provided to give multiple labels. See Details for further description.
}
  \item{col.clab, cex.clab}{Colour and expansion for contour labels.
}
  \item{pch, col, bg, cex}{Graphical parameters passed to \code{points}.
}
  \item{cex.label}{Expansion factor for point labels, passed to \code{text}.
}
  \item{pos, adj}{Specifies position/adjustment of point labels. Passed to \code{text}; 
  	see \code{text} for details.
}
  \item{add}{If TRUE the plot region is not cleared before plotting; points and contours 
  	are added to the present plot. Use \code{pch=NA} to suppress symbols if only added 
  	contours are required.
}
  \item{\dots}{Other parameters passed to \code{plot}.
}
}
\details{
A plot of standard deviations against locations is produced, together 
with optional confidence region(s) calculated (by default) by a 
method suggested in ISO 13528:2005.

If \code{s} is supplied, \code{x} is taken as a vector of locations and 
\code{s} a vector of standard deviations. \code{degfree} must be supplied
in this case.

If \code{g} is supplied and \code{s} is not, the locations and standard deviations
used are the means and standard deviations for each group. \code{degfree} is calculated 
from the median group size. Groups should, of course, be of the same size for accurate 
inference; however, using the median group size allows for some groups with missing values. 

If \code{s} and \code{g} are both supplied, \code{g} is ignored with a warning

If requested by \code{contours=TRUE}, confidence regions are drawn for each value
of \code{probs}. Contour location and shape are controlled by \code{basis}
which specifies the location and scale estimators used, and \code{method}, which 
specifies the method of calculation for the contours. Two methods are supported; 
one using the chi-squared distribution (\code{method="chisq"}) and one based on 
equal density countours (\code{method="density"}). The default, and the method 
recommended by the cited Standard, is \code{method="chisq"} and \code{basis="robust"}.

Both calculations for confidence regions require estimation of a location \eqn{\hat{\mu}}{X} 
and an estimate \eqn{\hat{\sigma}}{S} of the pooled within-group standard deviation or pooled 
estimate from \code{s}. If \code{basis="robust"}, \eqn{\hat{\mu}}{X} and \eqn{\hat{\sigma}}{S}
are calculated using \code{algA} and \code{algS} respectively. If \code{basis="classical"}, 
\eqn{\hat{\mu}}{X} and \eqn{\hat{\sigma}}{S} are the mean of the group means and the classical pooled 
standard deviation respectively. If \code{mu} or \code{sigma} are given, these are used
in place of the calculated  \eqn{\hat{\mu}}{X} and \eqn{\hat{\sigma}}{S} respectively.

If \code{method="chisq"}, contours for probability \eqn{p} are calculated as 

\deqn{s=\hat{\sigma}\exp\left ( \pm\frac{1}{\sqrt{2(n-1)}}\sqrt{\chi_{2,p}^2-n \left ( \frac{x-\hat{\mu}}{\hat{\sigma}}\right ) ^2}\right )}{s=S * exp( +-( 1/sqrt(2(n-1)) ) * sqrt(pchisq(p,2) - n * ( (x-X)/S ) ^2))}

for \eqn{x}{x} from %
\eqn{ \hat{\mu}-\hat{\sigma}\sqrt{\frac{\chi_{2,p}^2}{n}}}{X - S * sqrt(pchisq(p,2)/n)} to 
\eqn{\hat{\mu}+\hat{\sigma}\sqrt{\frac{\chi_{2,p}^2}{n}}}{X + S * sqrt(pchisq(p,2)/n)}.  


If \code{method="density"}, contours for probability \eqn{p} are calculated using 
Helmert's distribution to provide constant likelihood contours round the chosen mean 
and standard deviation. In the present implementation, these are found using
\code{uniroot} to find the mean corresponding to the required density at
given standard deviations. The density chosen is \eqn{d_{max}(1-p)} where \eqn{p} is
the probability and \eqn{d_{max}} the maximum density for Helmert's 
distribution for the requisite nunber of degrees of freedom. (See Kruskal 
(1946) for a description of Helmert's distribution and, for example, 
Pawitan (2001) for the rationale behind the choice of density
contour level.) This seems to give reasonably good results for \eqn{n \ge 3 }
but is anticonservative (particularly to high \eqn{s}) for \eqn{n = 2 }. 
 
Contours are by default labelled. Label locations can be specified using \code{pos.clab}. 
Options are \code{"top"}, \code{"topright"}, \code{"right"}, \code{"bottomright"}, 
\code{"bottom"}, \code{"bottomleft"}, \code{"left"} and \code{"topleft"}. A vector can be 
specified to give labels at more than one such location. 
Contour labels are usually placed approximately at the location(s) indicated and adjusted outward 
appropriately. For the special case of \code{method="density"} and  \code{degfree=1} (or where group 
sizes \eqn{n=1}{n=1}), for which the region is a maximu width at s=0, 
\code{"bottomright"} and \code{"bottomleft"} place labels immediately below the countour boundary 
at \eqn{s=0}{s=0} and, if specified, \code{"bottom"} is replaced with \code{c("bottomright", "bottomleft")}.

  		
\code{XSplot} is an alias for \code{xs.plot}.

}
\value{
A list with components:
\describe{
	\item{x, y}{respectively, the plotted locations and standard deviations. 
	(the names allow a simple call to \code{plot()}) 
	}
	\item{mu}{ The location and pooled SD estimates \eqn{\hat{\mu}}{X} and 
 	\eqn{\hat{\sigma}}{S} used to construct the confidence ellipsoids. 
	}
	\item{clist}{ A list of sets of coordinates for each confidence region. 
	}
}
}
\references{
ISO 13528:2005, Statistical methods for use in proficiency testing by interlaboratory
comparisons, International Organization for Standardization, Geneva (2005)

Y Pawitan, (2001) In all likelihood: Statistical Modelling and Inference Using 
Likelihood,Clarendon Press, Oxford, pp258-9

W Kruskal, American Mathematical Monthly 53, 435-438, (1946)

}
\author{S Ellison \email{s.ellison@lgcgroup.com}
}

\seealso{
\code{\link{algA}}, \code{\link{algS}}

\code{\link{axis}} for axis control; \code{\link{points}},  \code{\link{text}} for 
plotting parameters; \code{\link{sprintf}} for contour label formatting.

\code{\link{duewer.plot}} for an alternative plot for locations and associated 
standard errors or standard uncertainties; 

}
\examples{
	require(metRology)
	set.seed(1017)
	x <- rnorm(80)
	g <- gl(20,4)
	
	xs.plot(x,g)

	#Identical plot with precalculated s:
	X <- tapply(x,g,mean)
	S <- tapply(x,g,sd)
	xs.plot(X, s=S, degfree=3)

	#Specify different location and within-group SD estimates:
	xs.plot(X, s=S, degfree=3, mu=median(X), sigma=median(S))
	
	#Illustrate multiple contour labelling, point labels and further embellishment
	rv <- xs.plot(x,g, pos.clab=c("bottomleft", "bottomright"), labels=TRUE)
	abline(v=rv$mu, h=rv$s, col=2)
	
	
}

\keyword{ hplot }

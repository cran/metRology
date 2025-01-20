% file man/bootPDchisq.Rd
% (c) 2018 S Ellison
\name{ PDchisq-class }
\Rdversion{1.1}
\alias{ PDchisq-class }
\alias{ print.PDchisq}
\alias{ plot.PDchisq}
\alias{ barplot.PDchisq}
\alias{ bootPDchisq}

\title{ Methods for the object returned by \code{pdchisq}.}
\description{
  Print and plotting methods for the \code{PDchisq} object class returned by \code{\link{pdchisq}}.
}
\usage{
       \method{print}{PDchisq}(x, ...)

       \method{plot}{PDchisq}(x, type="h", ylab="Pair-difference chi-squared", ylim=NULL, ...)

       \method{barplot}{PDchisq}(height, ylab="Pair-difference chi-squared", names.arg=names(height), 
	crit.vals=TRUE, lty.crit=c(2,1), col.crit=2, lwd.crit=c(1,2), 
	probs=c(0.95, 0.99), n=length(height), ylim=NULL, ... )

	bootPDchisq(x, B=3000, probs=c(0.95, 0.99), 
		method=c("rnorm", "lhs"), keep=FALSE, labels=names(x), ...)
}

\arguments{
  \item{x, height}{
            Object of class \code{"PDchisq"}. 
	}
  \item{type}{
            The plot type. See \code{\link{plot.default}}. 
	}
  \item{ylab}{
	Label for vertical axis, passed to \code{barplot}
	}
  \item{names.arg}{
	Labels for individual bars in bar plot, passed to \code{barplot}. If \code{names(height)} 
	is \code{NULL}, bars are numbered. 
	}
  \item{crit.vals}{
	If \code{TRUE}, horizontal lines at critical values are added to the plot. 
	These are calculated from \code{link{qchisq}} based on supplied values of 
	\code{probs} and \code{n}. 
	}
  \item{lty.crit, col.crit, lwd.crit}{
	Vectors of line style parameters for plotted critical values, passed to 
	\code{\link{segments}}. Recycled to the length of \code{critical.values}
	in the supplied \code{bootPDchisq} object.}
   \item{probs}{vector of probabilities at which critical values are drawn (or, for the 
   	bootstrap, calculated and retained.}
   \item{n}{integer number of observations (NOT degrees of freedom) for critical value 
   	calculation; used (as \eqn{n-1}) by \code{\link{qchisq}}. See Details.}
   \item{B}{integer number of bootstrap replicates; passed to \code{boot.mtr.pairwise}.}
   \item{ylim}{Limits for y-axis. the default makes sure the axis begins at zero and includes 
        all values}
   \item{method}{simulation method; passed to \code{boot.mtr.pairwise}.}
   \item{keep}{logical, indicating whether to retain simulated values. 
   	Passed to \code{boot.mtr.pairwise}.}
   \item{labels}{labels for returned object; passed to \code{boot.mtr.pairwise}.}
   \item{\dots}{Parameters passed to other methods. }
}

\details{
	See \code{\link{pdchisq}} for the \dQuote{PDchisq} object description.

	The quantiles plotted by the barplot method are based on \code{qchisq}, 
	divided by \eqn{n-1}.
	Note that this assumes independence and is at best a guide; for more accurate
	quantiles, see \code{\link{bootPDchisq}}
	For the barplot method, the default critical values are \sQuote{single-observation}
	quantiles. For use as an outlier test, use probabilities adjusted for multiple comparison;
	for example, consider raising the default \code{probs} to the power \eqn{1/n}.
	
	Specifying \code{n} directly is principally intended for when the plotted values are not
   	the whole set of pair-difference chi-squared values calculated from a given data set. 
   	However, it can also be useful for limiting cases. Where a 
   	single value has very small uncertainty \eqn{s}, its pair-difference chi-squared
   	is distributed approximately as \eqn{\chi^2(n-1)}{chi-squared(n-1)}. At the other 
   	extreme, if a single value has very large uncertainty \eqn{s} compared to others, 
   	its pair-difference chi-squared is distributed approximately as \eqn{\chi^2(1)}{chi-squared(1)}.
   	
   	\code{bootPDchisq} generates a parametric bootstrap (Monte carlo simulation) 
   	      from an object of class \code{"PDchisq"}. This provides a case-specific test 
   	      of the null that all mean values are equal, with the uncertainties and/or
   	      covariances accurate. \code{bootPDchisq} is a wrapper for 
   	      \code{\link{boot.mtr.pairwise}}, which provides further information. 
}
\value{
	The \code{print} method returns the object, invisibly.
	
	The \code{plot} method returns NULL, invisibly.
	
	The \code{barplot} method returns the values at the midpoint of each bar.
	
	\code{bootPDchisq} returns an object of class \code{"bootMtrPairs"}; see
	\code{\link{bootMtrPairs-class}} for details.
}

\seealso{
  \code{\link{pdchisq}}, \code{\link{qchisq}}, \code{\link{bootMtrPairs-class}}.
}
\author{
S. L. R. Ellison \email{s.ellison@lgcgroup.com}
}
\examples{

  data(Pb)
  pdchisq.Pb<-pdchisq(Pb$value, Pb$u)  # Uses individual standard uncertainties
  names(pdchisq.Pb) <- as.character(Pb$lab)
  
  plot(pdchisq.Pb)

  barplot(pdchisq.Pb)

}
\keyword{ univar }

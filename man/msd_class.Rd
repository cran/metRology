% file man/bootMSD.Rd
% (c) 2017 S Ellison
\name{ MSD-class }
\Rdversion{1.1}
\alias{ MSD-class }
\alias{ print.MSD}
\alias{ plot.MSD}
\alias{ barplot.MSD}

\title{ Methods for the object returned by \code{msd}.}
\description{
  Print and plotting methods for the \code{MSD} object class returned by \code{\link{msd}}.
}
\usage{
       \method{print}{MSD}(x, ...)

       \method{plot}{MSD}(x, type="h", ylab="MSD", ylim=NULL, ...)

       \method{barplot}{MSD}(height, ylab="MSD", names.arg=names(height), 
	crit.vals=TRUE, lty.crit=c(2,1), col.crit=2, lwd.crit=c(1,2), 
	probs=c(0.95, 0.99), n=length(height), ylim=NULL, ... )

}

\arguments{
  \item{x, height}{
            Object of class \code{"MSD"}. 
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
	These are calculated by \code{link{qmsd}} based on supplied values of 
	\code{probs} and \code{n}. 
	}
  \item{lty.crit, col.crit, lwd.crit}{
	Vectors of line style parameters for plotted critical values, passed to 
	\code{\link{abline}}. }
   \item{probs}{vector of probabilities at which critical values are drawn.}
   \item{n}{integer number of observations for critical value calculation; passed to 
            \code{\link{qmsd}}.}
   \item{ylim}{Limits for y-axis. the default makes sure the axis begins at zero and includes 
        all values}
   \item{\dots}{Parameters passed to other methods. }
}

\details{
	See \code{\link{msd}} for the object description.

	For the barplot method, critical values are \sQuote{single-observation} quantiles. 
	For use as an outlier test, use probabilities adjusted for multiple comparison;
	for example, for the barplot method, consider raising the default \code{probs} 
	to the power \eqn{1/n}.
}
\value{
	The \code{print} method returns the object, invisibly.
	
	The \code{plot} method returns NULL, invisibly.
	
	The \code{barplot} methods return the values at the midpoint of each bar.
	
}

\seealso{
  \code{\link{msd}}, \code{\link{qmsd}}.
}
\author{
S. L. R. Ellison \email{s.ellison@lgcgroup.com}
}
\examples{

  data(Pb)
  msd.Pb<-msd(Pb$value, Pb$u)  # Uses individual standard uncertainties
  names(msd.Pb) <- as.character(Pb$lab)
  
  plot(msd.Pb)

  barplot(msd.Pb)

}
\keyword{ univar }

\name{algA}
\alias{algA}
\title{
Robust estimation of location and scale using Algorithm A
}
\description{
Algorithm A is an implementation of Huber's location and scale estimate
with iterated scale.
}
\usage{
	algA(x, k = 1.5, na.rm = FALSE, tol = .Machine$double.eps^0.25,
		maxiter = 25, verbose = FALSE)
}
\arguments{
  \item{x}{numeric vector or array of values.
  }
  \item{k}{Tuning factor; Winsorisation occurs ar k standard deviations.
}
  \item{na.rm}{a logical value indicating whether \code{NA} values should be
          removed before the computation proceeds.
}
  \item{tol}{Convergence tolerance Iteration continues until the relative 
  		change in estimated sd drops below \code{tol}.
}
  \item{maxiter}{Maximum number of iterations permitted.
}
  \item{verbose}{Controls information displayed during iteration; see Details.
}
}
\details{
Algorithm A is the robust estimate of location described in ISO 5725-5:1998. 
It proceeds by winsorisation and re-estimation of scale and location.

The argument \code{k} controls the point at which values are Winsorised 
and hence controls the efficiency. At \code{k=1.5}, the value chosen by 
ISO 5725, the estimator has asymptotic efficiency at the Normal of 0.964. 
With iterated estimate of scale and \code{k=1.5}, the estimator has a 
breakdown point of about 30% (Maronna et al, 2006, p61).

The convergence criterion for Algorithm A is not specified in ISO 5725-5:1998.
The criterion chosen here is reasonably stringent but the results will differ 
from those obtained using other choices. Use \code{verbose=2} to check the 
effect of different tolerance or maximum iteration count.

If \code{verbose} is non-zero, the current iteration number
and estimate are printed; if \code{verbose>1}, the current set
of truncated values \eqn{w}{w} is also printed.

}
\value{
  \item{mu}{Robust estimate of location}
  \item{s}{Robust estimate of scale}
}
\references{
ISO 5725-5:1998 Accuracy (trueness and precision) of measurement methods and 
results - Part 5: Alternative methods for the determination of the precision 
of a standard measurement method

Maronna R A, Martin R D, Yohai V J (2006) Robust statistics - theory and methods.
Jhn Wiley and Sons, West Sussex, England.
}
\author{
S L R Ellison \email{s.ellison@lgc.co.uk}
}
\note{
Algorithm A is identical to Huber's estimate with variable scale. 
The implementation here differs from \code{hubers} from MASS in:
   \itemize{
	\item hubers allows prior specification of fixed scale (which 
	provides higher breakdown if chosen properly) or location
	\item the option of verbose output in \code{algA}, 
	\item a maximum iteration option in \code{algA}  
	\item the convergence criterion; hubers converges on changes in \code{mu}, 
	whilst this implementation of Algorithm A converges on changes in \code{s}. 
	\item Internally, Algorithm A multiplies by a correction factor for
	standard deviation whilst \code{hubers} divides by a correction factor
	applied to the variance; the actual correction to \code{s} is identical.
    }
   The principal reasons for providing an implementation in the metRology 
   package are i) to ensure a close implementation of the cited Standard irrespective
   of other package developments (though the MASS implementation has proved very stable)
   and ii) to make the implementation easy to recognise for users of the ISO standard.
}

\section{Warning }{
Algorithm A uses the corrected median absolute deviation as the initial
estimate of scale; an error is returned if the resulting scale estimate is 
zero, which can occur with over 50\% of the data set equal. \code{huberM} in 
the robustbase package uses an alternative scale estimate in these 
circumstances.
}

\seealso{
\code{\link{algS}}, \code{\link[MASS]{hubers}}, \code{\link[robustbase]{huberM}}

}
\examples{
#Creosote example from ISO 5725-5:1998
#Means for each group are:
cm <-c(24.140, 20.155, 19.500, 20.300, 20.705, 17.570, 20.100, 20.940, 21.185)

algA(cm, verbose=TRUE)
	#Iteration 4 corresponds very closely to the ISO 5725 answer
	
}
\keyword{ robust }

\name{M-estimators}
\Rdversion{1.1}
\alias{M-estimators}
\alias{huber.estimate}
\alias{huber.estimate.default}
\alias{MM.estimate}
\alias{MM.estimate.default}
\title{
M- and MM-estimators for location.
}
\description{
Functions for calculating M- and MM-estimators for location given values 
and associated standard errors or standard uncertainties.
}
\usage{
MM.estimate(x, \dots)

\method{MM.estimate}{default}(x, u, c = 4.685, \dots)

huber.estimate(x, \dots)

\method{huber.estimate}{default}(x, u, k= 1.345, \dots)
}
\arguments{

  \item{x}{numeric vector of mean values for groups}
  \item{u}{numeric vector of standard deviations or standard uncertainties associated with the values \code{x}}
  \item{c, k}{Tuning parameters passed to other functions (see \code{rlm})}
  \item{\dots}{Parameters passed to other functions.}

}
\details{
These functions are wrappers for robust estimation using \code{rlm}. All simply 
call \code{rlm} with the formula \code{x~1} and weights \code{1/u^2}. 
}
\value{
An object of class \sQuote{loc.est}.
}
\references{
None, yet.
}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
\code{\link[MASS]{rlm}}, \code{\link{loc.est-class}}
}
\examples{


## Cd heat of vapourisation example 
## from Paule, R. C.  and Mandel, J. (1982) - see ?mpaule
x2<-c(27.044, 26.022, 26.340, 26.787, 26.796)
v<-c(3, 76, 464, 3, 14)*1e-3

MM.estimate(x2, sqrt(v))

huber.estimate(x2, sqrt(v))
}

\keyword{ univar }

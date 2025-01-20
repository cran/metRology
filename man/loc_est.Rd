\name{loc.est-class}
\Rdversion{1.1}
\alias{loc.est-class}
\alias{print.loc.est}
\title{The location estimate class}
\description{
The location estimate class contains output from a variety of estimators used in 
the metRology package.

A print method is provided.
}
\usage{
\method{print}{loc.est}(x, ...)
}

\arguments{
  \item{x}{An object of class \sQuote{loc.est}}
  \item{\dots}{Parameters passed to other functions. Currently unused.}
}
\details{
An object of class \sQuote{loc.est} is a list containing 
  \describe{
    \item{x}{Scalar estimate of location}
    \item{u}{Standard uncertainty (usually equivalent to standard error) of the location estimate.}
    \item{df}{Degrees of freedom associated with the location estimate (may be NA)}
    \item{xi}{Numeric vector of individual values contributing to the estimate}
    \item{ui}{Numeric vector of uncertainties initially associated with \code{xi}.}
    \item{dfi}{Numeric vector of degrees of freedom associated with \code{ui}.}
    \item{u.eff}{Numeric vector of \sQuote{effective uncertainties} in \code{xi} after 
    	any additional terms or adjustments are added (see below).}
    \item{w}{Numeric vector of weights associated with \code{xi} (see below).}
    \item{method}{Character string describing the method used to obtain the estimate.}
    \item{method.details}{An optional list of additional details provided by the particular method used.}
  }

The \sQuote{effective uncertainties} \code{u.eff} arise from some estimation methods (for example, 
Mandel-Paule). These typically involve either the estimation of an additional variance term, 
a scale adjustment to the output value uncertainty or (for example in the case of the 
arithmetic mean) replacement of the initial individual uncertainties with some single 
estimate based on the dispersion of values. These adjustments are usually equivalent to 
replacing the estimator used with a weighted mean using weights \eqn{1/u_{eff}^2}{1/u.eff^2}. 

The weight vector \code{w} is \emph{not} equivalent to  \eqn{1/u_{eff}^2}{1/u.eff^2}. Rather, it 
gives the ratio of prior weights \eqn{1/u_{eff}^2}{1/u.eff^2} to posterior weights, which combine 
prior weights with some additional weighting. Posterior weights arise in particular when using 
robust estimators, and are generally 1 otherwise. The returned location estimate in such cases
can be calculated as \code{sum(w*x/(u^2))/sum(w/(u^2)))}. 
%%Note that the associated standard error will usually
%%be estimated as \code{1/sum(w/(u.eff^2))} as most such estimators take the prior weights 
%%as known to a scale factor and estimate scale in conjunction with the location estimate.

\code{method.details} is an optional list that may contain anything from a short 
summary of a scale factor or additional variance to a complete object (e.g. an rlm 
object) returned by the function used to calculate the estimate.

}

\value{
The print method is called for its side effect; no value is returned.
}
\references{
None, yet
}
\author{
S. L. R. Ellison \email{s.ellison@lgcgroup.com}
}

\seealso{
\code{\link{mpaule}}
}
\examples{

## Cd heat of vapourisation example (see ?mpaule)
x2<-c(27.044, 26.022, 26.340, 26.787, 26.796)
v<-c(3, 76, 464, 3, 14)*1e-3
mp<-mpaule(x2, sqrt(v))

print(mp)

}
\keyword{ classes }

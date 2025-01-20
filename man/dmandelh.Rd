\name{Mandel-h}
\alias{dmandelh}
\alias{pmandelh}
\alias{qmandelh}
\alias{rmandelh}

\title{
     Mandel's h statistic.

}
\description{
     Density, distribution function, quantile function and random
     generation for Mandel's h statistic, a measure of relative
     deviation from a common mean.
}
\usage{
dmandelh(x, g, log = FALSE)
pmandelh(q, g, lower.tail = TRUE, log.p = FALSE)
qmandelh(p, g, lower.tail = TRUE, log.p = FALSE)
rmandelh(B, g)
}

\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{g}{number of means for which h is calculated.}
  \item{B}{ Number of observations. If 'length(B) > 1', the length is
          taken to be the number required.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x];
  	otherwise, P[X > x].}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\details{
 Mandel's h is calculated for a particular mean value \code{y[i]} in a set of 
 mean values \code{y} as 
 \code{h[i] = ( y[i] - mean(y) )/sd(y) )}
 
 The density, probabilities and quantiles can be derived from the beta distribution:
 (1+h*sqrt(g)/(g-1))/2 is distributed as Beta((g-2)/2, (g-2)/2).  
 
}
\value{
dmandelh returns the density at \code{x}, pmandelh the cumulative probability,  
qmandelh the quantiles for probability \code{p} and rmandelh returns \code{B} 
random values drawn from the distribution. 

Vector values of x, p, q and g are permitted, in which case the functions return vectors.
}

\references{
None.
}
\author{
S. L. R. Ellison, \email{s.ellison@lgcgroup.com}
}

\section{Warning }{
Note that \code{rmandelh} uses \code{B} and not \code{n} (as do most R random 
number functions) for number of random draws; this is for compatibility with 
the relevant functions for Mandel's k, for which \code{n} is conventionally
used for the number of replicates per group. Be careful when using named parameters!
}


\seealso{
	\code{\link{pmandelk}}
}
\examples{
	#Generate the 95% and 99% quantiles for comparison with tables in 
	#ISO 5725:1996 Part 2:
	n <- 3:30
	round(qmandelh(0.975, n), 2) #95% 2-tailed

	round(qmandelh(0.995, n), 2) #99% 2-tailed
	
}

\keyword{ distribution }


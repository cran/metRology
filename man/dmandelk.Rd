\name{Mandel-k}
\alias{dmandelk}
\alias{pmandelk}
\alias{qmandelk}
\alias{rmandelk}

\title{
     Mandel's k statistic.
}
\description{
     Density, distribution function, quantile function and random
     generation for Mandel's k statistic, a measure of relative 
     precision compared to a common variance.
}
\usage{
dmandelk(x, g, n, log = FALSE)
pmandelk(q, g, n, lower.tail = TRUE, log.p = FALSE)
qmandelk(p, g, n, lower.tail = TRUE, log.p = FALSE)
rmandelk(B, g, n)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{g}{number of groups for which k is calculated.}
  \item{n}{number of observations in each group of data for which k is calculated.}
  \item{B}{ Number of observations. If 'length(B) > 1', the length is
          taken to be the number required.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x];
  	otherwise, P[X > x].}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\details{
	Mandel's k for one of a set of \eqn{g} standard deviations \eqn{s} is calculated as
	\deqn{k=\frac{s_{ij}^2}{\sum_{i=1}^p{s_{ij}^2/p}}}
	
	Since the numerator is chi-squared(n-1), or Gamma((n-1)/2, 2),  and the denominator 
	can be written as the sum of the same quantity and a pooled variance with distribution
	Gamma((g-1)*(n-1)/2, 2), k is distributed as Beta((n-1)/2, (g-1)(n-1)/2).
	Quantiles, probabilities, density and random numbers can therefore be generated
	from the Beta distribution. For example, qmandelk is calculated as 
	\code{sqrt( g * qbeta( (n-1)/2, (g-1)*(n-1)/2))}.
	
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
S. L. R. Ellison, \email{s.ellison@lgc.co.uk}
}

\section{Warning }{
Note that \code{rmandelk} uses \code{B} and not \code{n} (as do most R random 
number functions) for number of random draws; this is because \code{n} is conventionally
used for the number of replicates per group. Be careful when using named parameters!
}

\seealso{
	\code{\link{pmandelh}}
}
\examples{
	#Generate the 95% and 99% quantiles for comparison with tables in 
	#ISO 5725:1996 Part 2:
	
	round(qmandelk(0.95, g=3:30, n=3), 2) #95% upper tail

	round(qmandelk(0.99, g=3:30, n=3), 2) #99% upper tail
	
}

\keyword{ distribution }

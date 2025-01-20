\name{Scaled t distribution}
\Rdversion{1.1}
\alias{dt.scaled}
\alias{pt.scaled}
\alias{qt.scaled}
\alias{rt.scaled}
\title{
Scaled and shifted t distribution.
}
\description{
Student's t distribution for 'df' degrees of freedom, shifted by 'mean' and
scaled by 'sd'.}
\usage{
dt.scaled(x, df, mean = 0, sd = 1, ncp, log = FALSE)
pt.scaled(q, df, mean = 0, sd = 1, ncp, lower.tail = TRUE, log.p = FALSE)
qt.scaled(p, df, mean = 0, sd = 1, ncp, lower.tail = TRUE, log.p = FALSE)
rt.scaled(n, df, mean = 0, sd = 1, ncp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
          taken to be the number required.}
  \item{df}{degrees of freedom (> 0, maybe non-integer).  \code{df = Inf} is
          allowed.}
  \item{mean}{mean value for the shifted, scaled distribution.}
  \item{sd}{Scale factor for the shifted, scaled distribution.}
  \item{ncp}{non-centrality parameter delta; currently except for \code{rt()},
          only for \code{abs(ncp) <= 37.62}. If omitted, use the central t
          distribution.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x]; otherwise, P[X > x].}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\details{
These are wrappers for the corresponding t distribution functions in package \code{stats}. 

The scaled, shifted t distribution has mean \code{mean} and variance \code{sd^2 * df/(df-2)}

The scaled, shifted t distribution is used for Monte Carlo evaluation when a value x has 
been assigned a standard uncertainty u associated with with df degrees of freedom; 
the corresponding distribution function for that is then \code{t.scaled} with
\code{mean=x}, \code{sd=u} and \code{df=df}.
}
\value{
     \code{dt.scaled} gives the density, \code{pt.scaled} gives the distribution function, 
     \code{qt.scaled} gives the quantile function, and \code{rt.scaled} generates random deviates.

     Invalid arguments will result in return value \code{NaN}, with a
     warning.

}

\author{
S. L. R. Ellison \email{s.ellison@lgcgroup.com}
}

\seealso{
\code{\link[stats]{TDist}}
}
\examples{

	u<-rt.scaled(20, df=5, mean=11, sd=0.7)
	
	qt.scaled(c(0.025,0.975), Inf, mean=10, sd=1) #10 +- 1.96*sd
	
	require(graphics)
	hist(rt.scaled(10000, df=4, mean=11, sd=0.7), breaks=50, probability=TRUE)
	x<-seq(0,25, 0.05)
	lines(x,dnorm(x,mean=11, sd=0.7), col=2)

}
\keyword{ distribution }


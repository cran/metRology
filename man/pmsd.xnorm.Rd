\name{.pmsd.xnorm}
\Rdversion{1.1}
\alias{.pmsd.xnorm}

\title{
Internal functions for package RKC.
}
\description{
The following functions are used internally by pmsd and qmsd.
}
\usage{
.pmsd.xnorm(q, x, n, sd = 1, scale = FALSE)
}

\arguments{
  \item{q}{vector of quantiles}
  \item{x}{scalar value}
  \item{n}{number of observations in data set to which msd is applied}
  \item{sd}{a standard deviation}
  \item{sd}{Standard deviation, used to scale q.}
  \item{scale}{If scale is \code{TRUE}, \code{sd} is divided by \code{sqrt(2)}.}
  }
\details{
Provides the distribution function for msd given a value x, using the cumulative 
probability for an order statistics givn in Mood, Graybill and Boes (1974) pp252ff.  
\code{pmsd} integrates (using \code{integrate}) over x to obtain the distribution 
function for msd for random x.}
\value{
A vector of probabilities}
\references{
Mood, Graybill and Boes (1974) pp252ff
}
\author{S Ellison \email{s.ellison@lgc.co.uk}}

\seealso{
\code{\link{pmsd}}, \code{\link{msd}}
}
\examples{

## The function is currently defined as
function(q, x, n, sd=1, scale=FALSE) {

if(scale) sd <- sd/sqrt(2)

pxnorm<-function(q,x,sd=1) ifelse(q>0, pnorm(x+q, 0, sd)-pnorm(x-q, 0, sd), 0) 

Fy<-rep(0, length(q))

n.med<-floor(n/2) #exact for odd samples, low for even
#Note that for n values there are n-1 differences,
#so an even-n case is an odd-median case
ph<-pxnorm(q,x,sd)

for(j in n.med:(n-1)) Fy <- Fy + choose(n-1,j) * (ph^j) * (1-ph)^(n-j-1)

return(Fy)

  }
}


\keyword{ internal }

% file RKV/man/pmsd.Rd
% (c) 2009 S Ellison
\name{ pmsd }
\Rdversion{1.1}
\alias{ pmsd }
\alias{ qmsd } 

\title{ Median scaled difference probabilities and quantiles }
\description{
  Cumulative lower tail probability and quantile for median of scaled differences.
}
\usage{
	pmsd(q, n, sd=1, scale=TRUE)
	qmsd(p, n, sd=1, scale=TRUE)
}
\arguments{
  \item{q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations from which msd was calculated. }
  \item{sd}{Standard deviation, used to scale q.}
  \item{scale}{If scale is \code{TRUE}, \code{sd} is divided by \code{sqrt(2)} prior to calculating
  	\code{p} or \code{q}.}
}
\details{
	
	\code{pmsd} and \code{qmsd} implement exact (for even \code{n}) or approximate (odd \code{n}) probabilities 
	and quantiles for the median scaled difference applied to a single observation in a normal
	distribution.

	\code{n} is the number of observations in the data set of interest and \emph{not} the degrees of 
	freedom or number of differences (msd for a value x[i] in a set of \code{n} observations 
	involves \code{n-1} scaled differences).
	
	The probabilities are calculated using quadrature integration over a distribution of an order 
	statistic, and may be quite slow (seconds for a vector of several hundred 100 values of q on an 
	Intel x86 machine running at 1GHz). \code{qmsd} is obtained even more slowly by root-finding 
	from \code{pmsd} using \code{uniroot}.
	
	\code{n}, \code{p} and \code{q} are recycled to the length of the longest, as necessary.

	Note that both functions are appropriate for the distribution of single values. If 
	seeking an outlier test, adjust \code{p} for \code{n} comparisons before applying \code{qmsd}.
}
\value{
	A vector of length \code{length(p)} or \code{length(q)} (or, if longer, \code{length(n)}) of 
	quantiles or probabilities respectively.
}
\author{S Ellison \email{s.ellison@lgc.co.uk}}
\seealso{
  \code{\link{msd}} for calculation of MSD values, and \code{\link{bootMSD}} for 
  a parametric bootstrap (MCS) method of obtaining \emph{p}-values and quantiles
  for the more general non-IID case.
}
\examples{

  data(Pb)
  msd(Pb$value)          # Uses mad(Pb$value) as scale estimate
  msd(Pb$value, Pb$u)    # Scales differences using standard uncertainties

}

\keyword{ distribution }

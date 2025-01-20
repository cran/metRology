% file RKV/man/pmsd.Rd
% (c) 2009 S Ellison
\name{ pmsd }
\Rdversion{1.1}
\alias{ pmsd }
\alias{ qmsd } 
\alias{ dmsd } 

\title{ Median scaled difference probabilities and quantiles }
\description{
  Cumulative lower tail probability and quantile for median of scaled differences.
}
\usage{
	dmsd(q, n, method=c('fast', 'exact', 'even', 'asymp'), max.odd=199)
	
	pmsd(q, n, lower.tail = TRUE, 
		method=c('fast', 'exact', 'even', 'asymp'), max.odd=199)
	           
	qmsd(p, n, lower.tail = TRUE, 
		method=c('fast', 'exact', 'even', 'asymp'), max.odd=199)
}
\arguments{
  \item{q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations from which msd was calculated. Unused (and can be missing)
           for \code{method="asymp"}}
  \item{lower.tail}{logical; if TRUE (the default), probabilities are P[X <= x];
          otherwise, P[X > x]. }
  \item{method}{Calculation method. See details.}
  \item{max.odd}{Highest odd \code{n} for which exact values are calculated.}
}
\details{
	
	\code{pmsd}, \code{dmsd} and \code{qmsd} return probabilities, densities and quantiles, respectively, 
	for the median scaled difference applied to a single observation in a standard normal
	distribution, where otehr values are also IID normal.

	\code{n} is the number of observations in the data set of interest and \emph{not} the degrees of 
	freedom or number of differences (msd for a value x[i] in a set of \code{n} observations 
	involves \code{n-1} scaled differences).

	\code{n}, \code{p} and \code{q} are recycled to the length of the longest, as necessary.
	
	\code{method} determines the method of calculation. 
	For \code{method="fast"}, probabilities are calculated using monotonic spline 
	interpolation on precalculated probabilities. \code{qmsd} with \code{method="fast"} is obtained 
	by root-finding on the corresponding spline function using \code{uniroot}, and densities are 
	estimated from the first derivative of the interpolating spline. This provides fast 
	calculation, and values for most practical probabilities are within 10^-6 of exact calculations. 
	For high probabilites and for low quantiles (below 0.48) at high \code{n}, fast quantile accuracy 
	is poorer due to the very low function gradients in this regions, but is still guaranteed 
	monotonic with \code{p}.

	For \code{method="exact"}, probabilities and densities are calculated using quadrature 
	integration for an order statistic. For odd \code{n}, this requires a double integral. Values for 
	odd \code{n} accordingly take about an order of magnitude longer to obtain than for even \code{n}. 
	This can be slow (seconds for a vector of several hundred values of \code{q} on an Intel x86
	machine running at 1-2GHz). \code{qmsd} with \code{method="exact"} is obtained by root-finding from 
	\code{pmsd(..., method="excat")} using \code{uniroot}, and is over an order of magnitude slower than 
	\code{pmsd} pmsd. 
	
	For \code{method="exact"}, asymptotic (large \eqn{n}) probabilities, densities and quantiles are 
	returned. \code{n} is unused and can be missing. 
	
	For \code{method="exact"}, odd \code{n} above \code{max.odd} are replaced with the next lower 
	even value.  This  provides a fair approximation for \code{n} above 30 (though the fast method is better)
	and a good approximation above the default of 199. Values of \code{max.odd} above 199 are not recommended 
	as integration can become unstable at high odd \code{n}; a warning is issued if \code{max.odd > 199}. 
	
	For \code{method="even"}, an exact calculation is performed with any odd \code{n} replaced with the 
	next lower even value. This is equivalent to setting \code{method="exact"} and \code{max.odd=0}. 
	This is provided for interest only; the \code{method="fast"} method provides a substantially better 
	approximation for odd \code{n} than \code{method="even"} and is faster.

	Note that these functions are appropriate for the distribution of single values. If 
	seeking an outlier test in a data set of size \eqn{N}, either adjust \code{p} for \eqn{N} 
	comparisons before applying \code{qmsd} to find a critical value, or adjust the returned 
	\emph{p}-values using, for example, Holm adjustment.
}
\value{
	A vector of length \code{length(p)} or \code{length(q)} (or, if longer, \code{length(n)}) of 
	cumulative probabilities, densities or quantiles respectively.
}
\author{S Ellison \email{s.ellison@lgcgroup.com}}
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

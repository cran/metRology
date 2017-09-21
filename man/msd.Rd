% file RKV/man/msd.Rd
% (c) 2009 S Ellison
\name{ msd }
\Rdversion{1.1}
\alias{ msd }

\title{ Median scaled difference }
\description{
  Generates median of scaled differences from each point in a data set to all other points..
}
\usage{
        msd(x, s=mad , \dots)
	
}
\arguments{
  \item{x}{Vector of observations}
  \item{s}{Either a function returning an estimate of scale for \code{x} or a vector of 
  	length \code{length(x)} of standard errors or standard uncertainties in \code{x}.}
  \item{\dots}{Parameters passed to \code{s} if \code{s} is a function.}
}
\details{
	
	For each observation \code{x[i]}, \code{msd} calculates the median of \code{|x[i]-x[j]|/sqrt(s[i]^2+s[j]^2), j!=i},
	that is, the median of differences divided by the estimated uncertainties of the distance.
	
	If \code{s} is a function, it is applied to \code{x} and replicated to length \code{length(x)}; if 
	a scalar, it is replicated to length \code{length(x)}. 
	
	The median scaled difference is a measure of how \sQuote{far} an individual observation is from 
	the majority of the other values in the data set. As a rule of thumb, values above
	2 are indicative of a suspect \code{(x[i], s[i])} data pair; that is, a value \code{x[i]} that 
	is remote from a large fraction of the remaining data given its associated standard 
	error or standard uncertainty \code{s[i]}.
	
}
\value{
	An object of class "MSD", consisting of a vector of length \code{length(x)} of median 
	scaled absolute deviations for each observation, with attributes:
        \item{names}{character vector of names, taken from \code{x}}
        \item{x}{values supplied as \code{x}}
        \item{s}{values supplied as \code{s}}
        
        
        Print and plotting methods are currently provided for the \code{"MSD"} class; 
        see \code{\link{MSD-class}}.
}

\seealso{
  \code{\link{pmsd}}, \code{\link{qmsd}}, \code{\link{MSD-class}}.
}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}
\examples{

  data(Pb)
  msd(Pb$value)  # Uses mad(Pb$value) as scale estimate
  msd(Pb$value, Pb$u)  # Scales differences using standard uncertainties

}

\keyword{ univar }

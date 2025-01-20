% file metRology/man/pdchisq.Rd
% (c) 2018 S Ellison
\name{ pdchisq }
\Rdversion{1.1}
\alias{ pdchisq }

\title{ Pair-difference chi-squared statistic  }
\description{
  Generates the pair-difference chi-squared statistic for each point in a data set, 
  summing over scaled differences frm other points.
}
\usage{
      
        pdchisq(x, s=sd, cov=NULL, cor = NULL, na.rm=FALSE, \dots)
	
}
\arguments{
  \item{x}{Vector of observations}
  \item{s}{Either a function returning an estimate of scale for \code{x} or a vector of 
  	length \code{length(x)} of standard errors or standard uncertainties in \code{x}.}
   \item{cov, cor}{
		Covariance or correlation matrix, respectively, describing the covariance structure 
		across \code{x}. 
	}
  \item{na.rm}{logical.  Controls whether missing values should be removed. 
  	NOTE: \code{na.rm} is passed to \code{s} as well as \code{sum}; 
  	\code{s} must accordingly accept an argument \code{na.rm}.}
  \item{\dots}{Parameters passed to \code{s} if \code{s} is a function.}
}
\details{
	
	For each observation \eqn{x_j}{x[j]}, \code{pdchisq} calculates the pairwise chi-squared statistic
	\deqn{\chi_{j \rm PD}^2 = \frac{1}{N-1} \sum_{i=1}^{N} (x_i - x_j)^2/(s_i^2 + s_j^2 - 2 {\rm cov}(x_i, x_j)) }{%
			sum( (x[i]-x[j] )^2 / (s[i]^2 + s[j]^2 - 2* cov[i,j]) }
	that is, the sum of squared differences divided by the estimated uncertainties of the distances.
	
	If \code{s} is a function, it is applied to \code{x} and replicated to length 
	\code{length(x)}. Note that \code{na.rm} is passed to \code{s}, so 
  	\code{s} must accept an argument \code{na.rm}. For some scale 
  	functions that may require \code{s} to be defined as a wrapper to 
  	avoid unused argument warnings. For example, if \code{foo} is a scale 
  	function that does not accept an \code{na.rm} argument, use 
  	\code{s=function(x, na.rm, ...) foo(x, ...) }.
  	
	
	If \code{s} is a scalar, it is replicated to length \code{length(x)}. 
	
	If \code{cov} is present, \code{s} is silently ignored.
	If \code{cor} is present, \code{cov} is constructed from \code{cor} and \code{s}
	
	The pair-difference chi-squared statistic is a measure of how \sQuote{far} an individual
	observation is from all the other values in the data set, taking account of uncertainties.
	
}
\value{
	An object of class "PDchisq", consisting of a vector of length \code{length(x)} of median 
	scaled absolute deviations for each observation, with attributes:
        
        \item{names}{character vector of names, taken from \code{x}}
        \item{x}{values supplied as \code{x}}
        \item{s}{values supplied as \code{s}}
        
        
        Print and plotting methods are provided for the \dQuote{PDchisq} class; 
        see \code{\link{PDchisq-class}}.
}

\seealso{
  \code{\link{PDchisq-class}}, \code{\link{bootPDchisq}}.
}
\references{
	R J Douglas, A G Steele (2006) Pair-difference chi-squared statistics for Key Comparisons, 
	Metrologia 43, 89-97
}
\author{
S. L. R. Ellison \email{s.ellison@lgcgroup.com}
}
\examples{

  data(Pb)
  pdchisq(Pb$value)  # Uses mad(Pb$value) as scale estimate
  pdchisq(Pb$value, Pb$u)  # Scales differences using standard uncertainties

}

\keyword{ univar }

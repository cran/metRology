\name{REML location estimate}
\alias{reml.loc}
\alias{reml.loc.default}

\title{
Restricted maximum likelihood estimate of location
}
\description{
Calculates REML estimate of location, with standard error, assuming 
a random-effects model
}
\usage{
	reml.loc(x, \dots, na.rm = FALSE)

	\method{reml.loc}{default}(x, s, n = NULL, groups = NULL, na.rm = FALSE,
	          tol=.Machine$double.eps^0.5, REML=TRUE, \dots)
}
\arguments{
  \item{x}{numeric vector of mean values for groups, or (if \code{groups} 
  	is given) of individual observations}
  \item{s}{numeric vector of length \code{length(x)} of standard deviations or 
  	standard uncertainties associated with the values \code{x}.
}
  \item{n}{integer giving the number of observations in each group. May be a vector 
  	of length \code{length(x)}. If \code{n} is \code{NULL}, \code{s} is interpreted as 
  	a vector of standard uncertainties or standard errors. \code{n} is recycled 
  	to \code{length(x)}
}
  \item{groups}{factor, or vetor which can be coerced to factor, of groups. If 
  	present, \code{x} is interpreted as a vector of individual observations 
  	and \code{s} and \code{n} ignored, if present, with a warning. 
}
  \item{na.rm}{logical: if \code{TRUE},  \code{NA} values are removed 
  	before processing.
}
  \item{tol}{numeric tolerance for convergence, used by \code{optimize()}.
}
  \item{REML}{logical: if \code{TRUE} (the default),  the function optimises the REML 
        criterion (see Details). If FALSE, the maximum likelihood criterion is used.
}
  \item{\dots}{Further parameters passed to \code{optimize()}.
}
}
\details{
\code{reml.loc} finds an excess variance \eqn{\tau^2}{tau^2} and location \eqn{\mu}{mu} that maximise the 
restricted maximum likelihood criterion.

The estimator assumes a model of the form
\deqn{x_i=\mu+b_i+e_i}{x[i]=mu+b[i]+e[i]}
in which \eqn{b_i}{b[i]} is drawn from \eqn{N(0, \tau^2)}{N(0, tau^2)} and
\eqn{e_i}{e[i]} is drawn from \eqn{N(0, \sigma_i^2)}{N(0, sigma[i]^2)}.

By default the function maximises the data-dependent part of the negative log restricted likelihood:

		\deqn{\frac{1}{2}  \left( \sum_{i=1}^{k}\frac{(x_i-mu)^2}{u_i^2}  + \sum_{i=1}^{k}log(u_i^2) + 
			 log\left(\sum_{i=1}^{k}(1/u_i^2)\right) \right)}

where \eqn{u_i=s_i^2 + \tau^2}{u[i]=s[i]^2 + tau^2} and \eqn{k} is the number of mean values.
If \code{REML=FALSE}, the final term is omitted to give the maximum likelihood criterion.


This implementation permits input in the form of:
\itemize{
    \item means \code{x} and standard errors \code{s}, in which case neither \code{n} nor 
    	\code{groups} are supplied;
    \item means \code{x}, standard deviations \code{s} and group size(s) \code{n}, 
    	standard errors then being calculated as \code{s/sqrt(n)}
    \item individual observations \code{x} with a groupinf factor \code{groups}, in 
    	which case standard errors are calculated from the groups using \code{tapply}.
}


}
\value{
  A loc.est object; see loc.est for details. In the returned object, individual
  values \code{xi} are always input means (calculated from groups and \code{n} as 
  necessary); \code{method.details} is returned as a list containing:
  \describe{
      \item{mu}{The estimated location.
      }
      \item{s}{The standard error in the location.
      }
      \item{tau}{The excess variance (as a standard deviation).
      }
      \item{REML}{Logical, giving the value of REML used.
      }
  }

}
\references{
None, but see documentation for the metafor package for a more general implementation of REML.
}
\author{
S L R Ellison \email{s.ellison@lgc.co.uk}
}
\seealso{
\code{\link{loc.est-class}}
}
\examples{

  #PCB measurements in a sediment from Key Comparison CCQM-K25
  #s are reported standard uncertainties
  pcb105 <- data.frame(x=c(10.21, 10.9, 10.94, 10.58, 10.81, 9.62, 10.8),
               s=c(0.381, 0.250, 0.130, 0.410, 0.445, 0.196, 0.093))
               		
  with( pcb105, reml.loc(x, s) )

}
\keyword{ univar }

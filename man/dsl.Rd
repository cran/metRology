\name{derSimonian-Laird}
\alias{dsl}
\alias{dsl.default}

\title{
derSimonian-Laird estimator
}
\description{
Calculates derSimonian-Laird estimate of location, with standard error, assuming 
a random-effects model
}
\usage{
	dsl(x, \dots, na.rm = FALSE)

	\method{dsl}{default}(x, s, n = NULL, groups = NULL, \dots, na.rm = FALSE)
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
  \item{\dots}{Further parameters passed to other methods.
}
}
\details{
\code{dsl} implements the derSimonian-Laird random-effects estimate of location,
using the implementation described by Jackson (2010).

The estimator assumes a model of the form
\deqn{x_i=\mu+b_i+e_i}{x[i]=mu+b[i]+e[i]}
in which \eqn{b_i}{b[i]} is drawn from \eqn{N(0, \tau^2)}{N(0, tau^2)} and
\eqn{e_i}{e[i]} is drawn from \eqn{N(0, \sigma_i^2)}{N(0, sigma[i]^2)}.

The estimator forms a direct calculation of \eqn{\tau}{tau}, and uses this to 
form revised estimates of standard error \eqn{\sqrt{s_i^2+\tau^2}}{sqrt(s[i]^2+tau^2)} 
in \code{x}, calculates weights as the inverse of these and in turn calculates a 
weighted mean, allowing for any calculated excess variance \eqn{\tau^2}{tau^2}.

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
  }

}
\references{
Jackson et al. (2010) J Stat Plan Inf 140, 961-970 
}
\author{
S L R Ellison \email{s.ellison@lgcgroup.com}
}
\seealso{
\code{\link{loc.est-class}}
}
\examples{

  #PCB measurements in a sediment from Key Comparison CCQM-K25
  #s are reported standard uncertainties
  pcb105 <- data.frame(x=c(10.21, 10.9, 10.94, 10.58, 10.81, 9.62, 10.8),
               s=c(0.381, 0.250, 0.130, 0.410, 0.445, 0.196, 0.093))
               		
  with( pcb105, dsl(x, s) )

}
\keyword{ univar }

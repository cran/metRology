\name{Triangular}
\Rdversion{1.1}
\alias{dtri}
\alias{ptri}
\alias{qtri}
\alias{rtri}
\title{
The triangular distribution.
}
\description{
Density, distribution function, quantile function and random
generation for the triangular distribution with range 'min' to 'max' 
and mode equal to 'mode'.
}
\usage{
dtri(x, min=-sqrt(6), max=sqrt(6), mode = (min + max)/2, 
		log = FALSE)
		
ptri(q, min=-sqrt(6), max=sqrt(6), mode = (min + max)/2, 
		lower.tail = TRUE, log.p = FALSE)
		
qtri(p, min=-sqrt(6), max=sqrt(6), mode = (min + max)/2, 
		lower.tail = TRUE, log.p = FALSE)
		
rtri(n, min=-sqrt(6), max=sqrt(6), mode = (min + max)/2)
}

\arguments{
  \item{x, q}{ Vector of quantiles.}
  \item{p}{ Vector of quantiles.}
  \item{n}{ Number of observations. If 'length(n) > 1', the length is
          taken to be the number required.}
  \item{min}{Vector of lower limits of distribution.}
  \item{max}{Vector of upper limits of distribution.}
  \item{mode}{Vector of modes}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x]; otherwise, P[X > x].}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}

\details{
The triangular distribution has density 
\deqn{f(x)=2*(x-min) / ((max-min)*(mode-min))  (min < x <= mode)}
\deqn{f(x)=2*(max-x) / ((max-min)*(max-mode))  (mode < x < max)}
and 0 elsewhere.

The mean is
\deqn{\frac{1}{3}(min + mode + max)}{(min + mode + max)/3}
and the variance is
\deqn{\frac{1}{18}(min^2 + mode^2 + max^2 - min*mode - min*max - mode*max)}

The default values of \code{min}, \code{max} and \code{mode} give a distribution with 
mean 0 and unit variance.

If \code{min>max}, \code{min} amd \code{max} will be silently interchanged. If 
mode is not within \code{[min, max]}, the functions return \code{NA}, with a warning.

\code{rtri} calls \code{runif(n, 0, 1)} to generate probabilities which are passed to 
\code{qtri}. 
}
\value{
A vector of densities, probabilities, quantiles or random deviates. 
\code{dtri} gives the density, \code{ptri} gives the distribution function, 
\code{qtri} gives the quantile function, and \code{rtri} generates random deviates.

}

\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}


\seealso{
\code{\link[stats]{runif}}
}
\examples{
     require(graphics)
     x<-seq(-3,3,0.02)

     par(mfrow=c(2,1))
     plot(x, dtri(x), type="l", main="Density")
     plot(x, ptri(x), type="l", main="p(X<x)")

     u <- rtri(5000)

     var(rtri(10000,-1,1))  # ~ = 1/6 = 0.167

}
\keyword{ distribution }

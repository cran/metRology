\name{data.ellipse}
\alias{data.ellipse}
\alias{summary.d.ellipse}
\alias{d.ellipse}
\alias{print.d.ellipse}
\title{
Construct data ellipses suitable for use with Youden plots.
}
\description{
Constructs and optionally plots a set of probability ellipses for a bivariate 
normal distribution with defined centre and covariance.
}
\usage{
	data.ellipse(cov, probs = 0.95, plot = TRUE, npoints = 100, ...)

	\method{summary}{d.ellipse}(object, ...)

	\method{print}{d.ellipse}(x,  ...)
}


\arguments{
    \item{cov}{
  Covariance and location object of class \code{cov.dellipse} as returned by \code{cov.dellipse()}
  }
    \item{probs}{
  A vector of probabilities at which ellipses will be constructed.
  }
    \item{plot}{
  Logical specifying whether the ellipses constructed will additionally be plotted. 
  If \code{TRUE}, the result is plotted using \code{plot.d.ellipse()}
  }
    \item{npoints}{
  Integer number of points for each quadrant of the ellipses returned.
  }
    \item{object}{
    Object of class \code{d.ellipse} (for summary method).
  }
    \item{x}{
    Object of class \code{d.ellipse} (for print method.
  }
    \item{\dots}{
  Arguments passed to other methods, particularly \code{plot.d.ellipse} and, for \code{print},
  \code{format} and \code{print.default}.
  }
}


\details{
\code{data.ellipse} constructs and returns one set of x, y coordinates for each value of 
\code{probs}, in a form that can be passed directly to \code{polygon}.

Ellipses are constructed from the upper \code{probs} quantile of the F distribution using

\code{T = sqrt( 2 * (n - 1) * qf(probs, 2, n - 1)/(n - 2))}

where \code{n} is the number of pairs used in forming the covariance matrix.
If the number of points is missing or NA, \code{Inf} is substituted.

Summary and print methods are provided. The summary method returns a list with the 
same names as class \code{d.ellipse}, each containing a default summary of the respective 
member of the \code{d.ellipse} object. The print method returns its argument invisibly.

}
\value{
An object of class \code{d.ellipse}, consisting of:

\item{ellipses}{A named list of ellipsoids named for each probability in \code{probs}. 
Each is a \eqn{4*npoints\times 2}{4*n by 2} matrix suitable for passing direct to 
\code{polygon}.}

\item{probs}{Numeric vector of probabilities as supplied by \code{probs}}

\item{cov}{Covariance object of class  \code{cov.dellipse} as provided in \code{cov}}

}
\references{
	ISO 13528:2005 Statistical methods for use in proficiency testing by 
	interlaboratory comparisons. (2005) International organization for 
	Standardizaton, Geneva

	
	Jackson, J. E. (1956) \emph{Quality control methods for two related variables}.
	Industrial Quality Control, Vol. 7, pp. 2-6 

	Jackson, J. E. (1959) \emph{Control Methods for Several Related Variables}.
	Technometrics, Vol. 1, pp. 359-377
}
\author{
	S L R Ellison
}

\seealso{
	\code{\link{cov.dellipse}}, \code{\link{plot.d.ellipse}}, \code{\link{polygon}} 
}
\examples{
data(chromium)
cov.Cr <- cov.dellipse(chromium)
dellipse.Cr <- data.ellipse(cov.Cr, plot=FALSE)
summary(dellipse.Cr)
}
\keyword{ dplot }
\keyword{ hplot }

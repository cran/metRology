\name{welch.satterthwaite}
\Rdversion{1.1}
\alias{w.s}
\alias{welch.satterthwaite}
\title{
Welch-Satterthwaite effective degrees of freedom}
\description{
Provides the Welch-Satterthwaite effective degrees of freedom given standard 
uncertainties and associated degrees of freedom.

w.s is an alias for welch.satterthwaite. 
}
\usage{
w.s(ui, df, ci = rep(1, length(ui)), uc=sqrt(sum((ci*ui)^2)))

welch.satterthwaite(ui, df, ci = rep(1, length(ui)), 
                       uc=sqrt(sum((ci*ui)^2)))
}

\arguments{
  \item{ui}{Standard uncertainties}
  \item{df}{Degrees of freedom}
  \item{ci}{Sensitivity coefficients \eqn{dy/dx_i}}
  \item{uc}{Combined standard uncertainty}
}
\details{
Implements the Welch-Satterthwaite equation as provided in the ISO Guide to the expression of
uncertainty in measurement (1995) (See JCGM 100:2008). This assumes that \code{uc} is the 
uncertainty in a measurement result \eqn{y}, where \eqn{y=f(x_1, x_2, \dots)}, \code{ci} are 
the partial derivatives \eqn{\partial y/\partial x_i}{dy/dx[i]} and \code{ui} is the standard uncertainty associated with \code{xi}.

The implementation assumes that the combined uncertainty \code{uc} is equal to
\code{sqrt(sum((ci*ui)^2)}. An independent estimate of \code{uc} can be provided.

The \code{ci} are 'sensitivity coefficients'; the default is 1, so that the \code{ui}
can be given either as standard uncertainties in the values of influence quantities \eqn{x_i}, 
together with the associated \code{ci}, or as contributions \code{ci*ui} to the uncertainty in \eqn{y}.

Correlation is not supported, because the Welch-Satterthwaite equation is only valid 
for independent variances.
}
\value{
The calculated effective degrees of freedom associated with \code{uc}.
}
\references{
JCGM 100 (2008) \emph{Evaluation of measurement data - Guide to the expression
of uncertainty in measurement}. \url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf}. 
(JCGM 100:2008 is a public domain copy of ISO/IEC \emph{Guide to the expression
of uncertainty in measurement} (1995) ). 

Satterthwaite, F. E. (1946), An Approximate Distribution of Estimates of Variance Components., 
Biometrics Bulletin \bold{2}, 110-114, doi:10.2307/3002019 

Welch, B. L. (1947), The generalization of "Student's" problem when several different 
population variances are involved., Biometrika \bold{34} 28-35

}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
None, yet.
}
\examples{

u <- c(0.1, 0.3, 0.2, 1.1)
ci <- c(1.0, 2.0, 3.0, 0.5)
degfree <- c(Inf,6,8,3)

w.s(ui=u,df=degfree, ci=ci)

}
\keyword{ univar }


\name{uncert}
\Rdversion{1.1}
\alias{uncert}
\alias{uncert.default}
\alias{uncert.expression}
\alias{uncert.function}
\alias{uncert.formula}
\title{
Uncertainty estimation functions
}
\description{
Functions for estimating measurement uncertainty from standard uncertainties
and either sensitivity coefficients or (for some methods) expressions or functions. 
Correlation is supported via either a correlation or covariance matrix.}
\usage{
   uncert(obj, \dots)

   \method{uncert}{default}(obj, c, method = c("GUM", "MC"), cor, 
            cov, distrib=NULL, distrib.pars=NULL, B=200, x=NULL, keep.x = TRUE, 
            u=obj, \dots)

   \method{uncert}{expression}(obj, x, u, method=c("GUM", "NUM", "kragten", "k2", "MC"), 
            cor, cov, distrib=NULL, distrib.pars=NULL, 
            B=200, delta=0.01, keep.x = TRUE, \dots) 

   \method{uncert}{function}(obj, x, u, method=c("NUM", "kragten", "k2", "MC"), 
            cor, cov, distrib=NULL, distrib.pars=NULL, 
            B=200, delta=0.01, keep.x = TRUE, \dots)

   \method{uncert}{formula}(obj, x, u, method=c("GUM", "NUM", "kragten", "k2", "MC"), 
            cor, cov, distrib=NULL, distrib.pars=NULL, 
            B=200, delta=0.01, keep.x = TRUE, \dots)

}
\arguments{
  \item{obj}{An R object used for method dispatch; see below. Methods currently exist for 
     numeric vector, expression, function, or formula}
  \item{u}{For the default method, a numeric vector of standard uncertainties. For 
     the formula or expression methods, a named list of standard uncertainties.Note that 
     for the default method, \code{u} is set to the value of \code{obj}, allowing 
     specification of either as the \emph{first} argument}
  \item{c}{A numeric vector of senstivity coefficients.}
  \item{x}{For the expression or formula methods, an R object which can be used 
     as an environment by \code{eval}.  For the function method, a list of 
     parameters supplied to FUN via \code{do.call}.}
  \item{method}{Method of uncertainty evaluation. The current list of methods is:
     \describe{
     \item{GUM}{First-order error propagation (the \dQuote{law of propagation of uncertainty}) 
         as implemented by the GUM.}
     \item{NUM}{Numerical differentiation using a simple small step size.}
     \item{kragten}{Numerical estimation of uncertainty following Kragten (Kragten (1994).}
     \item{k2}{A symmetric modification of Kragten's approach described by Ellison (Ellison (2005)).}
     \item{MC}{Monte Carlo simulation.}
     }
  }
  \item{cor, cov}{A (square, symmetric) correlation or covariance matrix, respectively.
     If neither is specified, \code{cor} is set to the identity matrix.}
  \item{distrib}{For \code{method="MC"}, a character vector of length 
  	\code{length(x)} or a named list of names of distribution functions associated 
  	with \code{u}. See Details for defaults. The list format may include 
  	user-specified functions. Silently ignored for other methods.}
  \item{distrib.pars}{For \code{method="MC"}, a list of 
  	lists of parameters describing the distributions associated with \code{u}
  	to be passed to the relevant distribution function. If \code{distrib} is present
  	but \code{distrib.pars} is not, neither are included in the return value
  	unless \code{method="MC"}. See Details for defaults when \code{method="MC"}.
  	Silently ignored for other methods.}
  \item{B}{Number of Monte Carlo replicates.}
  \item{delta}{Step size for numerical differentiation.}
  \item{keep.x}{For \code{method="MC"}, if \code{keepx=TRUE}, the simulated replicates of 
     \code{x} are included in the return object.}
  \item{\dots}{Additional parameters to be passed to a function (for the function method)
   or used in an expression (for expression or formula method).}
}
\details{The default \dQuote{GUM} method applies first-order error propagation principles to estimate 
a combined standard uncertainty from a set of sensitivity coefficients and \emph{either} a 
set of standard uncertainties and a correlation matrix (which defaults to an identity matrix) 
\emph{or} a covariance matrix. Both options use the same calculation, which is simply 
\code{ (t(c) \%*\% cov) \%*\% c }; standard uncertainties are  first combined with 
the correlation matrix provided to form the covariance matrix. Since the correlation 
matrix defaults to the identity matrix, the default is combination without 
correlation.

The default method takes \code{obj} as a vector of uncertainty contributions unless \code{u} 
is specified, in which case \code{u} is used. It is not necessary to specify both. 
The expression method requires \code{obj} to be a differentiable R expression which can 
be evaluated in the environment \code{x} to provide a numeric value.
For the function method, \code{obj} must be an R function which takes parameters from x and 
returns a numeric value.
For the formula method, \code{obj} must be a formula with no left-hand side (e.g. \code{~a*x+b*x^2}) 
which can be evaluated in the environment \code{x} to provide a numeric value.

The formula and expression methods first calculate derivatives for the expression or formula, 
evaluate them using the supplied values of \code{x} and then pass the resulting sensitivity 
coefficients, with supplied \code{u}, \code{cor} or \code{cov} to uncert.default.

The derivatives for the \dQuote{GUM} method (formula and expression methods only) are algorithmic 
derivatives (that is, algebraic or analytical derivatives) obtained using \code{deriv} 
applied to \code{expr} and \code{formula}.  

Numerical derivatives are computed in different ways depending on the method specified: 

- For \code{method="NUM"}, the derivatives are calculated as 
 \eqn{(f(x+delta*u)-f(x-delta*u))/(2*delta*u)}.

- For \code{method="kragten"}, derivatives are calculated as 
 \eqn{(f(x+u*sign(delta))-f(x))/u}.

- For \code{method="k2"}, derivatives are calculated as 
 \eqn{(f(x+u)-f(x-u))/(2*u)}.

\code{"NUM"} is likely to give a close approximation to analytical differentiation provided that 
\code{delta} is appreciably less than 1 but not so small as to give step sizes near machine 
precision. \code{"k2"} is equivalent to \code{"NUM"} with \code{delta=1.0}. Both will give zero coefficients
at stationary points (e.g minima), leading to under-estimation of uncertainty if 
the curvature is large. \code{"kragten"} uses a \emph{deliberately} one-sided (and large) step to 
avoid this problem; as a result, \code{"kragten"} is a poorer (sometimes much poorer) estimate of
the analytical differential but likely a better approximation to the truth.

Since these methods rely on \code{u}, if \code{u} is unspecified and \code{cov} is 
provided, \code{u} is extracted from \code{cov} (using \code{sqrt(diag(cov))}). It is
assumed that the row and column order in \code{cov} is identical to the order of named 
parameters in \code{x}.

Derivatives (and uncertainty contributions) are computed for all parameters in 
\code{x}. Additional parameters used in \code{FUN}, \code{expr} or \code{formula}
may be included in \code{\dots}; these will be treated as constants in the 
uncertainty calculation. 

If \code{distrib} is missing, or if it is a list with some members missing, the distribution 
is assumed Normal and \code{distrib$name} is set to \code{"norm"}. Similarly, if \code{distrib.pars}
or a member of it is missing, the default parameters for \code{x$name} are 
\code{list(mean=x$name, sd=u$name)}. If the list is not named, \code{names(x)} are used 
(so the list must be in order of \code{names(x)}).

If \code{method="MC"}, \code{uncert} calls \code{uncertMC}. Distributions and 
distribution parameters are required and B must be present and >1. See \code{\link{uncertMC}}
for details of distribution specification. 

For other evaluation methods, the distributions are silently ignored.

}

\value{
An object of class \sQuote{uncert} or, for \code{method="MC"} of class \sQuote{uncertMC}. 
See \code{\link{uncert-class}} and \code{\link{uncertMC-class}} for details.

  
}
\references{
JCGM 100 (2008) \emph{Evaluation of measurement data - Guide to the expression
of uncertainty in measurement}. \url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf}. 
(JCGM 100:2008 is a public domain copy of ISO/IEC \emph{Guide to the expression
of uncertainty in measurement} (1995) ). 

Kragten, J. (1994) Calculating standard deviations and confidence intervals with 
a universally applicable spreadsheet technique, Analyst, \bold{119}, 2161-2166.

	
Ellison, S. L. R. (2005) Including correlation effects in an improved spreadsheet 
calculation of combined standard uncertainties, Accred. Qual. Assur. \bold{10}, 338-343.

}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}
\seealso{
  \code{\link{uncert-class}},  \code{\link[stats]{deriv}}
  
  For \code{method="MC"} see \code{\link{uncertMC}} and \code{\link{uncertMC-class}}.
}
\examples{

  expr <- expression(a+b*2+c*3+d/2)
  x <- list(a=1, b=3, c=2, d=11)
  u <- lapply(x, function(x) x/10)
  u.expr<-uncert(expr, x, u, method="NUM")
  u.expr

  #Compare with default:
  uncert(u=c(0.1, 0.3, 0.2, 1.1), c=c(1.0, 2.0, 3.0, 0.5))
  
  #... or with function method
  f <- function(a,b,c,d) a+b*2+c*3+d/2
  u.fun<-uncert(f, x, u, method="NUM")
  u.fun

  #.. or with the formula method
  u.form<-uncert(~a+b*2+c*3+d/2, x, u, method="NUM")
  u.form
  
  #An example with correlation
  u.cor<-diag(1,4)
  u.cor[3,4]<-u.cor[4,3]<-0.5
  u.formc<-uncert(~a+b*2+c*3+d/2, x, u, method="NUM", cor=u.cor)
  u.formc
  
  #A Monte Carlo example
  #See uncertMC for a less linear example
  u.formc.MC<-uncert(~a+b*2+c*3+d/2, x, u, method="MC", cor=u.cor, B=200)
  u.formc.MC
}

\keyword{ univar }


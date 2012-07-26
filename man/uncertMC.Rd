\name{uncertMC}
\Rdversion{1.1}
\alias{uncertMC}
\title{
Monte Carlo evaluation of measurement uncertainty.
}
\description{
\code{uncertMC} estimates measurement uncertainty from a function,
expression or formula by Monte Carlo simulation.}
\usage{
uncertMC(expr, x, u, method = "MC", df, cor, cov, distrib, distrib.pars, 
   B = 200, keep.x = TRUE, vectorized=TRUE, ...)
}
\arguments{
  \item{expr}{An expression, function, or formula with no left-hand side (e.g. 
     \code{~a*x+b*x^2}) which can be evaluated in the environment \code{x} to
     provide a numeric value.}
  \item{x}{A named list or vector of parameters supplied to \code{expr}.}
  \item{u}{A named list or named vector of length \code{length(x)} of standard uncertainties.}
  \item{method}{Method of uncertainty evaluation. The only method currently supported
     by \code{uncertMC} is \code{"MC"}. If any other method is specified, control is 
     passed to \code{uncert}.}
  \item{df}{A named list or named vector of degrees of freedom. \code{df} can be
     a partial named list if not all distributions (see below) use degrees of 
     freedom.}
  \item{cor, cov}{Optional (square, symmetric) correlation or covariance matrices, respectively.
     If neither is specified, \code{uncertMC} assumes independent variables.}
  \item{distrib}{A character vector of length \code{length(x)} or a named list 
        of names of distribution functions associated with \code{u}. See Details
        for defaults. 
        }
  \item{distrib.pars}{A named list of lists of parameters describing the distributions 
  	associated with \code{u} to be passed to the relevant distribution function.
  	If \code{distrib} is present but \code{distrib.pars} is not, an attempt is made 
  	to set defaults based on other parameters; see Details.}
  \item{B}{Number of Monte Carlo replicates.}
  \item{keep.x}{If \code{TRUE}, the simulated replicates of \code{x} are included in 
     the return object.}
  \item{vectorized}{If \code{TRUE}, \code{expr} is assumed to take vector arguments. If \code{FALSE}, 
     \code{expr} is treated as if it takes scalar arguments. See Details for the difference.}
  \item{\dots}{Additional parameters to be passed to a function (for the function method)
   or used in an expression (for expression or formula method).}
}

\details{
Although most likely to be called by \code{uncert}, \code{uncertMC} may be called directly.

If any of \code{x}, \code{u}, \code{df}, \code{distrib} or \code{distrib.pars} are not lists,
they are coerced to lists. If \code{x} is not named, arbitrary names of the form 'Xn' 
are applied. If \code{u}, \code{df}, \code{distrib} or \code{distrib.pars} do not have
names, the names will be set to \code{names(x)} if they are of length exactly 
\code{length(x)}; if not, an error is returned.

For Monte Carlo evaluation, distributions and distribution parameters are needed but 
defaults are used if some or all are absent. If \code{distrib} is missing, or 
if it is a list with some members missing, the distribution is assumed Normal
and any missing member of \code{distrib} is set to "norm". 

Distributions are usually identified by the root of the distribution function name; for example 
to specify the Normal, \code{distrib$name="norm"}. At present, only the random value 
generator (e.g. \code{rnorm}) is used. Names of user-specified distributions functions can also be 
used, provided they have a random value generator named \code{r<dist>} where \code{<dist>}
is the abbreviated distribution. Parameters are passed to distribution functions using 
\code{do.call}, so the function must accept the parameters supplied in \code{distrib.pars}.
%% A function definition may also be included in the list, but 
%% will be assigned an arbitrary name for printing.

If \code{distrib.pars} or members of it are missing, an attempt is made to deduce 
appropriate distribution parameters from \code{x}, \code{u}, \code{df} and \code{distrib}. 
In doing so, the following assumptions and values apply for the respective distributions:
\describe{
   \item{norm}{\code{mean=x$name, sd=u$name}.}
   \item{unif}{\code{min=x-sqrt(3)*u, max=x+sqrt(3)*u}.}
   \item{tri}{\code{min=x-sqrt(6)*u, max=x+sqrt(6)*u, mode=x}.}
   \item{t, t.scaled}{\code{df=df, mean=x, sd=u}.}
}

If either \code{cor} or \code{cov} are present, a test is made to see if off-diagonal
elements are significant. If not, \code{uncertMC} treats the values as independent.
The test simply checks whether the sum of off-diagonal elements of \code{cor} (calculated 
from \code{cov} if \code{cov} is present) is bigger than 
\code{.Machine.double.eps*nrow^2}.

Correlation is supported as long as all correlated variables are normally distributed.
If correlation is present, \code{uncertMC} follows a two-stage simulation procedure. 
First, variables showing correlation are identified. Following a check that 
their associated \code{distrib} values are all \code{"norm"}, \code{mvrnorm} from 
the MASS library is called to generate the simulated \code{x} values for those variables.
Second, any remaining (i.e. independent) variables are simulated from their respective 
\code{distrib} and \code{distrib.pars}. 

Vectorisation makes a difference to execution speed. If \code{vectorize=TRUE}, MC evaluation 
uses \code{eval} using the simulated data as the evaluation environment; if not, \code{apply}
is used row-wise on the simulated input matrix. This makes an appreciable difference to 
execution speed (typically \code{eval} is faster by a factor of 5 or more) so the default 
assumes vectorised expressions. However, not all functions and expressions take vector arguments, 
especially user functions involving complicated arithmetic or numerical solutions. Use \code{vectorize=FALSE}
for functions or expressions that do not take vector arguments. 
Note: One common symptom of an expression that does not take vector arguments is
an R warning indicating that only the first element (typically of a parameter in \code{x}) is used. 
uncertMC may also return NA for \code{u} on attempting to take the sd of a single simulated point.

}
\value{
An object of class \code{uncertMC}. See \code{\link{uncertMC-class}} for details.
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
\code{\link{uncert}},  \code{\link{uncert-class}},   \code{\link{uncertMC-class}}
}
\examples{

  expr <- expression(a+b*2+c*3+d/2)
  x <- list(a=1, b=3, c=2, d=11)
  u <- lapply(x, function(x) x/10)
  u.MC<-uncertMC(expr, x, u, distrib=rep("norm", 4), method="MC")
  print(u.MC, simplify=FALSE)

  #An example with correlation
  u.cor<-diag(1,4)
  u.cor[3,4]<-u.cor[4,3]<-0.5
  u.formc.MC<-uncertMC(~a+b*2+c*3+d/2, x, u, cor=u.cor, keep.x=TRUE)
  u.formc.MC

  #A non-linear example
  expr <- expression(a/(b-c))
  x <- list(a=1, b=3, c=2)
  u <- lapply(x, function(x) x/20)
  set.seed(403)
  u.invexpr<-uncertMC(expr, x, u, distrib=rep("norm", 3), B=999, keep.x=TRUE )
  u.invexpr

  #Look at effect of vectorize
  system.time(uncertMC(expr, x, u, distrib=rep("norm", 3), B=9999, keep.x=TRUE ))
  system.time(uncertMC(expr, x, u, distrib=rep("norm", 3), B=9999, keep.x=TRUE, vectorize=FALSE))
  
}

\keyword{ univar }

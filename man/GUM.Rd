\name{GUM}

\alias{GUM}

\title{
Propagation of Measurement Uncertainty for Typical Metrology Applications 
Using the Methods Outlined in the GUM
}

\description{
A function for propagation of measurement uncertainty for typical metrology 
applications using the methods from the Joint Committee on Guides in Metrology (JCGM)
\emph{Guide to the Expression of Uncertainty in Measurement (GUM)}. This approach 
approximates the uncertainty of a function of random variables that define a 
measurement result by computing the uncertainty of the first-order Taylor series for
the measurement function. This function also serves as the primary computational tool
underlying the GUM uncertainty templates found in the metRology for Microsoft Excel user interface.
}

\usage{
GUM(var.name, x.i, u.i, nu.i, measurement.fnc, correlation = diag(length(var.name)),
shared.u.i = var.name, cl = 0.95, cov.factor = "Student's t", sig.digits.U = 2, ...)
}

\arguments{
  \item{var.name}{Character vector of input variable names.}
  \item{x.i}{Vector of input variable values.}
  \item{u.i}{Vector of standard uncertainties (i.e. standard errors) for each input variable value.}
  \item{nu.i}{Degrees of freedom associated with each standard uncertainty.}
  \item{measurement.fnc}{Character string specifying the functional relationship between
                         input variables that defines the output measurement result.}
  \item{correlation}{Matrix giving the correlation between the different input variable values. 
                     Default is to assume no correlation between input variable values.}
  \item{shared.u.i}{Character vector giving the relative relationship between the standard uncertainties
                    for each variable value. Groups of variables based on a common shared standard uncertainty
                    will all share the same variable name. The default is to assume all standard uncertainties
                    are assessed independently, resulting a value of shared.u.i that is identical to var.name.}
  \item{cl}{Nominal confidence level to be used to compute the expanded uncertainty of the output measurement result.
             Default value is 0.95.}
  \item{cov.factor}{Type of coverage factor to be used. The default is to use the value from the Student's t 
                    distribution with confidence level specified above and nu.eff effective degrees of freedom.}
  \item{sig.digits.U}{Number of significant digits to be reported in the expanded uncertainty of 
                      the measurement result. The measurement result will be rounded to the same
                      number of decimal places.}
  \item{\dots}{ Arguments passed to other functions. Currently unimplemented.}
}

\details{
Whenever possible, sensitivity coefficients are obtained analytically using the gradient attribute of the 
\code{\link{deriv}} function. In situations where some part of the measurement function is not found in
derivative table, sensitivity coefficients are obtained by numeric partial differentiation using the
\code{\link{grad}} function from the package numDeriv.
}

\value{
A list containing the 9 components:
\item{y}{Value of the measurement result obtained by evaluating the measurement function at the input variable values.}
\item{uc}{The combined standard uncertainty of the measurement result, y.}
\item{nu.eff}{The effective degrees of freedom associated with uc, computed using the Welch-Satterthwaite formula.}
\item{cl}{The nominal confidence level used to obtain the coverage factor, k.}
\item{k}{The coverage factor used to control the confidence level associated with the expanded uncertainty of the
         measurement result.}
\item{U}{The expanded uncertainty of the measurement result, computed as U=k*uc.}
\item{contributions}{Relative variance contributed to the standard uncertainty (uc) of the measurement result
                     from each input variable.}
\item{sensitivities}{Sensitivity coefficient associated with each input variable.}
\item{msgs}{Error and warning messages that point out potential problems with the inputs to the \code{GUM} function
            or with the interpretation of the function's output.}
}

\references{
Joint Committee on Guides in Metrology (JCGM), \emph{Evaluation of Measurement Data Guide to the Expression of
Uncertainty in Measurement}, http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf, 2008.
}

\author{
Hung-kung Liu \email{hung-kung.liu@nist.gov} and Will Guthrie \email{will.guthrie@nist.gov}
}

\seealso{
\code{\link{GUM.validate}} a function to assess the statistical performance of GUM uncertainty intervals
for the application of interest in terms of average attained coverage probability.
\code{\link{uncert}} for a family of functions focused on the study and comparison of different approaches
and numerical options in uncertainty analysis.
}

\examples{
## a simple uncertainty analysis for the product of two quantities
GUM(c("x1","x2"),c(2.3,1.1),c(0.030,0.015),c(5,9999),"x1*x2")

## example of the difference in the measurements of two standards, each 
## with a standard uncertainty based on a common value drawn from a control chart
## representative of the measurement process made using a check standard that 
## is comparable to the two individual standards under study
GUM(c("s1","s2"),c(45.3,46.0),c(0.26,0.26),c(134,134),"s1-s2",shared.u.i=c("s1","s1"))

## compare with results for equivalent, alternative specification of shared.u.i
GUM(c("s1","s2"),c(45.3,46.0),c(0.26,0.26),c(134,134),"s1-s2",shared.u.i=c("s2","s2"))
}

\keyword{univar}
\keyword{htest}

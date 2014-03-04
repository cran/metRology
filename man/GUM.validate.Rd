\name{GUM.validate}

\alias{GUM.validate}

\title{
Monte Carlo Check on the Statistical Performance of GUM Uncertainty Intervals 
Using Attained Coverage Probability
}

\description{
A function for assessing the statistical performance of measurement uncertainty intervals for 
particular metrology applications computed using the methods from the Joint Committee on 
Guides in Metrology (JCGM) \emph{Guide to the Expression of Uncertainty in Measurement (GUM)}. 
The validation is carried out using the input values as true values in a simulation that 
directly checks the attained coverage probability of the uncertainty
intervals produced using the \code{\link{GUM}} function.
}

\usage{
GUM.validate(var.name, x.i, u.i, nu.i, type, distribution, measurement.fnc, 
             correlation = diag(length(var.name)), shared.u.i = var.name, cl = 0.95, 
             cov.factor = "Student's t", sig.digits.U = 2)
}

\arguments{
  \item{var.name}{Character vector of input variable names.}
  \item{x.i}{Vector of input variable values.}
  \item{u.i}{Vector of standard uncertainties (i.e. standard errors) for each input variable value.}
  \item{nu.i}{Degrees of freedom associated with each standard uncertainty.}
  \item{type}{Character vector of values "A" and "B" indicating the methods used to evaluate the standard
              uncertainty of each input value. Standard uncertainties evaluated using statistical methods
              are denoted Type A in the GUM, while standard uncertainties evaluated using other means are
              denoted Type B.}
  \item{distribution}{Character vector of probability distributions associated with the potential values
                      taken on by each input variable. The current possible choices are "Normal" (i.e. Gaussian),
                      "Triangular", or "Rectangular" (i.e. Uniform).}
  \item{measurement.fnc}{Character string specifying the functional relationship between
                         input variables that defines the output measurement result.}
  \item{correlation}{Matrix giving the correlation between the different input variable values. 
                     Default is to assume no correlation between input variable values.}
  \item{shared.u.i}{Character vector giving the relative relationship between the standard uncertainties
                    for each variable value. Groups of variables based on a common shared standard uncertainty
                    share will all share the same variable name. The default is to assume all standard uncertainties
                    are assessed independently, resulting a value of shared.u.i that is identical to var.name.}
  \item{cl}{Nominal confidence level to be used to compute the expanded uncertainty of the output measurement result.
             Default value is 0.95.}
  \item{cov.factor}{Type of coverage factor to be used. The default is to use the value from the Student's t 
                    distribution with confidence level specified above and nu.eff effective degrees of freedom.}
  \item{sig.digits.U}{Number of significant digits to be reported in the expanded uncertainty of the measurement result.
                      The measurement result will be rounded to the same number of decimal places.}
}

\details{
Currently 1000 simulated sets of uncertainty data are used for the computation of the attained confidence level.
}

\value{
A Monte Carlo assessment of the attained coverage of expanded uncertainty intervals like those produced 
using the \code{\link{GUM}} function for the application of interest.
}

\references{
Joint Committee on Guides in Metrology (JCGM), \emph{Evaluation of Measurement Data Guide to the Expression of
Uncertainty in Measurement}, http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf, 2008.
}

\author{
Hung-kung Liu \email{hung-kung.liu@nist.gov} and Will Guthrie \email{will.guthrie@nist.gov}
}

\seealso{
\code{\link{GUM}} a function to compute GUM uncertainty intervals for general metrological applications.
}

\examples{
## a simple uncertainty analysis for the product of two quantities
GUM.validate(c("x1","x2"), c(2.3,1.1), c(0.030,0.015), c(5,9999),
               c("A","B"),c("Normal","Rectangular"),"x1*x2")
}

\keyword{univar}
\keyword{htest}

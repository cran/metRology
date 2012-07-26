% file metRology/man/GUM.H.1.Rd
% Will Guthrie
\name{GUM.H.1}
\docType{data}
\alias{GUM.H.1}
\title{
Example H.1 from the \emph{Guide to the Expression of Uncertainty in Measurement}
}
\description{
Calibration of an end gauge for length measurement. Notation is based on the presentation in chapter 3 of 
\emph{Data Modeling for Metrology and Testing in Measurement Science}.
}
\usage{
GUM.H.1
data(GUM.H.1)
}
\format{
A list containing 10 components that give the uncertainty budget, measurement function, and other quantities
needed to carry out an uncertainty analysis using the methods from the Guide to the Expression of Uncertainty
in Measurement:
\describe{
\item{var.name}{Character vector giving the name of each input variable included in the analysis.}
\item{unit}{Expression giving the units of each input variable}
\item{x.i}{Vector giving the reported value for each input variable.}
\item{u.i}{Vector giving the standard uncertainty associated with each reported value.}
\item{nu.i}{Vector giving the degrees of freedom associated with each standard uncertainty.}
\item{type}{Character vector indicating the method of evaluation (Type A or Type B) for each standard uncertainty.}
\item{distribution}{Character vector listing the probability distribution assumed to describe each variable value.}
\item{measurement.fnc}{Character string giving the measurement function for the output variable.}
\item{correlation}{Matrix giving the correlations between input variable values.}
\item{shared.u.i}{Character vector describing which standard uncertainties, if any, are based on a common 
      underlying standard uncertainty. A vector that is the same as var.name indicates that all standard
      uncertainties are based on independent assessments of standard uncertainty. A vector with fewer names
      than in var.names indicates that one or more variables are derived from a common uncertainty assessment.}
}
}
\details{
Details can be found in the GUM and in Data Modeling for Metrology and Testing in Measurement Science.
}
\source{
Joint Committee on Guides in Metrology (JCGM), \emph{Evaluation of Measurement Data Guide to the Expression of
Uncertainty in Measurement}, http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf, 2008.
}
\references{
Guthrie, W.F. et. al., "Three Statistical Paradigms for Assessment and Interpretation of Measurement Uncertainty" in
\emph{Data Modeling for Metrology and Testing in Measurement Science}, F. Pavese and A.B. Forbes, eds., Birkhauser,
Boston, 2009.
}
\keyword{datasets}

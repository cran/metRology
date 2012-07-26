% file metRology/man/Pb.Rd
% (c) 2009 S Ellison
\name{Pb}
\docType{data}
\alias{Pb}
\title{
Lead in wine
}
\description{
A data frame containing reported results for lead (in mg/kg) from 
CCQM Key Comparison CCQM-K30.
}
\usage{
	Pb
}
\format{
A data frame containing 11 reported results with uncertainty data:
	\describe{
		\item{lab}{Factor giving abbreviated laboratory identifier}
		\item{value}{The reported value for lead (mg/kg)}
		\item{u}{Standard uncertainty (mg/kg). The values in \code{Pb} were calculated from 
		the reported expanded uncertainty \code{U} and coverage factor \code{k} using \code{u=U/k}. }
		\item{k}{Coverage factor. Conventionally, the coverage factor is set to a suitable  
			quantile of Student's t based on the Welch-Satterthwaite effective degrees of freedom 
			or simply set to 2 for approximately 95\% confidence. In this data set, labs all quoted 
			\code{k} for approximately 95\% confidence. }
		\item{U}{Expanded uncertainty as reported by labs.}
		\item{method}{Factor indicating general measurement methodology:
			\describe{
				\item{IDMS}{Isotope dilution mass spectrometry}
				\item{ICP}{Inductively coupled plasma spectrometry}
				\item{GFAAS}{Graphite furnace atomic absorbtion spectrometry}
			}
		}
		\item{include}{logical; Whether the reported result was included in the calculation 
		of the Key Comparison Reference Value for the study. }
	}
}
\details{
The study involved circulation of a homogeneous set of samples of wine for analysis for lead (Pb) content
by a number of National Measurement Institutes. 

The Key Comparison Reference Value, or assigned value for the lead content, was
set at 2.99 mg/kg with expanded uncertainty 0.06 mg/kg.

}
\source{
Hearn, R., Santamaria-Fernandez, R. and Sargent, M. (2008) 
Final report on key comparison CCQM-K30: Determination of lead in wine.
\emph{Metrologia}
\bold{45}, 08001, 2008
}
\references{
See source.
 
}
\keyword{datasets}

\name{potassium}
\alias{potassium}
\docType{data}
\title{
Potassium data for two different materials included in an interlaboratory study
}
\description{
Potassium data for two different materials included in an interlaboratory study 
intended to provide data for certification of a reference material.
}
\usage{data("potassium")}
\format{
  A data frame with 25 observations on the following 2 variables.
  \describe{
    \item{\code{QC}}{Potassium concentrations (mg/kg) reported on a 
    	material used as a quality control material}
    \item{\code{RM}}{Potassium concentrations (mg/kg) reported on a 
    	candidate reference material material used as a quality control material}
  }
}
\details{
Potassium data for two different materials included in an interlaboratory study 
intended to provide data for certification of a crab tissue reference material. 
The study included a previously certified reference material (near end of stock) to 
serve as a quality control (QC) check. Laboratories were asked to report five replicate 
measurements on the candidate reference material and three for the QC material. Each 
row in the data set corresponds to the mean of replicate results reported by each laboratory. 

Inspection of the data suggests that one laboratory interchanged or mislabelled 
the test materials. The anomalous results appear as an outlier for both QC and RM, as well as
being visible as an off-diagonal outlier in a Youden plot - see \code{\link{youden.plot}}).
}
\source{
	Private communication - Pending publication
}
\examples{
data(potassium)
yplot(potassium)
}
\keyword{datasets}

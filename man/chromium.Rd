\name{chromium}
\alias{chromium}
\docType{data}
\title{
Chromium data for two different materials included in an interlaboratory study
}
\description{
Chromium data for two different materials included in an interlaboratory study 
intended to provide data for certification of a reference material.
}
\usage{data("chromium")}
\format{
  A data frame with 28 observations on the following 2 variables.
  \describe{
    \item{\code{QC}}{Chromium concentrations (ug/kg) reported on a 
    	material used as a quality control material}
    \item{\code{RM}}{Chromium concentrations (ug/kg) reported on a 
    	candidate reference material material used as a quality control material}
  }
}
\details{
Chromium data for two different materials included in an interlaboratory study 
intended to provide data for certification of a crab tissue reference material. 
The study included a previously certified reference material (near end of stock) to 
serve as a quality control (QC) check. Laboratories were asked to report five replicate 
measurements on the candidate reference material and three for the QC material. Each 
row in the data set corresponds to the mean of replicate results reported by each laboratory. 

Inspection of the data suggests that one laboratory interchanged or mislabelled 
the test materials; this is hard to see in univariate plots but relatively easy 
to see in a Youden plot (a type of pairwise scatter plot - see \code{\link{youden.plot}}).
}
\source{
	Private communication - Pending publication
}
\examples{
data(chromium)
yplot(chromium)
}
\keyword{datasets}

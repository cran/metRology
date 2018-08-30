% file metRology/man/RMstudy.Rd
% (c) 2011 S Ellison
\name{RMstudy}
\docType{data}
\alias{RMstudy}
\title{
Collaborative study results for metals in a reference material certification study
}
\description{
A data frame containing reported replicate results, in mg/L, for elements reported 
in an inter-laboratory certification study for a candiate drinking reference material.
}
\usage{
	data(RMstudy)
}
\format{
  A data frame with 145 observations on the following 9 variables.
  \describe{
    \item{\code{Lab}}{a factor with levels \code{Lab1} - \code{Lab29}}
    \item{\code{Arsenic}}{a numeric vector}
    \item{\code{Cadmium}}{a numeric vector}
    \item{\code{Chromium}}{a numeric vector}
    \item{\code{Copper}}{a numeric vector}
    \item{\code{Lead}}{a numeric vector}
    \item{\code{Manganese}}{a numeric vector}
    \item{\code{Nickel}}{a numeric vector}
    \item{\code{Zinc}}{a numeric vector}
  }
}
\details{
The data set includes results for eight of a total of 23 elements studied in an 
inter-laboratory study of a candidate reference material.  All observations are 
reported in ug/l. 
Laboratories were asked to report 5 replicate observations for each element. Replicate 
observations appear on separate rows. Most but not all laboratories reported five 
replicates, and some laboratories reported no results for some elements. The eight
elements included in the data set are those for which no more than three laboratories
failed to report any results. Missing observations are coded \code{NA}.

Laboratories were coded anonymously in order of receipt of results.

}
\source{
LGC limited, Teddington, UK (Private communication).
}
\keyword{datasets}

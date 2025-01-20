\name{cov.dellipse}
\alias{cov.dellipse}
\title{
Constructs a covariance and location object for use in plotting data ellipses. 
}
\description{
Constructs a covariance matrix and associated location using a variety 
of (possibly robust) estimators. The returned object is suitable for use 
by \code{\link{plot.d.ellipse}}.
}
\usage{
cov.dellipse(x, y = NULL, cov.method = c("spearman", "kendall", "pearson", 
                           "MCD", "OGK", "GK", "gk", "rgk", "mcd", "mve"), 
                           scalefn = NULL, locfn = NULL, cov.control = list())
}

\arguments{
  \item{x}{
	An R \code{numeric} object. Can be a vector (in which case \code{y} must be specified 
	and of the same length) or a two-column numeric matrix.
}
  \item{y}{
	A numeric vector of the same length as \code{x}. It is an error to provide \code{y} 
	in addition to a two-column matrix for \code{x}.
}
  \item{cov.method}{
A character value specifying the covariance method used.
}
  \item{scalefn}{
	A function that computes univariate scale and (optionally) location estimates from a 
	numeric vector. 
	If provided, \code{scalefn()} should return a single numeric value containing a scale 
	(standard deviation) estimate. For many covariance methods this can be a simple
	scale estimator. For cov.method "GK", scalefn must accept 
	an additional argument \code{mu.too}. When mu.too is true, \code{scalefn()} should 
	return a numeric vector of length 2 containing location and scale estimates. See 
	\code{\link[robustbase]{scaleTau2}}, \code{\link[robustbase]{s_Qn}},\code{\link[robustbase]{s_mad}}, or \code{\link[robustbase]{s_IQR}} 
	for examples to be used as \code{scalefn} argument.
}
  \item{locfn}{
	A function that computes univariate location estimates from a numeric vector. 
	If used, \code{locfn()} should return a single numeric value containing a location 
	(mean) estimate. 
}
  \item{cov.control}{
	A named list of arguments passed to the covariance calculation used. Note that this can 
	override \code{scalefn} and \code{locfn}; see below for details.
 }
}
\details{
	\code{cov.dellipse} is a wrapper for a range of covariance estimation methods found in 
	various packages. Its operation and defaults depend on the particular covariance 
	estimator specified by \code{cov.method}.  Details for each are as follows.
	\describe{
		\item{\code{spearman}, \code{kendall}}{By default, the median and mad are used as location and scale respectively, 
		                      and the covariance is calculated from the product of scale estimates and the 
		                      Spearman rank correlation or Kendall's tau respectively. 
		                      If either \code{scalefn} or \code{locfn}is supplied, \code{scalefn} is used for scale estimation and 
		                      \code{locfn} for location. For both \code{spearman} and \code{kendall}, \code{scalefn} is 
		                      only used as a scale estimator and need not take a \code{mu.too} argument.}
		\item{\code{pearson}}{By default, the mean and sd are used as location and scale respectively, 
		                      and the covariance is calculated from the product of scale estimates and the 
		                      Pearson correlation. 
		                      If either \code{scalefn} or \code{locfn}is supplied, \code{scalefn} is used for scale 
		                      estimation and \code{locfn} for location, making it possible (if not very sensible) to 
		                      use a combination of robust scale or location functions with the Pearson correlation coefficient.
		                      For this case,  \code{scalefn} is only used as a scale estimator and need 
		                      not take a \code{mu.too} argument.}
		\item{\code{MCD}, \code{mcd}}{Both compute the Minimum Covariance Determinant (MCD) estimator, a robust multivariate 
				      location and scale estimate with a high breakdown point, via the 'Fast MCD' or 'Deterministic MCD' 
				      ("DetMcd") algorithm. \code{"MCD"} uses the implementation \code{\link[robustbase]{covMcd}} in the robustbase package; 
				      \code{"mcd"} uses \code{\link[MASS]{cov.mcd}} in the MASS package.
				      Neither require or use \code{scalefn} or \code{locfn}. 
				      Note that these MCD implementations differ appreciably for small samples (at least to n=60). \code{MCD} 
				      includes consistency and finite sample correction whereas \code{mcd} apparently does not apply a finite 
				      sample correction. As a result, the \code{mcd} scales can be considerably smaller for modest 
				      data set sizes.}
		\item{\code{OGK}}{Computes the orthogonalized pairwise covariance matrix estimate described by Maronna and Zamar (2002),
				as implemented by the \code{\link[robustbase]{covOGK}} in the robustbase package. 
				By default, scale and location use \code{\link[robustbase]{scaleTau2}} from robustbase. Alternatives 
				can be specified either by providing \strong{both} \code{scalefn} and \code{locfn} or by including
				an argument \code{sigmamu} in \code{cov.control}, which is passed to \code{covOGK}. See 
				\code{\link[robustbase]{covOGK}} for a description of \code{sigmamu}.
				If \code{sigmamu} is not present in \code{cov.control} and both \code{scalefn} and \code{locfn} 
				are provided, scale and location are constructed from \code{scalefn} and \code{locfn}. If only one
				of these is provided, a warning is issued and \code{[robustbase]{scaleTau2}} is used.
				}
		\item{\code{GK}}{Computes a simple pairwise covariance estimate suggested by Gnanadesikan and Kettenring (1972),
				as implemented by the \code{\link[robustbase]{covGK}} in the robustbase package. 
				By default, scale and location use \code{\link[robustbase]{scaleTau2}} from robustbase. Alternatives 
				can be specified either by providing \code{scalefn} and \code{locfn} or by including
				an argument \code{scalefn} in \code{cov.control}, which is passed to \code{covGK}. See 
				\code{\link[robustbase]{covGK}} for a description of \code{scalefn}.
				If \code{scalefn} is not present in \code{cov.control}, scale and location are constructed from \code{scalefn} 
				and \code{locfn}. If \code{locfn} is omitted, \code{scalefn} is used if it takes an argument \code{mu.too} 
				and the median is used otherwise.} 
		\item{\code{gk}}{As \code{GK}, except that the variables are scaled to unit (robust) sd (using \code{scalefn}) before 
		                calculating the covariance (which is then rescaled). This can prevent large scale differences from 
		                masking outliers in a variable with small scale.} 
		\item{\code{rgk}}{Implements  Gnanadesikan and Kettenring's second covariance estimate 
			based on scaled variables \eqn{(Z_1, Z_2)}{(Z[1], Z[2])} and a robust correlation \eqn{\rho^*}{rho}
			calculated as 
			\deqn{\rho^*=(\hat{\sigma}_{+}^{*2} - \hat{\sigma}_{-}^{*2})/(\hat{\sigma}_{+}^{*2} - \hat{\sigma}_{-}^{*2})}{%
				rho=(sigma[+]^2 - sigma[-]^2)/(sigma[+]^2 - sigma[-]^2)}
			where \eqn{\hat{\sigma}_{+}^{*2}}{sigma[+]^2} and \eqn{\hat{\sigma}_{-}^{*2}}{sigma[-]^2} are robust variances of 
			\eqn{(Z_1+Z_2)}{Z[1]+Z[2]} and \eqn{(Z_1-Z_2)}{Z[1]-Z[2]} respectively, calculated using \code{scalefn}. 
			The advantage over \code{"gk"} and \code{"GK"} is that the correlation
			coefficient is guaranteed to be in \eqn{[-1,1]}, making for a positive definite covariance matrix. Scaling also
			helps prevent large scale differences from masking outliers in a variable with small scale.} 
		\item{\code{mve}}{Uses uses \code{\link[MASS]{cov.mve}} in the MASS package, which is based on the location and covariance matrix for 
				a minimum volume ellipsoid. The method neither requires nor uses \code{scalefn} or \code{locfn}.}
	}
}
\value{
An object of class \code{cov.dellipse}, which is a list with (at least) components
\describe{	
	\item{method}{Character string describing method; identical to \code{cov.method}}
	\item{cov}{2x2 covariance matrix}
	\item{cor}{2x2 correlation matrix}
	\item{center}{vector (length 2) specifying centre of ellipse}
	\item{scale}{vector, length 2, specifying scale estimates for each variable}
	\item{n.obs}{number of points (rows) used in the covariance estimate}
}
This list is intended to be consistent with that returned by \code{\link{cov.wt}}.

}

\references{
Maronna, R.A. and Zamar, R.H. (2002) Robust estimates of location and dispersion of 
high-dimensional datasets; Technometrics \bold{44(4)}, 307-317.

Gnanadesikan, R. and John R. Kettenring (1972) Robust estimates, residuals, 
and outlier detection with multiresponse data. Biometrics \bold{28}, 81-124
}
\author{
Stephen L R Ellison

}

\seealso{
 \code{\link[MASS]{cov.rob}} in MASS, \code{\link[robustbase]{covMcd}}, \code{\link[robustbase]{covOGK}} and 
 \code{\link[robustbase]{covGK}} in robustbase.

}
\examples{
data(potassium)
cov.dellipse(potassium) #Defaults to Spearman rank correlation

#With different method
cov.dellipse(potassium, cov.method="OGK") 

#Same as above but specifying control parameters
library(robustbase) #For scaleTau2
cov.dellipse(potassium, cov.method="OGK", cov.control=list(sigmamu=scaleTau2)) 
	
#With individually specified (mad) scale
cov.dellipse(potassium, cov.method="GK", scalefn=mad) 
	#Defaults to median for location because mad()
	#does not accept a mu.too argument

cov.dellipse(potassium, cov.method="GK", scalefn=scaleTau2) 
	#Defaults to specified scalefn for location because scaleTau2 
	#accepts mu.too=TRUE

}

\keyword{ robust }
\keyword{multivariate}


\name{buildCor}
\Rdversion{1.1}
\alias{buildCor}
\alias{buildCov}
\alias{updateCor}
\alias{updateCov}
\title{
Functions to build correlation and covariance matrices.
}
\description{
Functions to build and update correlation and covariance matrices using a compact
specification of off-diagonal terms.
}
\usage{
buildCor(n, cors)

updateCor(cor, cors, cor.names)

buildCov(s, covs, vars = s^2, cors, cov.names)

updateCov(cov, covs, cors, cov.names)

}
\arguments{
  \item{n}{scalar: number of rows/colums required in correlation matrix. 
}
  \item{cors, covs}{3-column matrix or data frame specification of individual 
  correlation or covariance terms. Can also be a vector of length 3. See Details.}
  \item{s, vars}{vector of standard deviations or variances respectively. One of 
     \code{s} or \code{vars} must be present.}
  \item{cor.names, cov.names}{vectors of names for the rows and columns of the returned 
      matrix. \code{cov.names} defaults to \code{names(s)} if \code{s} is named.}
  \item{cor, cov}{correlation or covariance matrix requiring amendment.}
}
\details{
For \code{buildCor}, the size of the returned correlation matrix is set using \code{n}; 
an \code{n} by \code{n} correlation matrix is returned. For \code{buildCov} the size is 
set to \code{length(s)} by \code{length(s)}.

Each row of \code{cors} specifies a correlation term \eqn{r_{ij}}{r[i,j]} in the form
\eqn{(i, j, r_{ij})}{(i, j, r[i,j]}. That is, the first two columns give the row 
and column in the desired correlation matrix, and the third gives the relevant 
correlation coefficient. On constructing or updating the correlation matrix, 
\eqn{r_{ij}}{r[i,j]} is set equal to \eqn{r_{ji}}{r[j,i]}, so it is only necessary 
to specify one of \eqn{r_{ij}}{r[i,j]} or \eqn{r_{ji}}{r[j,i]}.

\code{covs} specifies covariance terms in the same way except that the third column of 
\code{covs} must be a covariance.

If either \code{cors}  or \code{covs} is a vector of length 3, it is coerced 
to a matrix of three columns. 

If \code{cor.names} or \code{cov.names} are present, the matrix returned has 
dimnames set to the names supplied.

All four functions test for positive definite return values and generate a 
warning if not positive definite.
}
\value{
A square symmetric correlation or covariance matrix.

}
\references{
None.
}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
None.
}
\examples{
    
    #Duplicate correlation for example for uncert()
    buildCor(4, cors=c(3, 4, 0.5))
               
    
    #Multiple correlations
    r<-buildCor(3, cors=rbind( c(1,2,0.5), c(2,3,0.25) ) )
    r
    
    updateCor(r, cors=c(1,3,0.13)) #perhaps more realistic
    
    buildCov(1:3, cors=rbind( c(1,2,0.5), c(2,3,0.25),c(1,3,0.13) ) )
    
}

\keyword{ univar }

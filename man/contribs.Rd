\name{contribs}
\Rdversion{1.1}
\alias{contribs}
\title{
Extract contributions from an 'uncert' object.
}
\description{
Extracts the individual nonzero contributions 
to the combined uncertainty in an 'uncert' object.
}
\usage{
contribs(object, scope, as.sd = FALSE, keep.sign = TRUE, 
	simplify = TRUE, expand.dot=TRUE)
}
\arguments{
  \item{object}{
An object of class \code{uncert} returned by \code{\link{uncert}} or  \code{\link{uncertMC}}}
  \item{scope}{
An expression, one-sided formula or character vector describing the particular variables
for which contributions are desired. If missing, contributions for all variables
are returned.
}
  \item{as.sd}{ logical; controls whether values are returned in the form of 
  standard uncertainties or variance contributions. See Details.
}
  \item{keep.sign}{logical; controls whether the sign of the cobntributions is appended to
  the return value when \code{as.sd=TRUE}. See Details.
}
  \item{simplify}{ logical. If \code{simplify=FALSE} the contribution matrix itself 
  is returned. If \code{simplify=TRUE}, only the requested (by \code{scope}) 
  nonzero elements of the contribution matrix are returned, as a vector.  
  See Details for the treatment of off-diagonal terms
}
  \item{expand.dot}{logical; if \code{TRUE}, \sQuote{.} in a formula \code{scope} is
  expanded to all contributions including pairwise contributions. If \code{FALSE},
  the dot operator implies only the single-variable terms. See Details.
}

}
\details{
\code{contribs} calculates the contribution matrix \eqn{C} where 
\eqn{C_{i,j}=(c_iu_i)(c_ju_j)r_{i,j}}{C[i,j]=c[i]u[i]*c[j]u[j]*r[i,j]}. 
In general, these values are possibly negative (co)variance contributions 
to the variance (squared standard uncertainty) in \eqn{y}. In GUM notation 
(\sQuote{the GUM} is JCGM 100 (2008) - see references), the diagonal elements of 
\code{C} are squared standard uncertainties in \eqn{y}. The form of the 
return value depends on \code{simplify}, \code{as.sd} and \code{keep.sign}.

If \code{as.sd} is \code{FALSE} (the default), contributions \eqn{C_{ij}}{C[i,j]}
are returned unchanged. For the diagonal elements of \eqn{C} (contributions for 
individual individual terms), this form corresponds to squared standard uncertainties 
\eqn{u_i^2(y)}{(u[i]{y})^2} in GUM notation. 

If \code{as.sd=TRUE}, the magnitude of the value returned is \eqn{\sqrt{|C_{ij}|}}{sqrt( abs(C[i,j] ) )}.
For the diagonal elements of \eqn{C} this corresponds to standard uncertainties 
\eqn{u[i]{y}} in GUM notation.

If \code{as.sd=TRUE}, \code{keep.sign} controls whether the values are signed or 
returned as absolute values. If \code{keep.sign=TRUE}, the value returned is 
\eqn{sign(C_{ij}\sqrt{|C_{ij}|}{sign(C[i,j]sqrt( abs(C[i,j] ) )}}. If false, 
the absolute value is returned. Note that the sign is returned solely to indicate 
the direction of the original contribution. \code{keep.sign} has no effect if 
\code{as.sd=FALSE}.

If \code{simplify=FALSE} (the default), the requested elements of the contribution matrix 
\eqn{C} are returned as a matrix. If \code{simplify=FALSE}, the return value is a vector 
containing only those terms with nonzero values in the associated correlation matrix. 
The threshold for deciding a correlation is nonzero is that its magnitude is greater
than \code{2*.Machine$double.eps}.

Off-diagonal terms for the same pair of variables are summed, that is, for 
the pair \eqn{(C_{ij}, C_{ji}), j \neq i}{(C[i,j], C[j,i]), j != i}
the (single) value returned is \eqn{C_{i,j}+C_{j,i}=2C_{i,j}}{C[i,j]+C[j,i]=2C[i,j]}. 

The contributions returned can be limited to a chosen subset using \code{scope}; 
only the terms involving variables included in \code{scope} are returned. 
\code{scope} can be an expression, formula or character vector of variable names. 
If an expression or formula, only those contributions involving variables in 
the expression or formula are returned. 

Any variable names in \code{scope} which are not present in 
\code{row.names(object$budget)} are silently ignored except for 
the formula specification which will return an error.

If \code{simplify=FALSE}, the matrix returned always contains all contributions 
involving individual variables in \code{scope}. If \code{simplify=TRUE}, however, specifying
\code{scope} as a formula provides additional control over the returned contributions: 

If a formula, \code{scope} accepts the usual model formula operators \sQuote{.}, \sQuote{+}, \sQuote{-},
\sQuote{*} and \sQuote{^}, but the interpretation is not quite identical to \code{lm}.

First, if present, \sQuote{.} is taken by default as \sQuote{all contributions}, implying 
all single terms and all pairwise terms (like \sQuote{\code{.^2})} in other formula specifications). 
This can be disabled by specifying \code{expand.dot=FALSE}.

The negation operator \sQuote{-} removes terms, but removing a single variable also removes any 
associated covariance contributions. For example, \code{scope=~.-A} is expanded to all single 
and pairwise contributions to the uncertanty budget that do not involve \code{A}.

Interaction-like terms of the form \code{A:B} are interpreted as indicating the \emph{total} 
off-diagonal contribution, that is, \code{A:B} is equivalent to \code{B:A} and the associated 
value returned is based on \eqn{C_{i,j}+C_{j,i}}{C[i,j]+C[j,i]}. 

Cross-terms like \code{~A*B} are supported and expand, as usual, to \code{~A+B+A:B}. 

Unlike the two other scope specifications, single terms in the formula do \emph{not}
automatically imply off-diagonal terms; \code{A+B} will not return the off-diagonal contribution for 
\code{A} and \code{B}. Use \code{A*B} or \code{(A+B)^2} etc. to get off-diagonal contributions.  
Cross-terms of order above two are ignored so \code{A*B*C} safely returns only the set of 
individual and pairwise terms, but it is perhaps more precise to use \code{(A+B+C)^2}.

\code{I()} and other operators or functions are not supported.

}

\value{
A named vector or matrix of contributions. Names for off-diagonal contributions 
in the vector format are constructed from the names of the two contributing variables.

}
\references{
JCGM 100 (2008) \emph{Evaluation of measurement data - Guide to the expression
of uncertainty in measurement}. \url{http://www.bipm.org/utils/common/documents/jcgm/JCGM_100_2008_E.pdf}. 
(JCGM 100:2008 is a public domain copy of ISO/IEC \emph{Guide to the expression
of uncertainty in measurement} (1995) ). 
}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
  \code{\link{uncert-class}},  \code{\link{uncert}}.
}
\examples{
  #Example with negative correlation
  x <- list(a=1, b=3, c=2, d=11)
  u <- lapply(x, function(x) x/10)
  u.cor<-diag(1,4)
  u.cor[3,4]<-u.cor[4,3]<- -0.5
  u.form.c<-uncert(~a+b*2+c*3+d/2, x, u, method="NUM", cor=u.cor)

  contribs(u.form.c, simplify=FALSE)
  contribs(u.form.c)
  contribs(u.form.c, as.sd=TRUE)
  contribs(u.form.c, as.sd=TRUE, keep.sign=FALSE)

  contribs(u.form.c, scope=c("a", "c", "d") )

  #Effects of formula specification for scope:
  contribs(u.form.c, ~.)           #All contributions
  contribs(u.form.c, ~(a+b+c+d)^2) #same as ~.
  contribs(u.form.c, ~a+b+c+d )    #single-variable contributions only
  contribs(u.form.c, ~., expand.dot=FALSE )    # as ~a+b+c+d
  contribs(u.form.c, ~.-d)         #Drops d and c:d
  contribs(u.form.c, ~.-c:d)
  contribs(u.form.c, ~c+d)
  contribs(u.form.c, ~c*d)

}

\keyword{ univar }

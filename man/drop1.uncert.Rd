\name{drop1.uncert}
\Rdversion{1.1}
\alias{drop1.uncert}
\alias{drop1.uncertMC}
\alias{print.drop1.uncert}
\alias{plot.drop1.uncert}
\title{
Single variable deletions from the uncertainty budget.
}
\description{
\code{drop1} calculates revised combined uncertainty for single variable deletions from 
an object of class 'uncert'.
}
\usage{
\method{drop1}{uncert}(object, scope, simplify = TRUE, 
	which=c("\% Change", "var", "u", "var.change", "u.change"), \dots)

\method{drop1}{uncertMC}(object, scope, simplify = TRUE, 
	which=c("\% Change", "var", "u", "var.change", "u.change"), \dots)

#Print and plot methods 
\method{print}{drop1.uncert}(x, \dots, digits=2)

\method{plot}{drop1.uncert}(x, \dots, 
	which=c("\% Change", "var", "u", "var.change", "u.change"))
}

\arguments{
  \item{object}{An object of class \sQuote{uncert} or \sQuote{uncertMC}.}
  \item{scope}{character vector, expression or formula containing the list of 
    variables to be dropped. If missing, all variables in object$budget are taken
    as \code{}scope.}
  \item{simplify}{logical. If \code{simplify=TRUE}, the return value is simplified to 
  a named vector. If \code{FALSE}, all forms available (see \code{which}) are 
  returned in a data frame.}
  \item{which}{ logical; controls the form of information returned when \code{simplify=TRUE}. 
  Possible values are:
  \describe{
     \item{\code{"var"}}{The modified values of \eqn{u(y)^2}.}
     \item{\code{"u"}}{The modified values of \eqn{u(y)}.}
     \item{\code{"var.change"}}{The signed changes in \eqn{u(y)^2}.}
     \item{\code{"u.change"}}{The signed changes in \eqn{u(y)}.}
     \item{\code{"\% Change"}}{The percentage change in \eqn{u(y)}.}
  }
 }
  \item{x}{An object of class \sQuote{drop1.uncert} returned by \code{drop1.uncert}.
}
  \item{\dots}{Further objects passed to other functions.
}
  \item{digits}{number of digits used to format the output. See the \code{digits} argument
  of \code{\link{format}}.
}
}
\details{
By analogy with \code{drop1}, \code{drop1.uncert} perfoms single variable deletions from 
the uncertainty budget in \code{object}, calculates the resulting uncertainty and returns the 
results in the form requested by \code{simplify} and \code{which}.

\sQuote{Single variable deletion} of a variable \eqn{x_i}{x[i]} is equivalent to setting the uncertainty 
\eqn{u(x_i)}{u(x[i])} to zero. Note that this also sets covariance terms involving \eqn{x_i}{x[i]} 
to zero. \code{drop1.uncert} does not support the deletion of single \emph{terms} such as \eqn{cov(i, j)}.

In the case of \sQuote{uncertMC} objects, \code{drop1} currently requires \code{object$MC$x} to be 
present (i.e. \code{uncertMC} called with \code{keep.x=TRUE}). The \code{uncertMC} 
method does not support correlation.

For \code{which="var.change"}, \code{which="u.change"} and \code{which="\% Change"} the 
change on dropping a variable is negative if the uncertainty reduces on removing the variable.

The print method simply prints the output with a header formed from the \code{expr} attribute 
and with '\%' appended to the "\% Change" column.

The plot method produces a barplot of the chosen data column.  A plot for each value in \code{which}
is produced. Arguments in \sQuote{\ldots} are passed to barplot. If not already present in \sQuote{\ldots} 
a default main title and ylab are used. The \code{expr} attribute is shown as marginal text if not NA.
}

\value{
If \code{simplify=FALSE}, an object of class \sQuote{drop1.uncert}, consisting of a 
data frame with row names corresponding to \code{row.names(object$budget)}, columns 
corresponding to all possible values of \code{which} in the order \code{"var", "u", "var.change",
"u.change", "\% Change"}, and an attribute \code{expr} containing a copy of the \code{expr} 
value of the 'uncert' object to which \code{drop1.uncert} is applied.

If \code{simplify=TRUE}, the column of the above data frame corresponding to \code{which}
is returned as a vector with names \code{row.names(object$budget)}.
}
\references{
None.
}
\author{
S. L. R. Ellison, \email{s.ellison@lgcgroup.com}
}

\seealso{
  \code{\link{uncert}},  \code{\link{uncert-class}},  \code{\link{format}} for \code{digits}, 
  \code{\link{barplot}} for available plot parameters. 
}
\examples{
  \dontshow{
  #An example with negative correlation
  x <- list(a=1, b=3, c=2, d=11)
  u <- lapply(x, function(x) x/10)
  u.cor<-diag(1,4)
  u.cor[3,4]<-u.cor[4,3]<- -0.5
  u.form.c<-uncert(~a+b*2+c*3+d/2, x, u, method="NUM", cor=u.cor)
  }

  #Continuing the example from plot.uncert:
  require(graphics)
  
  d1<-drop1(u.form.c, simplify=FALSE)
  d1
  
  plot(d1)
  
  drop1(u.form.c)         #\% change only
  
  

}

\keyword{ univar }

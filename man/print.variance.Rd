% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variance.R
\name{print.variance}
\alias{print.variance}
\title{Special print method for variance class.}
\usage{
\method{print}{variance}(x, ...)
}
\arguments{
\item{x}{\code{"variance"} class object to be printed.}

\item{...}{further arguments to be passed to the next method.}
}
\value{
invisible(x).
}
\description{
Makes objects of class \code{"variance"} be printed differently.
}
\examples{
a <- as.variance(1)
print(a)
a <- annualize(a)
print(a)

}
\seealso{
\code{\link{print.volatility}}.
}

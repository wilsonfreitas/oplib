% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiroot.R
\name{print.multiroot}
\alias{print.multiroot}
\title{Print Bisection Method Result.}
\usage{
\method{print}{multiroot}(x, ...)
}
\arguments{
\item{x}{a list.}

\item{...}{additional arguments, will be ignored.}
}
\value{
\code{invisible(x)}.
}
\description{
Specific method for printing the multiroot class, created for multiroot
function, in a pleasing form.
}
\details{
The multiroot class is a list with slots root, iter, err.
}
\examples{
mt <- structure(list(root = 1, iter = 20, err = 1e-8), class = "multiroot")
print(mt)

}
\seealso{
\code{\link{multiroot}} for vectorial bisection method.
}

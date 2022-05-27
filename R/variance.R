#' variance Class Coercing
#'
#' Generic function for coercing to \code{"variance"}. There is a default
#' method and a method for \code{"volatility"}.
#'
#' @param x an object whose class will determine the method to be dispatched.
#' Value is the volatility or variance.
#' @param ... further arguments to be passed to the next method, such as:
#' @param annualized for default method, TRUE if variance value in x is annualized,
#' FALSE if variance is daily.
#' @param dib days in one year, day count convention.
#'
#' @return
#' All values of x as \code{"variance"} class, according to given daily or
#' annualized type.
#'
#' @seealso
#' \code{\link{as.volatility}}.
#'
#' @examples
#' as.variance(1)
#' as.variance(as.volatility(.5))
#' as.variance(0.2, TRUE, 252)
#' as.variance(as.volatility(100:300 / 100, TRUE))
#'
#' @name as.variance
#' @export
as.variance <- function(x, ...) UseMethod("as.variance")

#' @rdname as.variance
#' @export
as.variance.default <- function(x, annualized = FALSE, dib = NULL, ...) {
  structure(x, annualized = annualized, dib = dib, class = "variance")
}

#' @rdname as.variance
#' @export
as.variance.volatility <- function(x, ...) {
  y <- unclass(x)
  attributes(y) <- NULL
  p <- if (attr(x, "annualized")) 100 else 1
  structure((y / p)^2,
    annualized = attr(x, "annualized"), dib = attr(x, "dib"),
    class = "variance"
  )
}

#' Special print method for variance class.
#'
#' Makes objects of class \code{"variance"} be printed differently.
#'
#' @param x \code{"variance"} class object to be printed.
#' @param ... further arguments to be passed to the next method.
#'
#' @return
#' invisible(x).
#'
#' @seealso
#' \code{\link{print.volatility}}.
#'
#' @examples
#' a <- as.variance(1)
#' print(a)
#' a <- annualize(a)
#' print(a)
#'
#' @export
print.variance <- function(x, ...) {
  cat(
    if (attr(x, "annualized")) "Annual" else "Daily",
    "Variance",
    if (attr(x, "annualized")) paste(attr(x, "dib"), "days"),
    "\n"
  )
  y <- unclass(x)
  attributes(y) <- NULL
  print(y)
  invisible(x)
}
#' volatility Class Coercing
#'
#' Generic function for coercing to \code{"volatility"}. There is a default
#' method and a method for \code{"variance"}.
#'
#' @param x an object whose class will determine the method to be dispatched.
#' Value is the volatility or variance.
#' @param ... further arguments to be passed to the next method, such as:
#' @param annualized for default method, TRUE if volatility value in x is annualized,
#' FALSE if volatility is daily.
#' @param dib days in one year, day count convention.
#'
#' @return
#' All values of x as volatility class, according to given daily or annualized type.
#'
#' @seealso
#' \code{\link{as.variance}}.
#'
#' @examples
#' as.volatility(1)
#' as.volatility(as.variance(.5))
#' as.volatility(0.2, TRUE, 252)
#' as.volatility(as.variance(10:30 / 100, TRUE))
#'
#' @name as.volatility
#' @export
as.volatility <- function(x, ...) UseMethod("as.volatility")

#' @rdname as.volatility
#' @export
as.volatility.default <- function(x, annualized = FALSE, dib = NULL, ...) {
  structure(x, annualized = annualized, dib = dib, class = "volatility")
}

#' @rdname as.volatility
#' @export
as.volatility.variance <- function(x, ...) {
  y <- unclass(x)
  attributes(y) <- NULL
  p <- if (attr(x, "annualized")) 100 else 1
  structure(sqrt(y) * p,
    annualized = attr(x, "annualized"), dib = attr(x, "dib"),
    class = "volatility"
  )
}

#' Special print method for volatility class.
#'
#' Makes objects of class \code{"volatility"} be printed differently.
#'
#' @param x \code{"volatility"} class object to be printed.
#' @param ... further arguments to be passed to the next method.
#'
#' @return
#' invisible(x).
#'
#' @seealso
#' \code{\link{print.variance}}.
#'
#' @examples
#' a <- as.volatility(1)
#' print(a)
#' a <- annualize(a)
#' print(a)
#'
#' @export
print.volatility <- function(x, ...) {
  cat(
    if (attr(x, "annualized")) "Annual" else "Daily",
    "Volatility",
    if (attr(x, "annualized")) paste("(%)", attr(x, "dib"), "days"),
    "\n"
  )
  y <- unclass(x)
  attributes(y) <- NULL
  print(y)
  invisible(x)
}
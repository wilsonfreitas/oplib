#' Annualize Volatility or Variance
#'
#' Transforms variance or volatility to annual. There is no default
#' method, so only classes \code{"volatility"}, \code{"variance"}, are supported.
#'
#' @param x an object whose class will determine the method to be dispatched.
#' Volatility or variance, will become annual.
#' @param ... further arguments to be passed to the next method.
#' @param dib days in one year, day count convention.
#'
#' @return
#' Same class with annual volatility or variance.
#'
#' @seealso
#' \code{\link{daily}} to transform variance or volatility to daily basis.
#'
#' @examples
#' annualize(as.variance(.5))
#' annualize(as.volatility(.25))
#'
#' @name annualize
#' @export
annualize <- function(x, ...) UseMethod("annualize")

#' @rdname annualize
#' @export
annualize.volatility <- function(x, dib = 252, ...) {
  y <- unclass(x)
  attributes(y) <- NULL
  if (attr(x, "annualized")) {
    return(x)
  } else {
    y <- y * sqrt(dib) * 100
  }
  as.volatility(y, annualized = TRUE, dib = dib)
}

#' @rdname annualize
#' @export
annualize.variance <- function(x, dib = 252, ...) {
  y <- unclass(x)
  attributes(y) <- NULL
  if (attr(x, "annualized")) {
    return(x)
  } else {
    y <- y * dib
  }
  as.variance(y, annualized = TRUE, dib = dib)
}

#' Daily Volatility or Variance
#'
#' Transforms variance or volatility to daily. There is no default
#' method, so only classes \code{"volatility"}, \code{"variance"}, are supported.
#'
#' @param x an object whose class will determine the method to be dispatched.
#' Volatility or variance, will become daily.
#' @param ... further arguments to be passed to the next method.
#'
#' @return
#' Same class with daily volatility or variance.
#'
#' @seealso
#' \code{\link{annualize}} to transform variance or volatility to annual basis.
#'
#' @examples
#' daily(annualize(as.variance(.5)))
#' daily(annualize(as.volatility(.25)))
#'
#' @name daily
#' @export
daily <- function(x, ...) UseMethod("daily")

#' @rdname daily
#' @export
daily.volatility <- function(x, ...) {
  y <- unclass(x)
  attributes(y) <- NULL
  if (!attr(x, "annualized")) {
    return(x)
  } else {
    y <- y / (sqrt(attr(x, "dib")) * 100)
  }
  as.volatility(y, annualized = FALSE, dib = NULL)
}

#' @rdname daily
#' @export
daily.variance <- function(x, ...) {
  y <- unclass(x)
  attributes(y) <- NULL
  if (!attr(x, "annualized")) {
    return(x)
  } else {
    y <- y / attr(x, "dib")
  }
  as.variance(y, annualized = FALSE, dib = NULL)
}
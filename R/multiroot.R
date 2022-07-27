#' Vectorial bisection method varying parameters from one function.
#'
#' Searches the interval for a root (i.e., zero) of the function func_ with
#' respect to its first argument. Other arguments may be vectorial, so
#' \code{multiroot} will return a vector with roots.
#'
#' \code{multiroot} was created aiming implied volatility calculation for
#' several options.
#'
#' @param func_ objective function for root finding.
#' @param interval a vector containing the end-points of the interval to be
#' searched for all the roots.
#' @param ... additional named or unnamed arguments to be passed to func_,
#' may be vectorial, and this is the difference between multiroot and uniroot.
#' @param tolerance the desired accuracy (convergence tolerance).
#' @param maxiter the maximum number of iterations.
#' @param check function for checking convergence in \code{\link{any}} or
#' \code{\link{all}} of the roots. Since it is bisection method, they are
#' the same.
#'
#' @return
#' Vector with roots.
#'
#' @seealso
#' \code{\link{uniroot}} for R built-in one dimensional root finding.
#'
#' @examples
#' f <- function(x, y) x - y
#' multiroot(f, c(-3, 10), 1:3)
#' g <- function(sigma, premium, type, spot, strike, time, rate, yield) {
#'   bsmprice(type, spot, strike, time, rate, yield, sigma) - premium
#' }
#' multiroot(g, c(1e-8, 10), 5, "call", 30, 31, c(1, 1.1, 1.2, 1.3), 0.15, 0)
#' @export
multiroot <- function(func_, interval, ...,
                      tolerance = .Machine$double.eps,
                      maxiter = 100,
                      check = function(...) any(..., na.rm = TRUE)) {
  s_lower <- sign(func_(interval[1], ...))
  s_upper <- sign(func_(interval[2], ...))
  ix <- s_upper * s_lower != -1

  func <- function(a, ...) func_(a, ...) * s_upper
  initial_guess <- sum(interval) / 2
  x <- initial_guess * rep(1, length(s_upper))
  lower <- interval[1] * rep(1, length(s_upper))
  upper <- interval[2] * rep(1, length(s_upper))

  iter <- 0
  err <- func(x, ...)
  err[ix] <- NA
  stopifnot(!all(is.na(err)))
  while (check(abs(err) > tolerance) && iter < maxiter) {
    idx.err <- err < 0
    idx.err <- ifelse(is.na(idx.err), FALSE, idx.err)
    lower[idx.err] <- x[idx.err]
    x[idx.err] <- (upper[idx.err] + x[idx.err]) / 2
    idx.err <- !idx.err
    upper[idx.err] <- x[idx.err]
    x[idx.err] <- (lower[idx.err] + x[idx.err]) / 2
    err <- func(x, ...)
    err[ix] <- NA
    x[ix] <- NA
    iter <- iter + 1
    stopifnot(!all(is.na(err)))
  }
  res <- list(root = x, iter = iter, err = err)

  structure(res, class = "multiroot")
}

#' Print Bisection Method Result.
#'
#' Specific method for printing the multiroot class, created for multiroot
#' function, in a pleasing form.
#'
#' The multiroot class is a list with slots root, iter, err.
#'
#' @param x a list.
#' @param ... additional arguments, will be ignored.
#'
#' @return
#' \code{invisible(x)}.
#'
#' @seealso
#' \code{\link{multiroot}} for vectorial bisection method.
#'
#' @examples
#' mt <- structure(list(root = 1, iter = 20, err = 1e-8), class = "multiroot")
#' print(mt)
#'
#' @export
print.multiroot <- function(x, ...) {
  cat("root:", x$root, "\n")
  cat("iter:", x$iter, "\n")
  cat("err:", x$err, "\n")
  invisible(x)
}
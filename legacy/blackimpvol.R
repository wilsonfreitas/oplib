#' Computes implied volatility inverting Black-76 formula.
#'
#' Uses bisection method to calculate implied volatility, given other parameters
#' for Black-76 formula for pricing options on futures.
#'
#' @param option_prices future option price.
#' @param type 'call' for call option, any other value for put option.
#' @param future current future price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param tolerance the tolerance for the approximation.
#' @param maxiter maximum number of iterations for bisection method.
#' @param check function for checking the convergence of all bisections.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with approximations implied volatility. If the algorithm converges
#' before \code{maxiter} is reached, Black-76 formula calculated with this
#' volatility should not differ from \code{option_prices} by more than
#' \code{tolerance}.
#'
#' @seealso
#' \code{\link{bsmimpvol}} for options on stocks, with Black-Scholes-Merton model.
#'
#' @examples
#' blackimpvol(c(3.66765, 5.441491), c("call", "put"), 50, 52, 1, 0.12)
#' blackimpvol(c(2, 3, 4, 5, 6, 7), "put", 50, 52, 1, 0.12)
#' blackimpvol(2000:6000 / 1000, "call", 50, 52, 1, 0.12)
#'
#' @export
blackimpvol <- function(option_prices, type, future, strike, time, rate,
                        tolerance = .Machine$double.eps, maxiter = 100, check = any) {
  f <- function(sigma, ...) {
    blackprice(sigma = sigma, ...) - option_prices
  }

  res <- multiroot(
    func_ = f, interval = c(1e-8, 10), tolerance = tolerance,
    maxiter = maxiter, check = check, type = type, future = future,
    strike = strike, rate = rate, time = time
  )
  res$root
}
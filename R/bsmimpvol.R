#' Computes implied volatility inverting Black-Scholes-Merton formula.
#'
#' Uses bisection method to calculate implied volatility, given other parameters
#' for Black-Scholes-Merton formula for pricing options on stocks.
#'
#' @param option_prices stock option price.
#' @param type 'call' for call option, any other value for put.
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
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
#' An approximation for the implied volatility. If the algorithm converges before
#' \code{maxiter} is reached, Black-Scholes-Merton formula calculated with this
#' volatility should not differ from \code{option_prices} by more than
#' \code{tolerance}.
#'
#' @examples
#' bsmimpvol(6, "call", 50, 52, 1, 0.1, 0)
#' bsmimpvol(c(3, 4, 5, 6), "call", 50, 52, 1, 0.1, 0)
#' bsmimpvol(1.5, "put", 40, 38, 0.5, 0.15, 0, tolerance = 1e-8)
#'
#' @export
bsmimpvol <- function(option_prices,
                      type, spot, strike, time, rate, yield,
                      tolerance = .Machine$double.eps,
                      maxiter = 100,
                      check = any) {
  f <- function(sigma, ...) {
    bsmprice(sigma = sigma, ...) - option_prices
  }

  res <- multiroot(
    func_ = f, interval = c(1e-8, 10), tolerance = tolerance,
    maxiter = maxiter, check = check, type = type, spot = spot,
    strike = strike, rate = rate, time = time, yield = yield
  )
  res$root
}
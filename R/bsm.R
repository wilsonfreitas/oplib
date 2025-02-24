#' Computes option price via Black-Scholes-Merton formula.
#'
#' Black-Scholes-Merton formula implementation, pricing european options.
#'
#' @param type 'call' for call option, any other value for put.
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param sigma volatility of the stock price.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with price of european options on stocks, according to
#' Black-Scholes-Merton model.
#'
#' @examples
#' bsmprice(c("call", "put"), c(50, 49, 48, 47), 45, 0.25, 0.13, 0.01, 0.2)
#' bsmprice("call", 50, 45, 0.25, 0.13, 0.01, 1:1000 / 1000)
#' bsmprice("put", 20, 30, 0:2000 / 1000, 0.15, 0, 0.25)
#'
#' @importFrom stats pnorm
#' @export
bsmprice <- function(type, spot, strike, time, rate, yield, sigma) {
  df <- data.frame(type, spot, strike, time, rate, yield, sigma)
  with(df, {
    d1 <- bsmd1(spot, strike, time, rate, yield, sigma)
    d2 <- d1 - sigma * sqrt(time)

    ifelse(
      tolower(as.character(type)) == "call",
      spot * exp((-yield) * time) *
        pnorm(d1) - strike * exp(-rate * time) * pnorm(d2),
      strike * exp(-rate * time) *
        pnorm(-d2) - spot * exp((-yield) * time) * pnorm(-d1)
    )
  })
}

#' Computes d1 from Black-Scholes-Merton model.
#'
#' d1 formula implementation, from Black-Scholes-Merton model.
#'
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param sigma volatility of the stock price.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with d1.
#'
#' @examples
#' bsmd1(c(50, 49, 48, 47), 45, 0.25, 0.13, 0.01, 0.2)
#' bsmd1(50, 45, 0.25, 0.13, 0.01, 1:1000 / 1000)
#' bsmd1(20, 30, 0:2000 / 1000, 0.15, 0, 0.25)
#'
#' @export
bsmd1 <- function(spot, strike, time, rate, yield, sigma) {
  (log(spot / strike) + (rate - yield + sigma * sigma / 2) * time) /
    (sigma * sqrt(time))
}

#' Computes the desired greek letter from Black-Scholes-Merton model.
#'
#' Greeks formula implementation, from Black-Scholes-Merton model.
#'
#' @param greeks desired greek letter: delta, gamma, vega, rho or theta.
#' @param type call for call option, any other value for put.
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param sigma volatility of the stock price.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with desired greeks for input data.
#'
#' @examples
#' bsmgreeks("delta", "put", 20, 30, 0:2000 / 1000, 0.15, 0, 0.25)
#'
#' @export
bsmgreeks <- function(greeks = c("delta", "gamma", "vega", "rho", "theta"),
                      type, spot, strike, time, rate, yield, sigma) {
  df <- data.frame(type, spot, strike, time, rate, yield, sigma)
  grks <- match.arg(greeks)
  with(df, {
    if (grks == "delta") {
      bsmdelta(type, spot, strike, time, rate, yield, sigma)
    } else if (grks == "gamma") {
      bsmgamma(type, spot, strike, time, rate, yield, sigma)
    } else if (grks == "vega") {
      bsmvega(type, spot, strike, time, rate, yield, sigma)
    } else if (grks == "rho") {
      bsmrho(type, spot, strike, time, rate, yield, sigma)
    } else if (grks == "theta") {
      bsmtheta(type, spot, strike, time, rate, yield, sigma)
    } else {
      stop("Unknown greek:", grks)
    }
  })
}

#' @rdname bsmgreeks
#'
#' @section bsmdelta:
#' Delta greek letter measures the rate of change of the theoretical option
#' value with respect to changes in the underlying asset's price.
#'
#' @examples
#' bsmdelta(c("call", "put"), c(50, 49, 48, 47), 45, 0.25, 0.13, 0.01, 0.2)
#' bsmdelta("call", 50, 45, 0.25, 0.13, 0.01, 1:1000 / 1000)
#' bsmdelta("put", 20, 30, 0:2000 / 1000, 0.15, 0, 0.25)
#'
#' @export
bsmdelta <- function(type, spot, strike, time, rate, yield, sigma) {
  data.frame(type, spot, strike, time, rate, yield, sigma) |> with({
    d1 <- bsmd1(spot, strike, time, rate, yield, sigma)
    ifelse(
      tolower(as.character(type)) == "call",
      exp((-yield) * time) * pnorm(d1),
      exp((-yield) * time) * (pnorm(d1) - 1)
    )
  })
}

#' @rdname bsmgreeks
#'
#' @section bsmvega:
#' Vega greek letter measures sensitivity of volatility.
#'
#' @examples
#' bsmvega(c("call", "put"), c(50, 49, 48, 47), 45, 0.25, 0.13, 0.01, 0.2)
#' bsmvega("call", 50, 45, 0.25, 0.13, 0.01, 1:1000 / 1000)
#' bsmvega("put", 20, 30, 0:2000 / 1000, 0.15, 0, 0.25)
#'
#' @importFrom stats dnorm
#' @export
bsmvega <- function(type, spot, strike, time, rate, yield, sigma) {
  data.frame(type, spot, strike, time, rate, yield, sigma) |> with({
    d1 <- bsmd1(spot, strike, time, rate, yield, sigma)
    spot * exp((-yield) * time) * dnorm(d1) * sqrt(time)
  })
}

#' @rdname bsmgreeks
#'
#' @section bsmgamma:
#' Gamma greek letter measures the rate of change in the delta with respect
#' to changes in the underlying price.
#'
#' @examples
#' bsmgamma(c("call", "put"), c(50, 49, 48, 47), 45, 0.25, 0.13, 0.01, 0.2)
#' bsmgamma("call", 50, 45, 0.25, 0.13, 0.01, 1:1000 / 1000)
#' bsmgamma("put", 20, 30, 0:2000 / 1000, 0.15, 0, 0.25)
#'
#' @export
bsmgamma <- function(type, spot, strike, time, rate, yield, sigma) {
  data.frame(type, spot, strike, time, rate, yield, sigma) |> with({
    d1 <- bsmd1(spot, strike, time, rate, yield, sigma)
    exp((-yield) * time) * dnorm(d1) / (spot * sigma * sqrt(time))
  })
}

#' @rdname bsmgreeks
#'
#' @section bsmtheta:
#' Theta measures the sensitivity of the value of the derivative to the
#' passage of time.
#
#' @examples
#' bsmtheta(c("call", "put"), c(50, 49, 48, 47), 45, 0.25, 0.13, 0.01, 0.2)
#' bsmtheta("call", 50, 45, 0.25, 0.13, 0.01, 1:1000 / 1000)
#' bsmtheta("put", 20, 30, 0:2000 / 1000, 0.15, 0, 0.25)
#'
#' @export
bsmtheta <- function(type, spot, strike, time, rate, yield, sigma) {
  data.frame(type, spot, strike, time, rate, yield, sigma) |> with({
    d1 <- bsmd1(spot, strike, time, rate, yield, sigma)
    d2 <- d1 - sigma * sqrt(time)
    Theta1 <- -(spot * exp((-yield) * time) * dnorm(d1) * sigma) /
      (2 * sqrt(time))
    ifelse(
      tolower(as.character(type)) == "call",
      Theta1 - (-yield) * spot * exp((-yield) * time) * pnorm(+d1) -
        rate * strike * exp(-rate * time) * pnorm(+d2),
      Theta1 + (-yield) * spot * exp((-yield) * time) * pnorm(-d1) +
        rate * strike * exp(-rate * time) * pnorm(-d2)
    )
  })
}

#' @rdname bsmgreeks
#'
#' @section bsmrho:
#' Rho measures sensitivity to the interest rate: it is the derivative
#' of the option value with respect to the risk free interest rate.
#'
#' @examples
#' bsmrho(c("call", "put"), c(50, 49, 48, 47), 45, 0.25, 0.13, 0.01, 0.2)
#' bsmrho("call", 50, 45, 0.25, 0.13, 0.01, 1:1000 / 1000)
#' bsmrho("put", 20, 30, 0:2000 / 1000, 0.15, 0, 0.25)
#'
#' @export
bsmrho <- function(type, spot, strike, time, rate, yield, sigma) {
  data.frame(type, spot, strike, time, rate, yield, sigma) |> with({
    d1 <- bsmd1(spot, strike, time, rate, yield, sigma)
    d2 <- d1 - sigma * sqrt(time)
    CallPut <- bsmprice(type, spot, strike, time, rate, yield, sigma)
    ifelse(
      tolower(as.character(type)) == "call",
      ifelse(
        rate - yield != 0,
        time * strike * exp(-rate * time) * pnorm(d2),
        -time * CallPut
      ),
      ifelse(
        rate - yield != 0,
        -time * strike * exp(-rate * time) * pnorm(-d2),
        -time * CallPut
      )
    )
  })
}

#' Computes a measure of moneyness for Black-Scholes-Merton model.
#'
#' Moneyness formula implementation, for Black-Scholes-Merton model.
#'
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with moneyness.
#'
#' @examples
#' bsmmoneyness(0:49, 1:50, 3, 0.1, 0)
#' bsmmoneyness(10, 10, 1, 0.12, 0:1000 / 1000)
#' bsmmoneyness(10, 10, 1, 0:1000 / 1000, 0)
#'
#' @export
bsmmoneyness <- function(spot, strike, time, rate, yield) {
  1 - spot * exp(-yield * time) / (strike * exp(-rate * time))
}

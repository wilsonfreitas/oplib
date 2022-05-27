#' Computes option price via Black-76 formula.
#'
#' Black-76 formula implementation, pricing european options.
#'
#' @param type 'call' for call option, any other value for put.
#' @param future current future price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param sigma volatility of the future price.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with price of european options on futures, according to Black-76 model.
#'
#' @seealso
#' \code{\link{bsmprice}} for pricing options on stocks, with
#' Black-Scholes-Merton model.
#'
#' @examples
#' blackprice(c("put", "call"), c(50, 49, 48, 47), 45, 0.25, 0.13, 0.3)
#' blackprice("p", 100:500 / 10, 45, 0.25, 0.13, 0.17)
#' blackprice("call", 50, 100:900 / 10, 2, 0.1, 0.18)
#'
#' @export
blackprice <- function(type, future, strike, time, rate, sigma) {
  df <- data.frame(type, future, strike, time, rate, sigma)
  with(df, {
    d1 <- (log(future / strike) + (sigma * sigma / 2) * time) / (sigma * sqrt(time))
    d2 <- d1 - sigma * sqrt(time)
    ifelse(tolower(as.character(type)) == "call",
      {
        future * exp((-rate) * time) * pnorm(d1) - strike * exp(-rate * time) * pnorm(d2)
      },
      {
        strike * exp(-rate * time) * pnorm(-d2) - future * exp((rate) * time) * pnorm(-d1)
      }
    )
  })
}

#' Computes the desired greek letter from Black-76 model.
#'
#' Greeks formula implementation, from Black-76 model.
#'
#' @param greeks desired greek letter: 'delta', 'gamma', 'vega', 'rho' or 'theta'.
#' @param type 'call' for call option, any other value for put.
#' @param future current underlying forward price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
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
#' blackgreeks("delta", "put", 20, 30, 0:2000 / 1000, 0.15, 0.25)
#'
#' @export
blackgreeks <- function(greeks = c("delta", "gamma", "vega", "rho", "theta"),
                        type, future, strike, time, rate, sigma) {
  df <- data.frame(type, future, strike, time, rate, sigma)
  greeks <- match.arg(greeks)
  with(df, {
    if (greeks == "delta") blackdelta(type, future, strike, time, rate, sigma)
    if (greeks == "gamma") blackgamma(type, future, strike, time, rate, sigma)
    if (greeks == "vega") blackvega(type, future, strike, time, rate, sigma)
    if (greeks == "rho") blackrho(type, future, strike, time, rate, sigma)
    if (greeks == "theta") blacktheta(type, future, strike, time, rate, sigma)
  })
}

#' @rdname blackgreeks
#'
#' @section blackdelta:
#' Delta measures the rate of change of the theoretical option value with respect to changes in the underlying asset's price.
#'
#' @seealso
#' \code{\link{bsmdelta}} for delta from Black-Scholes-Merton model.
#'
#' @examples
#' blackdelta(c("call", "put"), c(50, 49, 48, 47), c(45, 40), 0.25, 0.13, 0.2)
#' blackdelta("call", 50, 45, 0.25, 0.13, 1:1000 / 1000)
#' blackdelta("put", 20, 30, 0:2000 / 1000, 0.15, 0.25)
#'
#' @export
blackdelta <- function(type, future, strike, time, rate, sigma) {
  d1 <- (log(future / strike) + (sigma * sigma / 2) * time) / (sigma * sqrt(time))
  if (tolower(as.character(type)) == "call") {
    result <- exp((-rate) * time) * pnorm(d1)
  } else {
    result <- exp((-rate) * time) * (pnorm(d1) - 1)
  }
  result
}

#' @rdname blackgreeks
#'
#' @section blackvega:
#' Vega measures sensitivity of volatility.
#'
#' @seealso
#' \code{\link{bsmvega}} for vega from Black-Scholes-Merton model.
#'
#' @examples
#' blackvega(c("call", "put"), c(50, 49, 48, 47), c(45, 40), 0.25, 0.13, 0.2)
#' blackvega("call", 50, 45, 0.25, 0.13, 1:1000 / 1000)
#' blackvega("put", 20, 30, 0:2000 / 1000, 0.15, 0.25)
#'
#' @export
blackvega <- function(type, future, strike, time, rate, sigma) {
  d1 <- (log(future / strike) + (sigma * sigma / 2) * time) / (sigma * sqrt(time))
  result <- future * exp((-rate) * time) * dnorm(d1) * sqrt(time)
  result
}

#' @rdname blackgreeks
#'
#' @section blackgamma:
#' Gamma measures the rate of change in the delta with respect to changes in the underlying price.
#'
#' @seealso
#' \code{\link{bsmgamma}} for gamma from Black-Scholes-Merton model.
#'
#' @examples
#' blackgamma(c("call", "put"), c(50, 49, 48, 47), c(45, 40), 0.25, 0.13, 0.2)
#' blackgamma("call", 50, 45, 0.25, 0.13, 1:1000 / 1000)
#' blackgamma("put", 20, 30, 0:2000 / 1000, 0.15, 0.25)
#'
#' @export
blackgamma <- function(type, future, strike, time, rate, sigma) {
  d1 <- (log(future / strike) + (sigma * sigma / 2) * time) / (sigma * sqrt(time))
  result <- exp((-rate) * time) * dnorm(d1) / (future * sigma * sqrt(time))
  result
}

#' @rdname blackgreeks
#'
#' @section blacktheta:
#' Theta measures the sensitivity of the value of the derivative to the passage of time.
#'
#' @seealso
#' \code{\link{bsmtheta}} for theta from Black-Scholes-Merton model.
#'
#' @examples
#' blacktheta(c("call", "put"), c(50, 49, 48, 47), c(45, 40), 0.25, 0.13, 0.2)
#' blacktheta("call", 50, 45, 0.25, 0.13, 1:1000 / 1000)
#' blacktheta("put", 20, 30, 0:2000 / 1000, 0.15, 0.25)
#'
#' @export
blacktheta <- function(type, future, strike, time, rate, sigma) {
  d1 <- (log(future / strike) + (sigma * sigma / 2) * time) / (sigma * sqrt(time))
  d2 <- d1 - sigma * sqrt(time)
  Theta1 <- -(future * exp((-rate) * time) * dnorm(d1) * sigma) / (2 * sqrt(time))
  if (tolower(as.character(type)) == "call") {
    result <- Theta1 - (-rate) * future * exp((-rate) * time) * pnorm(+d1) - rate * strike * exp(-rate * time) * pnorm(+d2)
  } else {
    result <- Theta1 + (-rate) * future * exp((-rate) * time) * pnorm(-d1) + rate * strike * exp(-rate * time) * pnorm(-d2)
  }
  result
}

#' @rdname blackgreeks
#'
#' @section blackrho:
#' Rho measures sensitivity to the interest rate: it is the derivative of the option value with respect to the risk free interest rate.
#'
#' @seealso
#' \code{\link{bsmrho}} for rho from Black-Scholes-Merton model.
#'
#' @examples
#' blackrho(c("call", "put"), c(50, 49, 48, 47), c(45, 40), 0.25, 0.13, 0.2)
#' blackrho("call", 50, 45, 0.25, 0.13, 1:1000 / 1000)
#' blackrho("put", 20, 30, 0:2000 / 1000, 0.15, 0.25)
#'
#' @export
blackrho <- function(type, future, strike, time, rate, sigma) {
  d1 <- (log(future / strike) + (sigma * sigma / 2) * time) / (sigma * sqrt(time))
  d2 <- d1 - sigma * sqrt(time)
  CallPut <- bsmprice(type, future, strike, time, rate, rate, sigma)
  result <- -time * CallPut
  result
}

#' Computes d1 from Black-76 model.
#'
#' d1 formula implementation, from Black-76 model.
#'
#' @param future current future price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
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
#' @seealso
#' \code{\link{bsmd1}} for d1 from Black-Scholes-Merton model.
#'
#' @examples
#' blackd1(c(50, 49, 48, 47), 45, 0.25, 0.13, 0.2)
#' blackd1(50, 45, 0.25, 0.13, 1:1000 / 1000)
#' blackd1(20, 30, 0:2000 / 1000, 0.15, 0.25)
#'
#' @export
blackd1 <- function(future, strike, time, rate, sigma) {
  (log(future / strike) + (sigma * sigma / 2) * time) / (sigma * sqrt(time))
}


#' Computes a measure of moneyness for Black-76 model.
#'
#' Moneyness formula implementation, for Black-76 model.
#'
#' @param future current future price.
#' @param strike the strike price.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with moneyness.
#'
#' @seealso
#' \code{\link{bsmmoneyness}} for one measure for moneyness in
#' Black-Scholes-Merton model.
#'
#' @examples
#' blackmoneyness(2:51, 1:50)
#' blackmoneyness(10, 50:150 / 10)
#' blackmoneyness(0:600 / 10, 30)
#'
#' @export
blackmoneyness <- function(future, strike) {
  1 - future / strike
}


blackstrike <- function(future, delta, sigma, t) {
  inv <- qnorm(delta)
  ee <- sigma * sqrt(t)
  strike <- future / (exp((inv - 0.5 * ee) * ee))
  strike
}
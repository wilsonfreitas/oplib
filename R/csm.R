#' Computes w from modified Corrado-Su model.
#'
#' w formula implementation from modified Corrado-Su model.
#'
#' @param sigma volatility of the stock price.
#' @param time time to option expiration in years.
#' @param mu3 skewness.
#' @param mu4 kurtosis.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with w.
#'
#' @seealso
#' \code{\link{csmd}} for d from modified Corrado-Su model.\\
#' \code{\link{csmq3}} for Q_3 from modified Corrado-Su model.\\
#' \code{\link{csmq4}} for Q_4 from modified Corrado-Su model.\\
#' \code{\link{csmprice}} for modified Corrado-Su option pricing formula.
#'
#' @examples
#' csmw(0.25, 1.5, 0, -100:100 / 100)
#' csmw(0.25, 1.5, -100:100 / 100, 3)
#' csmw(0:100 / 100, 1.5, 0, 3)
#'
#' @export
csmw <- function(sigma, time, mu3, mu4) {
  (mu3 / 6) * (sigma^3) * (sqrt(time^3)) + (mu4 / 24) * (sigma^4) * (time^2)
}

#' Computes d from modified Corrado-Su model.
#'
#' d formula implementation from modified Corrado-Su model.
#'
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param sigma volatility of the stock price.
#' @param w constant from modified Corrado-Su model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with d.
#'
#' @seealso
#' \code{\link{csmw}} for w from modified Corrado-Su model.\\
#' \code{\link{csmq3}} for Q_3 from modified Corrado-Su model.\\
#' \code{\link{csmq4}} for Q_4 from modified Corrado-Su model.\\
#' \code{\link{csmprice}} for modified Corrado-Su option pricing formula.
#'
#' @examples
#' csmd(14, 120:160 / 10, 0.5, 0.15, 0.1, 0.18, csmw(0.18, 0.5, 0, 3))
#' csmd(14, 13, 0.5, 0.15, 0.1, 0.18, csmw(0.18, 0.5, -1000:1000 / 100, 3))
#' csmd(14, 13, 0.5, 0.15, 0.1, 0.18, csmw(0.18, 0.5, 0, -1000:1000 / 100))
#'
#' @export
csmd <- function(spot, strike, time, rate, yield, sigma, w) {
  (log(spot / strike) + (rate - yield + 0.5 * sigma^2) * time - log(1 + w)) /
    (sigma * sqrt(time))
}


#' Computes Q_3 from modified Corrado-Su model.
#'
#' Q_3 formula implementation from modified Corrado-Su model.
#'
#' @param spot current stock price.
#' @param sigma volatility of the stock price.
#' @param time time to option expiration in years.
#' @param d1 d value from modified Corrado-Su model.
#' @param w constant from modified Corrado-Su model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with Q_3.
#'
#' @seealso
#' \code{\link{csmw}} for w from modified Corrado-Su model.\\
#' \code{\link{csmd}} for d from modified Corrado-Su model.\\
#' \code{\link{csmq4}} for Q_4 from modified Corrado-Su model.\\
#' \code{\link{csmprice}} for modified Corrado-Su option pricing formula.
#'
#' @examples
#' spot <- 25
#' sigma <- 0.2
#' time <- 0.5
#' strike <- 26
#' rate <- 0.12
#' yield <- 0
#' skewness <- 0
#' kurtosis <- 3
#' w <- csmw(sigma, time, skewness, kurtosis)
#' d <- csmd(spot, strike, time, rate, yield, sigma, w)
#' csmq3(spot, sigma, time, d, w)
#' csmq3(200:320 / 10, sigma, time, d, w)
#'
#' @export
csmq3 <- function(spot, sigma, time, d1, w) {
  q3 <- spot * sigma * sqrt(time) / (6 * (1 + w))
  q3 * (2 * sigma * sqrt(time) - d1) * dnorm(d1)
}


#' Computes Q_4 from modified Corrado-Su model.
#'
#' Q_4 formula implementation from modified Corrado-Su model.
#'
#' @param spot current stock price.
#' @param sigma volatility of the stock price.
#' @param time time to option expiration in years.
#' @param d1 d value from modified Corrado-Su model.
#' @param w constant from modified Corrado-Su model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with Q_4.
#'
#' @seealso
#' \code{\link{csmw}} for w from modified Corrado-Su model.\\
#' \code{\link{csmd}} for d from modified Corrado-Su model.\\
#' \code{\link{csmq3}} for Q_3 from modified Corrado-Su model.\\
#' \code{\link{csmprice}} for modified Corrado-Su option pricing formula.
#'
#' @examples
#' spot <- 25
#' sigma <- 0.2
#' time <- 0.5
#' strike <- 26
#' rate <- 0.12
#' yield <- 0
#' skewness <- 0
#' kurtosis <- 3
#' w <- csmw(sigma, time, skewness, kurtosis)
#' d <- csmd(spot, strike, time, rate, yield, sigma, w)
#' csmq4(spot, sigma, time, d, w)
#' csmq4(200:320 / 10, sigma, time, d, w)
#'
#' @export
csmq4 <- function(spot, sigma, time, d1, w) {
  q4 <- (d1^2 - 3 * d1 * sigma * sqrt(time) + 3 * sigma^2 * time - 1) *
    dnorm(d1)
  q4 <- q4 * spot * sigma * sqrt(time) / (24 * (1 + w))
  return(q4)
}


#' Computes option price with modified Corrado-Su model.
#'
#' Implementation of modified Corrado-Su formula for option pricing.
#'
#' @param type 'call' for call option, any other value for put.
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param sigma volatility of the stock price.
#' @param mu3 skewness.
#' @param mu4 kurtosis.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with option prices.
#'
#' @seealso
#' \code{\link{csmw}} for w from modified Corrado-Su model.\\
#' \code{\link{csmd}} for d from modified Corrado-Su model.\\
#' \code{\link{csmq3}} for Q_3 from modified Corrado-Su model.\\
#' \code{\link{csmq4}} for Q_4 from modified Corrado-Su model.
#'
#' @examples
#' spot <- 25
#' sigma <- 0.2
#' time <- 0.5
#' strike <- 26
#' rate <- 0.12
#' yield <- 0
#' skewness <- 1
#' kurtosis <- 4
#' csmprice(c("call", "put"), spot, strike, time, rate, yield, sigma, 0, 3) ==
#'   bsmprice(c("call", "put"), spot, strike, time, rate, yield, sigma)
#' csmprice(
#'   "put", spot, strike, time, rate, 0:100 / 100, sigma, skewness, kurtosis
#' )
#'
#' @export
csmprice <- function(type, spot, strike, time, rate, yield, sigma, mu3, mu4) {
  c_bs <- bsmprice("call", spot, strike, time, rate, yield, sigma)
  w <- csmw(sigma, time, mu3, mu4)
  d1 <- csmd(spot, strike, time, rate, yield, sigma, w)
  q3 <- csmq3(spot, sigma, time, d1, w)
  q4 <- csmq4(spot, sigma, time, d1, w)
  callp <- c_bs + mu3 * q3 + (mu4 - 3) * q4
  df <- data.frame(callp, type, spot, yield, time, rate, strike)
  with(df, {
    ifelse(
      tolower(as.character(type)) == "call",
      callp,
      callp - spot * exp(-yield * time) + strike * exp(-rate * time)
    )
  })
}

#' Computes Vega from modified Corrado-Su model.
#'
#' This function computes vega greek for an option based on modified Corrado-Su
#' pricing model.
#'
#' @param type 'call' for call option, any other value for put.
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param sigma volatility of the stock price.
#' @param mu3 skewness.
#' @param mu4 kurtosis.
#'
#' @return
#' Vector with options Vega.
#'
#' @examples
#' csmvega("call", 14, 13, 0.1, 0.15, 0, 0.5, 0.18, 3.53)
#' csmvega("put", 14, 13, 0.1, 0.15, 0, 0.5, 0.18, 2.53)
#'
#' @export
csmvega <- function(type, spot, strike, time, rate, yield, sigma, mu3, mu4) {
  dsig <- 0.0001
  csmu <- csmprice(
    type, spot, strike, time, rate, yield, sigma + dsig, mu3, mu4
  )
  csmd <- csmprice(
    type, spot, strike, time, rate, yield, sigma - dsig, mu3, mu4
  )
  return((csmu - csmd) / (2 * dsig))
}

#' Computes Delta from modified Corrado-Su model.
#'
#' This function computes delta greek for an option with the modified Corrado-Su
#' pricing model.
#'
#' @param type 'call' for call option, any other value for put.
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param sigma volatility of the stock price.
#' @param mu3 skewness.
#' @param mu4 kurtosis.
#'
#' @return
#' Vector with options Delta.
#'
#' @examples
#' csmdelta("call", 14, 13, 0.1, 0.15, 0, 0.5, 0.18, 3.53)
#' csmdelta("put", 14, 13, 0.1, 0.15, 0, 0.5, 0.18, 2.53)
#'
#' @export
csmdelta <- function(type, spot, strike, time, rate, yield, sigma, mu3, mu4) {
  dsig <- 0.001
  csmu <- csmprice(
    type, spot + dsig, strike, time, rate, yield, sigma, mu3, mu4
  )
  csmd <- csmprice(
    type, spot - dsig, strike, time, rate, yield, sigma, mu3, mu4
  )
  return((csmu - csmd) / (2 * dsig))
}

#' Unconstrained Gram-Charlier Distribution Parametrization
#'
#' Creates the unconstrained parametrization for Gram-Charlier Distribution.
#'
#' @param x a numeric value
#' @param mu3p a numeric value
#' @param mu4p a numeric value
#'
#' @name uncons-gram-charlier
#' @importFrom stats approxfun
NULL

#' @rdname uncons-gram-charlier
#' @export
d_poly <- function(x) {
  2 * (x^3 - 3 * x)^2 - 1 * (x^2 - 1) * (x^4 - 6 * x^2 + 3)
}

#' @rdname uncons-gram-charlier
#' @export
mu3 <- function(x) {
  -12 * (x^3 - 3 * x) / d_poly(x)
}

#' @rdname uncons-gram-charlier
#' @export
mu4 <- function(x) {
  24 * (x^2 - 1) / d_poly(x) + 3
}

#' @rdname uncons-gram-charlier
#' @export
create_uncons_regionD <- function() {
  # calculate limits
  a <- -3
  start <- sqrt(3)
  end <- 100
  x <- seq(start^a, end^a, length.out = 100)
  x <- x ^ (1 / a)
  # calculate region D
  dm <- cbind(x = x, mu3 = mu3(x), mu4 = mu4(x))
  rd <- dm[order(dm[, "mu4"]), ]
  rd_curve <- approxfun(rd[, "mu4"], rd[, "mu3"])
  ff <- function(x, a, b) a + (b - a) / (1 + exp(-x))
  function(mu3p, mu4p) {
    mu4 <- ff(mu4p, 3, 7)
    mu3_l <- rd_curve(mu4)
    if (is.na(mu3_l)) {
      c(0, mu4)
    } else {
      mu3_u <- -mu3_l
      mu3 <- ff(mu3p, mu3_l, mu3_u)
      c(mu3, mu4)
    }
  }
}

#' @rdname uncons-gram-charlier
#' @export
uncons_regionD <- create_uncons_regionD()

#' Corrado-Su Model Objective Function in Implied Volatility
#'
#' Objective function to Corrado-Su modelo that minimizes the error in implied
#' volatility.
#'
#' @param par a numeric value
#' @param type option type: `Call` or `Put`
#' @param spot a numeric value
#' @param strike a numeric value
#' @param rate a numeric value
#' @param yield a numeric value
#' @param time a numeric value
#' @param vol a numeric value
#' @param weights a numeric value
#'
#' @export
csm_obj_min_vol <- function(par,
                            type, spot, strike, rate, yield, time,
                            vol, weights = 1) {
  sigma <- par[1]
  params <- uncons_regionD(par[2], par[3])
  mu3 <- params[1]
  mu4 <- params[2]
  pf <- csmprice(type, spot, strike, time, rate, yield, sigma, mu3, mu4)
  yf <- bsmimpvol(pf, type, spot, strike, time, rate, yield)
  sum(((yf - vol) * weights)^2, na.rm = TRUE)
}

#' Corrado-Su Model Objective Function in Price
#'
#' Objective function to Corrado-Su modelo that minimizes the error in price.
#'
#' @param par a numeric value
#' @param type option type: `Call` or `Put`
#' @param spot a numeric value
#' @param strike a numeric value
#' @param rate a numeric value
#' @param yield a numeric value
#' @param time a numeric value
#' @param price a numeric value
#' @param weights a numeric value
#'
#' @export
csm_obj_min_price <- function(par,
                              type, spot, strike, rate, yield, time, price,
                              weights = 1) {
  sigma <- par[1]
  params <- uncons_regionD(par[2], par[3])
  mu3 <- params[1]
  mu4 <- params[2]
  yf <- csmprice(type, spot, strike, time, rate, yield, sigma, mu3, mu4)
  sum(((yf - price) * weights)^2, na.rm = TRUE)
}

#' @rdname csm_obj_min_price
#' @export
csm_obj_grad_price <- function(par,
                               type, spot, strike, rate, yield, time, price,
                               weights = 1) {
  par_ <- uncons_regionD(par[2], par[3])
  par <- c(par[1], par_[1], par_[2])

  #numeric vega for Corrado-Su mod
  dcs.dsig <- csmvega(
    type, spot, strike, time, rate, yield,
    par[1], par[2], par[3]
  )

  #other derivatives for the gradient
  csmw_ <- csmw(par[1], time, par[2], par[3])
  dmod <- csmd(spot, strike, time, rate, yield, par[1], csmw_)
  q3 <- csmq3(spot, par[1], time, dmod, csmw_)
  q4 <- csmq4(spot, par[1], time, dmod, csmw_)

  premium <- csmprice(type, spot, strike, time, rate, yield, par[1], par[2], par[3])
  grad1 <- sum(2 * dcs.dsig * (premium -  price) / (weights ^ 2))
  grad2 <- sum(2 * q3 * (premium -  price) / (weights ^ 2))
  grad3 <- sum(2 * q4 * (premium -  price) / (weights ^ 2))

  ret <- c(grad1, grad2, grad3)
  ret
}

#' Corrado-Su Model Fit in Implied Volatility
#'
#' Function to fit Corrado-Su modelo that minimizes the error in implied
#' volatility.
#'
#' @param par a numeric value
#' @param type option type: `Call` or `Put`
#' @param spot a numeric value
#' @param strike a numeric value
#' @param rate a numeric value
#' @param yield a numeric value
#' @param time a numeric value
#' @param vol a numeric value
#' @param weights a numeric value
#' @param ... additional arguments passed to `optim`
#'
#' @importFrom stats optim
#' @export
csm_fit_min_vol <- function(par,
                            type, spot, strike, rate, yield, time,
                            vol, weights = 1, ...) {
  res <- optim(
    par,
    fn = csm_obj_min_vol, gr = NULL,
    type, spot, strike, rate, yield, time, vol, weights,
    lower = c(1e-1, -Inf, -Inf), upper = Inf,
    method = "L-BFGS-B", ...
  )

  c(res$par[1], uncons_regionD(res$par[2], res$par[3]))
}

#' Corrado-Su Model Fit in Price
#'
#' Function to fit Corrado-Su modelo that minimizes the error in price.
#'
#' @param par a numeric value
#' @param type option type: `Call` or `Put`
#' @param spot a numeric value
#' @param strike a numeric value
#' @param rate a numeric value
#' @param yield a numeric value
#' @param time a numeric value
#' @param price a numeric value
#' @param weights a numeric value
#' @param ... additional arguments passed to `optim`
#'
#' @importFrom stats optim
#' @export
csm_fit_min_price <- function(par,
                              type, spot, strike, rate, yield, time,
                              price, weights = 1, ...) {
  res <- optim(
    par,
    fn = csm_obj_min_price, gr = csm_obj_grad_price,
    method = "L-BFGS-B",
    type, spot, strike, rate, yield, time, price, weights,
    ...
  )

  c(res$par[1], uncons_regionD(res$par[2], res$par[3]))
}

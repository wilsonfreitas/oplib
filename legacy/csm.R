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
  dsig <- 0.001
  csmu <- csmprice(
    type, spot, strike, time, rate, yield, sigma + dsig, mu3, mu4
  )
  csmd <- csmprice(
    type, spot, strike, time, rate, yield, sigma - dsig, mu3, mu4
  )
  return((csmu - csmd) / (2 * dsig))
}


#' Computes Objective Function for Optimizing Modified Corrado-Su Parameters.
#'
#' This is an auxiliary function that provides the objective function used in
#' the optimization to obtain the Modified Corrado-Su parameters.
#'
#' @param par a 3-dimension vector with the parameters
#' @param spot current stock price.
#' @param strike the strike price.
#' @param rate the risk-free interest rate.
#' @param time time to option expiration in years.
#' @param y.data known prices for options.
#' @param sy.data known standard deviation for prices for options.
#'
#' @return
#' Normalized error between prices obtained from data and the ones obtained
#' from model.
#'
#' @export
f.optim.csm <- function(par, spot, strike, rate, time, y.data, sy.data) {
  sigma <- par[1]
  mu3 <- par[2]
  mu4 <- par[3]

  yf <- csmprice(spot, strike, sigma, rate, time, mu3, mu4)
  return(sum(((yf - y.data) / sy.data)^2))
}

#' Computes Inequality Constraint for Optimizing Modified Corrado-Su Parameters.
#'
#' This is an auxiliary function that provides the inequality constraint used in
#' the optimization to obtain the Modified Corrado-Su parameters.
#'
#' @param par a 3-dimension vector with the parameters
#' @param spot current stock price.
#' @param strike the strike price.
#' @param rate the risk-free interest rate.
#' @param time time to option expiration in years.
#' @param y.data known prices for options.
#' @param sy.data known standard deviation for prices for options.
#'
#' @return
#' Constraint to be considered in the optimization.
#'
#' @export
g_ineq <- function(par, spot, strike, rate, time, y.data, sy.data) {
  mu3 <- par[2]
  mu4 <- par[3]
  return(abs(mu3) - skew_gram(mu4))
}

#' Computes Gradient of Inequality Constraint for Optimizing Modified
#' Corrado-Su Parameters.
#'
#' This is an auxiliary function that provides the gradient for inequality
#' constraint
#' used in the optimization to obtain the Modified Corrado-Su parameters.
#'
#' @param par a 3-dimension vector with the parameters
#' @param spot current stock price.
#' @param strike the strike price.
#' @param rate the risk-free interest rate.
#' @param time time to option expiration in years.
#' @param y.data known prices for options.
#' @param sy.data known standard deviation for prices for options.
#'
#' @return
#' Three-dimension vector with the gradient of the inequality constraint.
#'
#' @export
grad_ineq <- function(par, spot, strike, rate, time, y.data, sy.data) {
  return(c(0, 1, 1))
}

#' Computes Gradient of the Objective Function for Optimizing Modified
#' Corrado-Su Parameters.
#'
#' This is an auxiliary function that provides the gradient of the objective
#' function used in the optimization to obtain the Modified Corrado-Su
#' parameters.
#'
#' @param par a 3-dimension vector with the parameters
#' @param type 'call' for call option, any other value for put.
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param y.data known prices for options.
#' @param sy.data known standard deviation for prices for options.
#'
#' @return
#' Matrix with the gradient components of the objective function.
#'
#' @export
grad <- function(par, type, spot, strike, time, rate, yield, y.data, sy.data) {
  # numeric vega for Corrado-Su mod
  dcs.dsig <- csmvega(
    type, spot, strike, time, rate, yield,
    par[1], par[2], par[3]
  )

  # other derivatives for the gradient
  csmw <- csmw(par[1], time, par[2], par[3])
  dmod <- csmd(spot, strike, 0, par[1], time, csmw)
  q3 <- csmq3(spot, par[1], time, dmod, csmw)
  q4 <- csmq4(spot, par[1], time, dmod, csmw)

  premium <- csmprice(spot, strike, par[1], rate, time, par[2], par[3])
  grad1 <- sum(2 * dcs.dsig * (premium - y.data) / (sy.data^2))
  grad2 <- sum(2 * q3 * (premium - y.data) / (sy.data^2))
  grad3 <- sum(2 * q4 * (premium - y.data) / (sy.data^2))

  return(cbind(grad1, grad2, grad3))
}


#' Computes Objective Function for Optimizing Bi-Modal Corrado-Su Parameters.
#'
#' This is an auxiliary function that provides the objective function used in
#' the optimization to obtain the Bi-Modal Corrado-Su parameters.
#'
#' @param par a 7-dimension vector with the parameters
#' @param spot current stock price.
#' @param strike the strike price.
#' @param rate the risk-free interest rate.
#' @param time time to option expiration in years.
#' @param y.data known prices for options.
#' @param sy.data known standard deviation for prices for options.
#'
#' @return
#' Normalized error between prices obtained from data and the ones obtained
#' from model.
#'
#' @export
f.optim.2csm <- function(par, spot, strike, rate, time, y.data, sy.data) {
  sigma1 <- par[1]
  mu31 <- par[2]
  mu41 <- par[3]
  sigma2 <- par[4]
  mu32 <- par[5]
  mu42 <- par[6]
  w <- par[7]

  yf1 <- csmprice(spot, strike, sigma1, rate, time, mu31, mu41)
  yf2 <- csmprice(spot, strike, sigma2, rate, time, mu32, mu42)
  yf <- yf1 * w + yf2 * (1 - w)
  return(sum(((yf - y.data) / sy.data)^2))
}

#' @rdname f.optim.2csm
#' @export
csm.2opt <- function(par, spot, strike, rate, time) {
  sigma1 <- par[1]
  mu31 <- par[2]
  mu41 <- par[3]
  sigma2 <- par[4]
  mu32 <- par[5]
  mu42 <- par[6]
  w <- par[7]

  yf1 <- csmprice(spot, strike, sigma1, rate, time, mu31, mu41)
  yf2 <- csmprice(spot, strike, sigma2, rate, time, mu32, mu42)
  return(yf1 * w + yf2 * (1 - w))
}

#' Computes Inequality Constraint for Optimizing Bi-Modal Corrado-Su Parameters.
#'
#' This is an auxiliary function that provides the inequality constraint used in
#' the optimization to obtain the Bi-Modal Corrado-Su parameters.
#'
#' @param par a 7-dimension vector with the parameters
#' @param spot current stock price.
#' @param strike the strike price.
#' @param rate the risk-free interest rate.
#' @param time time to option expiration in years.
#' @param y.data known prices for options.
#' @param sy.data known standard deviation for prices for options.
#'
#' @return
#' Constraint to be considered in the optimization.
#'
#' @export
g_ineq2 <- function(par, spot, strike, rate, time, y.data, sy.data) {
  mu31 <- par[2]
  mu41 <- par[3]
  mu32 <- par[5]
  mu42 <- par[6]

  tst1 <- abs(mu31) - skew_gram(mu41)
  tst2 <- abs(mu32) - skew_gram(mu42)
  if (tst1 > 0 || tst2 > 0) {
    # print('tst > 0')
    return(1)
  } else {
    # print('tst < 0')
    return(-1)
  }
}

#' Computes Gradient of Inequality Constraint for Optimizing Bi-Modal
#' Corrado-Su Parameters.
#'
#' This is an auxiliary function that provides the gradient for inequality
#' constraint
#' used in the optimization to obtain the Bi-Modal Corrado-Su parameters.
#'
#' @param par a 7-dimension vector with the parameters
#' @param spot current stock price.
#' @param strike the strike price.
#' @param rate the risk-free interest rate.
#' @param time time to option expiration in years.
#' @param y.data known prices for options.
#' @param sy.data known standard deviation for prices for options.
#'
#' @return
#' Three-dimension vector with the gradient of the inequality constraint.
#'
#' @export
grad_ineq2 <- function(par, spot, strike, rate, time, y.data, sy.data) {
  return(c(0, 1, 1, 0, 1, 1, 0))
}

#' Computes Gradient of the Objective Function for Optimizing Bi-Modal
#' Corrado-Su Parameters.
#'
#' This is an auxiliary function that provides the gradient of the objective
#' function used in the optimization to obtain the Bi-Modal Corrado-Su
#' parameters.
#'
#' @param par a 7-dimension vector with the parameters
#' @param type 'call' for call option, any other value for put.
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param y.data known prices for options.
#' @param sy.data known standard deviation for prices for options.
#'
#' @return
#' Matrix with the gradient components of the objective function.
#'
#' @export
grad2 <- function(par, type, spot, strike, time, rate, yield, y.data, sy.data) {

  # numeric vega for Corrado-Su mod
  dcs.dsig1 <- csmvega(
    type, spot, strike, time, rate,
    yield, par[1], par[2], par[3]
  )
  dcs.dsig2 <- csmvega(
    type, spot, strike, time, rate,
    yield, par[4], par[5], par[6]
  )

  # other derivatives for the gradient
  csmw1 <- csmw(par[1], time, par[2], par[3])
  dmod1 <- csmd(spot, strike, 0, par[1], time, csmw1)
  q31 <- csmq3(spot, par[1], time, dmod1, csmw1)
  q41 <- csmq4(spot, par[1], time, dmod1, csmw1)

  csmw2 <- csmw(par[4], time, par[5], par[6])
  dmod2 <- csmd(spot, strike, 0, par[4], time, csmw2)
  q32 <- csmq3(spot, par[4], time, dmod2, csmw2)
  q42 <- csmq4(spot, par[4], time, dmod2, csmw2)

  premium1 <- csmprice(spot, strike, par[1], rate, time, par[2], par[3])
  premium2 <- csmprice(spot, strike, par[4], rate, time, par[5], par[6])

  grad1 <- sum(2 * dcs.dsig1 * par[7] * (premium1 - y.data) / (sy.data^2))
  grad2 <- sum(2 * q31 * par[7] * (premium1 - y.data) / (sy.data^2))
  grad3 <- sum(2 * q41 * par[7] * (premium1 - y.data) / (sy.data^2))

  grad4 <- sum(2 * dcs.dsig2 * par[7] * (premium2 - y.data) / (sy.data^2))
  grad5 <- sum(2 * q32 * par[7] * (premium2 - y.data) / (sy.data^2))
  grad6 <- sum(2 * q42 * par[7] * (premium2 - y.data) / (sy.data^2))

  grad7 <- sum((premium1 - premium2) * (premium2 - y.data) / (sy.data^2))
  # print('grd')
  # print(cbind(grad1, grad2, grad3, grad4, grad5, grad6, grad7))

  return(cbind(grad1, grad2, grad3, grad4, grad5, grad6, grad7))
}

#' Modified Corrado-Su std. Deltas
#'
#' Finds and returns a dataframe with the strikes and implied
#' volatilites for the provided set of deltas or for the
#' default deltas.
#'
#' @param type Option type, "call" for calls or anything else for puts
#' @param spot underlying spot price at a given date
#' @param time time in years from a given date to a maturity date
#' @param rate interest free rate
#' @param yield dividend yield
#' @param sigma volatility of the stock price
#' @param mu3 skewness
#' @param mu4 kurtosis
#' @param deltas deltas for which strikes and vols should be calculated, 1-99.
#' Defaults to usual delta values.
#'
#' @return
#' Dataframe with columns delta, strike and impvol
#' @export
csm_deltas <- function(type, spot, time, rate, yield, sigma, mu3, mu4,
                       deltas = DELTAS) {
  ft <- spot * exp((rate - yield) * time)
  eps <- 1e-8

  strk_cs <- sapply(deltas / 100, function(delta) {
    uniroot(
      f = csm_strikes, interval = c(1, spot * 10), tol = eps,
      delta = delta, type = type, spot = spot, time = time,
      rate = rate, yield = yield, sigma = sigma, mu3 = mu3, mu4 = mu4, ft = ft
    )$root
  })

  impvol_cs <- blackimpvol(
    csmprice(type, spot, strk_cs, time, rate, yield, sigma, mu3, mu4),
    type, ft, strk_cs, time, rate
  )


  data.frame(delta = deltas, strike = strk_cs, impvol = impvol_cs)
}

csm_strikes <- function(x, delta, type, spot, strike, time, rate, yield, sigma,
                        mu3, mu4, ft) {
  csp <- csmprice(type, spot, x, time, rate, yield, sigma, mu3, mu4)
  if (csp < 0) {
    return(999)
  }
  varx <- blackimpvol(csp, type, ft, x, time, rate)
  x - blackstrike(ft, delta, varx, time)
}
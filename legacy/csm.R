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
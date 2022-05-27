#' z auxiliary function for SABR(Stochastic Alpha Beta Rho) volatility model.
#'
#' SABR model formula is broken in auxiliary functions, and z is one
#' of them.
#'
#' @param fut current future price.
#' @param strike the strike price.
#' @param alpha parameter from model.
#' @param beta parameter from model.
#' @param nu parameter from model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with the z formula values.
#'
#' @seealso
#' \code{\link{x.func}} for x auxiliary function from SABR model.
#' \code{\link{sabrimpvol}} for volatility fit with SABR model.
#'
#' @examples
#' z.func(10, 9, 0.5, 0.5, 0.5)
#' z.func(50:150 / 10, 9, 0.5, 0.5, 0.5)
#' z.func(50, 45, 0.4325, 0.65846, 1.5)
#'
#' @export
z.func <- function(fut, strike, alpha, beta, nu) {
  return((nu / alpha) * ((fut * strike)^(0.5 * (1 - beta))) * log(fut / strike))
}

#' x auxiliary function for SABR(Stochastic Alpha Beta Rho) volatility model.
#'
#' SABR model formula is separated in auxiliary functions, and x is one
#' of them.
#'
#' @param z other auxiliary function for SABR.
#' @param rho parameter from model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with the x formula values.
#'
#' @seealso
#' \code{\link{z.func}} for z auxiliary function from SABR model.
#' \code{\link{sabrimpvol}} for volatility fit with SABR model.
#'
#' @examples
#' x.func(z.func(10, 9, 0.5, 0.5, 0.5), 0.5)
#' x.func(z.func(50, 45, 0.4325, 0.65846, 1.5), -0.15678)
#' x.func(-100:100 / 100, 0)
#'
#' @export
x.func <- function(z, rho) {
  fct <- 1 - 2 * rho * z + z^2
  fct <- ifelse(fct < 0, 0, fct)
  fct <- sqrt(fct)
  # den <- 1 - rho
  log.arg <- (fct + z - rho) / (1 - rho)
  log.arg <- ifelse(log.arg < 0, 0, log.arg)
  # TODO: verificar caso em que 1 - rho == 0
  return(log(log.arg))
}

#' Volatility fit with SABR(Stochastic Alpha Beta Rho) model.
#'
#' SABR model for volatility market data fit implementation.
#'
#' @param spot current stock price.
#' @param rate the risk-free interest rate.
#' @param strike the strike price.
#' @param bizdays time to option expiration in years, based in working days.
#' @param alpha parameter from model.
#' @param beta parameter from model.
#' @param rho parameter from model.
#' @param nu parameter from model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with fit volatility, according to SABR model.
#'
#' @seealso
#' \code{\link{z.func}} for z auxiliary function from SABR model.
#' \code{\link{x.func}} for x auxiliary function from SABR model.
#'
#' @examples
#' sabrimpvol(5, 0.16, 4.5, 1, 0.5, 0.5, 0.5, 0.5)
#' sabrimpvol(5, 0.16, 4.5, 1, 0.5, 0:100 / 100, 0.5, 0.5)
#' sabrimpvol(5, .16, 4.5, 1, .5, .5, -100:100 / 100, .5)
#'
#' @export
sabrimpvol <- function(spot = NULL, rate = NULL, fut = NULL, strike, bizdays, alpha, beta, rho, nu) {
  if (is.null(fut)) {
    fut <- spot * exp(rate * bizdays)
  }

  z <- z.func(fut, strike, alpha, beta, nu)
  x <- x.func(z, rho)

  dn1 <- (fut * strike)^(0.5 * (1 - beta))
  dn2 <- 1 + (1 / 24) * ((1 - beta)^2) * ((log(fut / strike))^2) +
    (1 / 1920) * ((1 - beta)^4) * (log(fut / strike))^4

  term1 <- (alpha / (dn1 * dn2)) * (z / x)

  fct <- (1 / 24) * ((1 - beta)^2) * (alpha^2) / ((fut * strike)^(1 - beta)) +
    (1 / 4) * (rho * beta * nu * alpha) / ((fut * strike)^(0.5 * (1 - beta))) +
    (1 / 24) * (2 - 3 * (rho^2)) * nu^2

  term2 <- 1 + (fct * bizdays)

  impvol <- term1 * term2
  if (any(is.nan(impvol))) impvol[is.nan(impvol)] <- 0
  return(impvol)
}

#' Numerical derivative of implied volatility formula from SABR method with
#' respect to alpha.
#'
#' Symmetric Difference Quotient method of numerical differentiation applied
#' in the volatility formula from SABR method, with respect to alpha.
#'
#' @param spot current stock price.
#' @param rate the risk-free interest rate.
#' @param strike the strike price.
#' @param bizdays time to option expiration in years, based in working days.
#' @param alpha parameter from model.
#' @param beta parameter from model.
#' @param rho parameter from model.
#' @param nu parameter from model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with the derivatives.
#'
#' @seealso
#' \code{\link{sabrimpvol}} for SABR volatility formula.
#'
#' @examples
#' dsabrimpvol.dalpha(5, 0.16, 4.5, 1, 0.5, 0.5, 0.5, 0.5)
#' dsabrimpvol.dalpha(5, 0.16, 4.5, 1, 0.5, 0:100 / 100, 0.5, 0.5)
#' dsabrimpvol.dalpha(5, .16, 4.5, 1, .5, .5, -100:100 / 100, .5)
#'
#' @export
dsabrimpvol.dalpha <- function(spot, rate, strike, bizdays, alpha, beta, rho, nu) {
  d.alpha <- 1e-6
  impvol.up <- sabrimpvol(spot, rate, strike, bizdays, alpha + d.alpha, beta, rho, nu)
  impvol.down <- sabrimpvol(spot, rate, strike, bizdays, alpha - d.alpha, beta, rho, nu)
  return((impvol.up - impvol.down) / (2 * d.alpha))
}


#' SABR model fit
#'
#' Finds and returns SABR model parameters alpha, beta, nu and rho for
#' a given set of strikes and sigmas.
#'
#' @param spot current stock price. If not provided will try to use the future parameter.
#' @param rate the risk-free interest rate. If not provided will try to use the future parameter.
#' @param future the future price. If not provided will calculate future from spot and rate.
#' @param tau time to option expiration in years, based in working days.
#' @param strikes options strikes
#' @param sigmas options sigmas
#' @param initial_guess named vector in order alpha, nu, rho, beta of initial parameters for optimization (optional)
#'
#' @return
#' List with fitted parameters alpha, beta, nu, rho
#'
#' @seealso
#' \code{\link{sabrimpvol}} for SABR volatility
#'
#' @export
sabr_fit <- function(future = NULL, strikes, tau, sigmas, spot = NULL, rate = NULL, initial_guess = NULL) {
  if (is.null(future) && (is.null(spot) || is.null(rate))) stop("Either future price or spot and rate must be provided.")
  if (is.null(future)) future <- spot * exp(rate * tau)
  if (is.null(initial_guess)) initial_guess <- c(alpha = ALPHA, nu = NU, rho = RHO, beta = BETA)
  fit <- nloptr::lbfgs(
    x0 = initial_guess,
    fn = sabr_optim(future, strikes, tau, sigmas),
    lower = c(0.0001, 0.0001, -1, 0),
    upper = c(20, Inf, 0.99, 1)
  )
  list(alpha = fit$par[1], nu = fit$par[2], rho = fit$par[3], beta = fit$par[4])
}

sabr_optim <- function(future, strikes, tau, sigmas) {
  function(par) {
    alpha <- par[1]
    nu <- par[2]
    rho <- par[3]
    beta <- par[4]
    sabr_sigs <- sapply(strikes, function(strike) {
      sabrimpvol(
        fut = future, strike = strike, bizdays = tau,
        alpha = alpha, beta = beta, rho = rho, nu = nu
      )
    })
    sum((sigmas - sabr_sigs)^2)
  }
}


#' SABR model std. Deltas
#'
#' Finds and returns a dataframe with the strikes and implied
#' volatilites for the provided set of deltas or for the
#' default deltas.
#'
#' @param alpha a SABR model fit alpha parameter
#' @param beta a SABR model fit beta parameter
#' @param rho a SABR model fit rho parameter
#' @param nu a SABR model fit nu parameter
#' @param S underlying spot price at a given date
#' @param tau time in years from a given date to a maturity date
#' @param r interest free rate
#' @param q yield rate
#' @param future future price. Can be given instead of spot/rate/yield.
#' @param deltas deltas for which strikes and vols should be calculated. Defaults to usual
#' delta values.
#'
#' @return
#' Dataframe with columns delta, strike and impvol
#'
#' @export
sabr_deltas <- function(alpha, beta, rho, nu, S = NULL, tau, r = NULL,
                        future = NULL, q = NULL, deltas = DELTAS) {
  eps <- 1e-8
  if (is.null(future)) {
    ft <- S * exp((r - q) * tau)
  } else {
    ft <- future
  }

  if (is.null(ft)) stop("Future price or spot, rate and yield must be provided.")

  strk_sabr <- sapply(deltas / 100, function(delta) {
    uniroot(
      f = sabr_strikes, interval = c(1, ft * 10), tol = eps,
      t = tau, delta = delta, ft = ft, alpha = alpha,
      beta = beta, rho = rho, nu = nu
    )$root
  })

  ## multiroot não funciona bem, varx <- sabrimpvol(...) não usa vetor de deltas
  # strk_sabr <- multiroot(func_ = sabr_strikes, interval = c(1, S*10), tolerance = eps,
  # t = tau, delta = deltas, ft = ft, r = r, spot = S, alpha = alpha,
  # beta = beta, rho = rho, nu = nu)$root


  impvol_sabr <- sabrimpvol(
    fut = ft, strike = strk_sabr, bizdays = tau,
    alpha = alpha, beta = beta, rho = rho, nu = nu
  )

  data.frame(delta = deltas, strike = strk_sabr, impvol = impvol_sabr)
}


sabr_strikes <- function(x, delta, ft, t, alpha, beta, rho, nu) {
  varx <- sabrimpvol(
    fut = ft, strike = x, bizdays = t,
    alpha = alpha, beta = beta, rho = rho, nu = nu
  )
  x - blackstrike(ft, delta, varx, t)
}
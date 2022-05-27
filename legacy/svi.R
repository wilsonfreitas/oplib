#' SVI Variance
#'
#' Calculates the variance as given by the SVI model for
#' the provided parameters.
#'
#' @param a 'a' in the SVI model (~~ instantaneous variance)
#' @param b 'b' in the SVI model (~~ speed of variance mean reversion)
#' @param m 'm' in the SVI model (~~ long term expected variance)
#' @param rho 'rho' in the SVI model (~~ correlation of variance and price)
#' @param x 'x' in the SVI model, given by ln(K/F) where K is an options strike price and F its underlying's future
#' @param sigma 'sigma' in the SVI model (~~ volatility of variance)
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Volatility value as given by the SVI model
#'
#' @seealso
#' \code{\link{SVI_fit}} for SVI model fitting.
#' \code{\link{SVI_deltas}} for calculating strikes and imp. vols. for given deltas
#'
#' @export
SVI_var <- function(a, b, m, rho, x, sigma) {
  a + b * (rho * (x - m) + sqrt((x - m)^2 + sigma^2))
}

#' SVI fit
#'
#' Fits a SVI model for the given smile as defined by strikes, futures
#' and their respective variances. Time to maturity should also be provided
#' for the non-arbitrage restriction. Returns the fitted parameters.
#'
#' @param strikes options strikes
#' @param future options underlying future If NULL will try to use spot and rate.
#' @param variances variances for the given strikes
#' @param ttm time to maturity in years
#' @param spot underlying spot priced, used in conjunction with rate if future price is not provided
#' @param rate risk free rate, used for calculating a future price in conjunction with spot if not provided
#'
#' @return
#' A list with the SVI model parameters a, b, x, m, rho and sigma.
#'
#' @seealso
#' \code{\link{SVI_var}} for SVI model variance calculation
#' \code{\link{SVI_deltas}} for calculating strikes and imp. vols. for given deltas
#'
#' @export
SVI_fit <- function(strikes, future = NULL, variances, ttm, spot = NULL, rate = NULL, initial_guess = NULL) {
  if (is.null(future) && (is.null(spot) || is.null(rate))) stop("Either future price or spot and rate must be provided.")
  if (is.null(future)) future <- spot * exp(rate * ttm)

  x <- log(strikes / future)
  tol <- length(variances) * min(variances)^2

  # SVI restrictions
  # 0 <= a <= max(variances)
  # 0 <= b <= 4/(ttm * (1 + abs(rho)))
  # -1 <= rho <= 1
  # min(x) <= m <= max(x)
  # 0 < sigma < 10
  # a, b, m, rho, sigma
  eval_g0 <- function(par) {
    c(par[2] - 4 / (ttm * (1 + abs(par[4]))))
  }


  lo <- c(0, 0, min(x), -1, 1e-12)
  hi <- c(max(c(variances, SVI.A)), 4 / (ttm * (1 + 1e-4)), max(x), 1, 10)

  if (is.null(initial_guess)) {
    initial_guess <- c(
      a = SVI.A, b = SVI.B,
      m = min(x), rho = SVI.RHO, sigma = SVI.SIGMA
    )

    fit <- nloptr::nloptr(
      x0 = initial_guess,
      eval_f = SVI_minsq(variances, x),
      lb = lo,
      ub = hi,
      eval_g_ineq = eval_g0,
      opts = list(
        algorithm = "NLOPT_GN_ISRES",
        print_level = 0,
        xtol_rel = 1e-4,
        maxeval = 2000
      )
    )

    initial_guess <- fit$solution
  }

  fit <- nloptr::nloptr(
    x0 = initial_guess,
    eval_f = SVI_minsq(variances, x),
    lb = lo,
    ub = hi,
    eval_g_ineq = eval_g0,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      print_level = 0,
      xtol_rel = 1e-4,
      maxeval = 2000
    )
  )

  curr_param <- fit$solution

  list(
    a = curr_param[1], b = curr_param[2], m = curr_param[3],
    rho = curr_param[4], sigma = curr_param[5]
  )
}


SVI_minsq <- function(variances, x) {
  function(...) {
    par <- c(...)
    a <- par[1]
    b <- par[2]
    m <- par[3]
    rho <- par[4]
    sigma <- par[5]
    SVI <- SVI_var(a, b, m, rho, x, sigma)
    sum((variances - SVI)^2)
  }
}



#' SVI fit std. Deltas
#'
#' Finds and returns a dataframe with the strikes and implied
#' volatilites for the provided set of deltas or for the
#' default deltas.
#'
#' @param a a SVI model fit a parameter
#' @param b a SVI model fit b parameter
#' @param m a SVI model fit m parameter
#' @param rho a SVI model fit rho parameter
#' @param sigma a SVI model fit sigma parameter
#' @param future future price
#' @param t time to maturity
#' @param deltas deltas for which strikes and vols should be calculated. Defaults to usual
#' delta values.
#'
#' @return
#' Dataframe with columns delta, strike and impvol
#'
#' @export
SVI_deltas <- function(a, b, m, rho, sigma, future, t, deltas = DELTAS) {
  eps <- 1e-8

  strk_svi <- sapply(deltas / 100, function(delta) {
    uniroot(
      f = svi_strikes, interval = c(1, future * 10), tol = eps,
      t = t, delta = delta, a = a, b = b,
      m = m, rho = rho, sigma = sigma, ft = future
    )$root
  })

  impvol_svi <- sqrt(SVI_var(a, b, m, rho, log(strk_svi / future), sigma))

  data.frame(delta = deltas, strike = strk_svi, impvol = impvol_svi)
}


svi_strikes <- function(x, delta, a, b, m, rho, sigma, ft, t) {
  varx <- sqrt(SVI_var(a, b, m, rho, log(x / ft), sigma))
  x - blackstrike(ft, delta, varx, t)
}
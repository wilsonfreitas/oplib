#' SVI Variance
#'
#' Calculates the variance as given by the SVI model for
#' the provided parameters.
#'
#' @param a 'a' in the SVI model (~~ instantaneous variance)
#' @param b 'b' in the SVI model (~~ speed of variance mean reversion)
#' @param m 'm' in the SVI model (~~ long term expected variance)
#' @param rho 'rho' in the SVI model (~~ correlation of variance and price)
#' @param x 'x' in the SVI model, given by ln(K/F) where K is an options strike
#'        price and F its underlying's future
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
#' \code{\link{svi_fit}} for SVI model fitting.
#'
#' @export
svi_var <- function(a, b, m, rho, x, sigma) {
  a + b * (rho * (x - m) + sqrt((x - m)^2 + sigma^2))
}

# Default SVI initial parameters
SVI.A <- 1e-12
SVI.B <- 1e-8
SVI.SIGMA <- 1e-8
SVI.RHO <- -0.8
SVI.M <- 1e-8
SVI.ERR <- 0.02

#' SVI fit
#'
#' Fits a SVI model for the given smile as defined by strikes
#' and their respective variances. Time to maturity should also be provided
#' for the non-arbitrage restriction.
#'
#' @param variance variances for the given strikes
#' @param strike options strikes
#' @param time time to maturity in years
#' @param spot underlying spot price
#' @param rate risk free rate
#' @param yield yield rate
#' @param initial_guess Initial guess used in the optimization process
#'
#' @return
#' A list with the SVI model parameters a, b, x, m, rho and sigma.
#'
#' @seealso
#' \code{\link{svi_var}} for SVI model variance calculation
#'
#' @export
svi_fit <- function(variance,
                    spot,
                    strike,
                    time,
                    rate,
                    yield = 0,
                    initial_guess = NULL) {
  future <- spot * exp(rate * time)

  x <- log(strike / future)
  tol <- length(variance) * min(variance)^2

  # SVI restrictions
  # 0 <= a <= max(variance)
  # 0 <= b <= 4/(time * (1 + abs(rho)))
  # -1 <= rho <= 1
  # min(x) <= m <= max(x)
  # 0 < sigma < 10
  # a, b, m, rho, sigma
  eval_g0 <- function(par) {
    c(par[2] - 4 / (time * (1 + abs(par[4]))))
  }

  lo <- c(0, 0, min(x), -1, 1e-12)
  hi <- c(max(c(variance, SVI.A)), 4 / (time * (1 + 1e-4)), max(x), 1, 10)

  if (is.null(initial_guess)) {
    initial_guess <- c(
      a = SVI.A, b = SVI.B,
      m = min(x), rho = SVI.RHO, sigma = SVI.SIGMA
    )

    fit <- nloptr::nloptr(
      x0 = initial_guess,
      eval_f = SVI_minsq(variance, x),
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
    eval_f = SVI_minsq(variance, x),
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

  list(
    a = fit$solution[1],
    b = fit$solution[2],
    m = fit$solution[3],
    rho = fit$solution[4],
    sigma = fit$solution[5]
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
    SVI <- svi_var(a, b, m, rho, x, sigma)
    sum((variances - SVI)^2)
  }
}
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
#' @param deltas deltas for which strikes and vols should be calculated.
#'        Defaults to usual delta values.
#'
#' @return
#' Dataframe with columns delta, strike and impvol
#'
#' @export
svi_deltas <- function(a, b, m, rho, sigma, future, t, deltas = DELTAS) {
  eps <- 1e-8

  strk_svi <- sapply(deltas / 100, function(delta) {
    uniroot(
      f = svi_strikes, interval = c(1, future * 10), tol = eps,
      t = t, delta = delta, a = a, b = b,
      m = m, rho = rho, sigma = sigma, ft = future
    )$root
  })

  impvol_svi <- sqrt(svi_var(a, b, m, rho, log(strk_svi / future), sigma))

  data.frame(delta = deltas, strike = strk_svi, impvol = impvol_svi)
}

svi_strikes <- function(x, delta, a, b, m, rho, sigma, ft, t) {
  varx <- sqrt(svi_var(a, b, m, rho, log(x / ft), sigma))
  x - blackstrike(ft, delta, varx, t)
}
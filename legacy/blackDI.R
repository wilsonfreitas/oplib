#' Computes DI future option price via Black-Scholes formula.
#'
#' Black-Scholes formula implementation, pricing DI future options.
#'
#' @param type 'call' for call option, 'put' for put option. May be abbreviated.
#' @param short_price unit price, based on the date of option expiration.
#' @param long_price unit price, based on the date of underlying asset expiration.
#' @param strike the strike price.
#' @param short_t_cur time to option expiration in years, based in current days.
#' @param long_t_cur time to underlying asset expiration in years, based in
#' current days.
#' @param short_t_biz time to option expiration in years, based in working days.
#' @param long_t_biz time to underlying asset expiration in years, based in
#' working days.
#' @param sigma volatility of the underlying asset price.
#'
#' @return
#' DI future option price.
#'
#' @seealso
#' \code{\link{bsmprice}} for pricing options on stocks, with
#' Black-Scholes-Merton model. \\
#' \code{\link{blackprice}} for pricing options
#' on futures, with Black-76 model.
#'
#' @examples
#' blackDI("call", 110000, 100000, 0.14, 1, 2, 1, 2, 0.15)
#' blackDI("p", 100000, 90000, 0.08, 1, 2, 1, 2, 0.15)
#'
#' @export

blackDI <- function(type = c("call", "put"), short_price, long_price, strike,
                    short_t_cur, long_t_cur, short_t_biz, long_t_biz, sigma) {
  type <- match.arg(type)
  df <- data.frame(
    short_price, long_price, strike,
    short_t_cur, long_t_cur, short_t_biz, long_t_biz, sigma
  )
  df <- within(df, {
    Kast <- ((1 + strike)^(long_t_biz - short_t_biz) - 1) * 1 / (long_t_cur - short_t_cur)
    Szero <- (short_price / long_price - 1) * 1 / (long_t_cur - short_t_cur)
    delta <- (long_price * (long_t_cur - short_t_cur)) / (1 + Kast * (long_t_cur - short_t_cur))
    d1 <- (log(Szero / Kast) + (sigma * sigma) / 2 * short_t_cur) / (sigma * (short_t_cur)^0.5)
    d2 <- d1 - sigma * (short_t_cur)^0.5
  })
  with(
    df,
    ifelse(type == "call",
      delta * (Szero * pnorm(d1) - Kast * pnorm(d2)),
      delta * (-Szero * pnorm(-d1) + Kast * pnorm(-d2))
    )
  )
}
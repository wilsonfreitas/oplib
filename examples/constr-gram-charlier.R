library(rb3)
library(bizdays)
library(tidyverse)
devtools::load_all()

d_poly <- function(x) {
  2 * (x^3 - 3 * x)^2 - 1 * (x^2 - 1) * (x^4 - 6 * x^2 + 3)
}

mu3 <- function(x) {
  -12 * (x^3 - 3 * x) / d_poly(x)
}

mu4 <- function(x) {
  24 * (x^2 - 1) / d_poly(x) + 3
}

regionD <- function(x) {
  dm <- cbind(x = x, mu3 = mu3(x), mu4 = mu4(x))
  idx <- order(dm[, "mu4"])
  dm[idx, ]
}

adj_seq <- function(start, end, a = -3, length.out = 100) {
  x <- seq(start^a, end^a, length.out = length.out)
  x ^ (1 / a)
}

create_uncons_regionD <- function() {
  x <- adj_seq(sqrt(3), 100)
  rd <- regionD(x)
  rd_curve <- approxfun(rd[, "mu4"], rd[, "mu3"])
  ff <- function(x, a, b) a + (b - a) / (1 + exp(-x))
  function(mu3p, mu4p) {
    mu4 <- ff(mu4p, 3, 7)
    mu3_l <- rd_curve(mu4)
    mu3_u <- -mu3_l
    mu3 <- ff(mu3p, mu3_l, mu3_u)
    c(mu3, mu4)
  }
}

uncons_regionD <- create_uncons_regionD()

# ----

refdate <- preceding(Sys.Date() - 1, "Brazil/B3")
ch <- cotahist_get(refdate, "daily")
yc <- yc_get(refdate)
op <- cotahist_equity_options_superset(ch, yc)

symbol_ <- "B3SA3"
op1 <- op |>
  filter(
    symbol.underlying == symbol_
  )

maturities <- unique(op1$maturity_date) |> sort()
close_underlying <- op1$close.underlying[1]

op_vol <- op1 |>
  filter(maturity_date %in% maturities[1]) |>
  mutate(
    biz_days = bizdays(
      refdate, following(maturity_date, "Brazil/B3"), "Brazil/B3"
    ),
    time_to_maturity = biz_days / 252,
    rate = log(1 + r_252),
    bsm_impvol = bsmimpvol(
      close, type, close.underlying, strike, time_to_maturity, rate, 0
    ),
    delta = bsmdelta(
      type, close.underlying, strike, time_to_maturity, rate, 0, bsm_impvol
    ),
    adj_delta = ifelse(str_to_lower(type) == "call", delta, 1 + delta)
  ) |>
  select(
    symbol, volume,
    type, close.underlying, strike, time_to_maturity, rate,
    biz_days, close, high, low, bsm_impvol, delta, adj_delta
  )

# ----

f.optim.csm <- function(par, type, spot, strike, rate, time, y.data, sy.data) {
  sigma <- par[1]
  params <- uncons_regionD(par[2], par[3])
  params <- c(par[2], par[3])
  yf <- csmprice(type, spot, strike, time, rate, 0, sigma, params[1], params[2])
  return(sum(((yf - y.data) / sy.data)^2))
}

# ----

comb <- expand.grid(
  mu3 = seq(-1.2, 1.2, length.out = 10),
  mu4 = seq(3, 7, length.out = 10)
)

res <- with(op_vol, {
  optim(
    par = c(0.5, 0, 3), fn = f.optim.csm, gr = NULL,
    type, close.underlying, strike, rate, time_to_maturity, close, 1
  )
})

op_vol_f <- op_vol |>
  mutate(
    theo_price = csmprice(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      res$par[1], res$par[2], res$par[3]
    ),
    csm_impvol = bsmimpvol(
      theo_price, type, close.underlying, strike, time_to_maturity, rate, 0
    )
  )

View(op_vol_f)

op_vol_f |>
  ggplot(aes(x = strike, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = csm_impvol), colour = "black") +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")
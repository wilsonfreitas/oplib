
library(rb3)
library(bizdays)
library(tidyverse)

f.optim.csm <- function(par, type, spot, strike, rate, time, y.data, sy.data) {
  sigma <- par[1]
  mu3 <- par[2]
  mu4 <- par[3]
  yf <- csmprice(type, spot, strike, time, rate, 0, sigma, mu3, mu4)
  return(sum(((yf - y.data) / sy.data)^2))
}

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

# no weights ----

res <- with(op_vol, {
  optim(
    par = c(0.3, 1, 4), fn = f.optim.csm, gr = NULL,
    type, close.underlying, strike, rate, time_to_maturity, close, 1
  )
})

op_csm <- op_vol |>
  mutate(
    theo_price = csmprice(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      res$par[1], res$par[2], res$par[3]
    ),
    csm_impvol = bsmimpvol(
      theo_price, type, close.underlying, strike, time_to_maturity, rate, 0
    )
  )

op_csm |>
  ggplot(aes(x = strike, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = csm_impvol), colour = "black") +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")

op_csm |>
  mutate(error = close - theo_price) |>
  ggplot(aes(x = strike, y = error)) +
  geom_point()

op_csm |>
  mutate(error = bsm_impvol - csm_impvol) |>
  ggplot(aes(x = strike, y = error)) +
  geom_point()

# volume weights ----

res <- with(op_vol, {
  optim(
    par = c(0.3, 1, 4), fn = f.optim.csm, gr = NULL,
    type, close.underlying, strike, rate, time_to_maturity, close, 1 / volume,
    lower = c(0, -Inf, 3), upper = Inf, method = "L-BFGS-B"
  )
})

op_csm <- op_vol |>
  mutate(
    theo_price = csmprice(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      res$par[1], res$par[2], res$par[3]
    ),
    csm_impvol = bsmimpvol(
      theo_price, type, close.underlying, strike, time_to_maturity, rate, 0
    )
  )

op_csm |>
  ggplot(aes(x = strike, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = csm_impvol), colour = "black") +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")

op_csm |>
  mutate(error = close - theo_price) |>
  ggplot(aes(x = strike, y = error)) +
  geom_point()

op_csm |>
  mutate(error = bsm_impvol - csm_impvol) |>
  ggplot(aes(x = strike, y = error)) +
  geom_point()

# min vol error ----

f_obj_csm_min_vol <- function(par,
                              type, spot, strike, rate, time,
                              y.data, sy.data) {
  print(par)
  sigma <- par[1]
  mu3 <- par[2]
  mu4 <- par[3]
  pf <- csmprice(type, spot, strike, time, rate, 0, sigma, mu3, mu4)
  yf <- bsmimpvol(pf, type, spot, strike, time, rate, 0)
  sum(((yf - y.data) / sy.data)^2, na.rm = TRUE)
}

res <- with(op_vol, {
  optim(
    par = c(0.3, 1, 4), fn = f_obj_csm_min_vol, gr = NULL,
    type, close.underlying, strike, rate, time_to_maturity, bsm_impvol, 1,
    lower = c(1e-2, -Inf, 3), upper = Inf, method = "L-BFGS-B",
    control = list(trace = 3)
  )
})

op_csm <- op_vol |>
  mutate(
    theo_price = csmprice(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      res$par[1], res$par[2], res$par[3]
    ),
    csm_impvol = bsmimpvol(
      theo_price, type, close.underlying, strike, time_to_maturity, rate, 0
    )
  )

op_csm |>
  ggplot(aes(x = strike, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = csm_impvol), colour = "black") +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")

op_csm |>
  mutate(error = close - theo_price) |>
  ggplot(aes(x = strike, y = error)) +
  geom_point()

op_csm |>
  mutate(error = bsm_impvol - csm_impvol) |>
  ggplot(aes(x = strike, y = error)) +
  geom_point()

# high-low weights ----

res <- with(op_vol, {
  optim(
    par = c(0.3, 1, 4), fn = f.optim.csm, gr = NULL,
    type, close.underlying, strike, rate, time_to_maturity, close,
    ifelse(high - low == 0, 1, high - low),
    lower = c(1e-6, -Inf, 3), upper = Inf, method = "L-BFGS-B"
  )
})

op_csm <- op_vol |>
  mutate(
    theo_price = csmprice(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      res$par[1], res$par[2], res$par[3]
    ),
    csm_delta = csmdelta(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      res$par[1], res$par[2], res$par[3]
    ),
    adj_csm_delta = ifelse(csm_delta < 0, 1 + csm_delta, csm_delta),
    csm_impvol = bsmimpvol(
      theo_price, type, close.underlying, strike, time_to_maturity, rate, 0
    )
  )

op_csm |>
  ggplot(aes(x = strike, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = csm_impvol), colour = "black") +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")

op_csm |>
  ggplot(aes(x = adj_delta, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_line(aes(x = adj_csm_delta, y = csm_impvol), colour = "black") +
  geom_vline(xintercept = 0.5, alpha = 0.5, size = 1) +
  labs(
    x = "Delta", y = "Implied Volatility",
    title = str_glue(
      "Implied Volatility {symbol_} - Maturity: {format(maturities[1])}"
    ),
    caption = "Source: B3 (data imported using \U1F4E6 rb3) - wilsonfreitas"
  ) +
  theme(legend.position = "bottom")

op_csm |>
  mutate(error = close - theo_price) |>
  ggplot(aes(x = strike, y = error)) +
  geom_point()

op_csm |>
  mutate(error = bsm_impvol - csm_impvol) |>
  ggplot(aes(x = strike, y = error)) +
  geom_point()

# ----

op1 |>
  ggplot(aes(x = strike, y = theo_price, colour = type, size = volume)) +
  geom_point() +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  facet_wrap(type ~ biz_days) +
  theme(legend.position = "none")

op1 |>
  ggplot(aes(x = adj_delta, y = impvol, colour = type, size = volume)) +
  geom_point() +
  geom_vline(xintercept = 0.5, alpha = 0.5, size = 1) +
  facet_wrap(~biz_days) +
  theme(legend.position = "bottom") +
  labs(
    x = "Delta", y = "Implied Volatility",
    title = str_glue("Equity Options Volatility - {symbol_} {format(refdate)}")
  )

with(op_vol |> filter(type == "Put"), {
  optim(
    par = c(0.3, 1, 4), fn = f.optim.csm, gr = NULL,
    type, close.underlying, strike, rate, time_to_maturity, close, 1
  )
})$par

# high-low weights ----
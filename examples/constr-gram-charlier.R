library(rb3)
library(bizdays)
library(tidyverse)
library(oplib)

# ----

refdate <- getdate("last bizday", Sys.Date(), "Brazil/B3")
ch <- cotahist_get(refdate, "daily")
yc <- yc_get(refdate)
op <- cotahist_equity_options_superset(ch, yc)

symbol_ <- "BBDC4"
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
    vega = bsmvega(
      type, close.underlying, strike, time_to_maturity, rate, 0, bsm_impvol
    ),
    adj_delta = ifelse(str_to_lower(type) == "call", delta, 1 + delta)
  ) |>
  select(
    symbol, volume,
    type, close.underlying, strike, time_to_maturity, rate,
    biz_days, close, high, low, bsm_impvol, vega, delta, adj_delta
  )

# ----

typo <- c("Call")

params <- with(op_vol |> filter(type %in% typo, !is.na(adj_delta)), {
  csm_fit_min_price(
    par = c(0.1, 0, 3),
    type, close.underlying, strike, rate, 0, time_to_maturity, close, 1,
    control = list(trace = 3)
  )
})

strike_rng <- range(op_vol$strike)

gen_data <- tibble(
  type = "Call",
  close.underlying = op_vol$close.underlying[1],
  strike = seq(strike_rng[1], strike_rng[2], length.out = 100),
  time_to_maturity = op_vol$time_to_maturity[1],
  rate = op_vol$rate[1]
)

gen_data <- gen_data |>
  mutate(
    theo_price = csmprice(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      as.numeric(params[1]), as.numeric(params[2]), as.numeric(params[3])
    ),
    csm_impvol = bsmimpvol(
      theo_price, type, close.underlying, strike, time_to_maturity, rate, 0
    ),
    delta = csmdelta(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      as.numeric(params[1]), as.numeric(params[2]), as.numeric(params[3])
    )
  )

op_vol_f <- op_vol |>
  filter(type %in% typo) |>
  mutate(
    theo_price = csmprice(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      as.numeric(params[1]), as.numeric(params[2]), as.numeric(params[3])
    ),
    csm_impvol = bsmimpvol(
      theo_price, type, close.underlying, strike, time_to_maturity, rate, 0
    )
  )

op_vol_f |>
  mutate(error = (close - theo_price) / close) |>
  ggplot(aes(x = strike, y = error, size = volume, colour = type)) +
  geom_point(alpha = 0.5)

op_vol_f |>
  ggplot(aes(x = strike, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(x = strike, y = csm_impvol), data = gen_data, colour = "black") +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")

op_vol_f |>
  ggplot(aes(x = delta, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(x = delta, y = csm_impvol), data = gen_data, colour = "black") +
  theme(legend.position = "none")

op_vol_f |>
  ggplot(aes(x = strike, y = close, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = theo_price), colour = "black") +
  theme(legend.position = "none") +
  facet_wrap(. ~ type)

library(rb3)
library(bizdays)
library(tidyverse)
library(DEoptim)
devtools::load_all()

# ----

refdate <- getdate("last bizday", Sys.Date(), "Brazil/B3")
ch <- cotahist_get(refdate, "daily")
yc <- yc_get(refdate)
op <- cotahist_equity_options_superset(ch, yc)

symbol_ <- "BBAS3"
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
    adj_delta = ifelse(str_to_lower(type) == "call", delta, 1 + delta),
    price_error = sqrt(0.25 * log(high / low) / log(2))
  ) |>
  select(
    symbol, volume,
    type, close.underlying, strike, time_to_maturity, rate,
    biz_days, close, high, low, bsm_impvol, vega, delta, adj_delta, price_error
  )

# min vol error ----

typo <- c("Call")

op_vol <- op_vol |>
  filter(type %in% typo, !is.na(delta))

# BBAS3 1 vencimento 2022-07-27, fica melhor com weights = vega, no ATM, OTM zoa
params <- with(op_vol, {
  csm_fit_min_vol(
    par = c(0.5, 0, 0), type, close.underlying, strike, rate, 0,
    time_to_maturity, bsm_impvol, vega
  )
})

op_vol_f <- op_vol |>
  mutate(
    theo_price = csmprice(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      params[1], params[2], params[3]
    ),
    csm_impvol = bsmimpvol(
      theo_price, type, close.underlying, strike, time_to_maturity, rate, 0
    )
  )

op_vol_f |>
  mutate(error = (bsm_impvol - csm_impvol)) |>
  summarise(se = sum(error ^2))

op_vol_f |>
  mutate(error = (bsm_impvol - csm_impvol) / bsm_impvol) |>
  ggplot(aes(x = adj_delta, y = error, colour = type)) +
  geom_point(alpha = 0.5) +
  facet_wrap(. ~ type)

op_vol_f |>
  ggplot(aes(x = strike, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_point(aes(y = csm_impvol), colour = "black") +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")

# deoptim ----

csm_fit_deoptim_min_vol <- function(par,
                                    type, spot, strike, rate, yield, time,
                                    vol, weights = 1, ...) {
  fn <- function(par) {
    csm_obj_min_vol(
      par, type, spot, strike, rate, yield, time, bsm_impvol, weights
    )
  }
  DEoptim(
    fn = fn, lower = c(1e-1, -40, -10), upper = c(1, 40, 10), ...
  )

  c(
    res$optim$bestmem[1],
    uncons_regionD(res$optim$bestmem[2], res$optim$bestmem[3])
  ) |> unname()
}

params <- with(op_vol, {
  csm_fit_deoptim_min_vol(
    par, type, close.underlying, strike, rate, 0,
    time_to_maturity, bsm_impvol, 1,
    control = DEoptim.control(itermax = 100)
  )
})

op_vol_f <- op_vol |>
  mutate(
    theo_price = csmprice(
      type, close.underlying, strike, time_to_maturity, rate, 0,
      params[1], params[2], params[3]
    ),
    csm_impvol = bsmimpvol(
      theo_price, type, close.underlying, strike, time_to_maturity, rate, 0
    )
  )

op_vol_f |>
  mutate(error = (bsm_impvol - csm_impvol) / bsm_impvol) |>
  ggplot(aes(x = adj_delta, y = error, colour = type)) +
  geom_point(alpha = 0.5) +
  facet_wrap(. ~ type)

op_vol_f |>
  ggplot(aes(x = strike, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_point(aes(y = csm_impvol), colour = "black") +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")
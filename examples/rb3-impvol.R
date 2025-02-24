library(rb3)
library(bizdays)
library(tidyverse)
library(oplib)

refdate <- getdate("last bizday", Sys.Date(), "Brazil/B3")
ch <- cotahist_get(refdate, "daily")
yc <- yc_get(refdate)

symbol_ <- "PETR4"
op1 <- cotahist_options_by_symbol_superset(symbol_, ch, yc)

maturities <- unique(op1$maturity_date) |> sort()
close_underlying <- op1$close.underlying[1]

op_vol <- op1 |>
  filter(maturity_date %in% maturities[2]) |>
  mutate(
    biz_days = bizdays(
      refdate, following(maturity_date, "Brazil/ANBIMA"), "Brazil/ANBIMA"
    ),
    time_to_maturity = biz_days / 252,
    rate = log(1 + r_252),
    impvol = bsmimpvol(
      close, type, close.underlying, strike, time_to_maturity, rate, 0
    ),
    delta = bsmdelta(
      type, close.underlying, strike, time_to_maturity, rate, 0, impvol
    ),
    adj_delta = ifelse(str_to_lower(type) == "call", delta, 1 + delta)
  ) |>
  select(
    symbol, volume,
    type, close.underlying, strike, time_to_maturity, rate, impvol,
    delta, adj_delta, biz_days, volume
  )

with(
  op1 |>
    filter(maturity_date %in% maturities[2]) |>
    mutate(
      biz_days = bizdays(
        refdate, following(maturity_date, "Brazil/ANBIMA"), "Brazil/ANBIMA"
      ),
      time_to_maturity = biz_days / 252,
      rate = log(1 + r_252)
    ),
  bsmimpvol(
    close, type, close.underlying, strike, time_to_maturity, rate, 0
  )
)

op_vol |>
  filter(!is.na(impvol)) |>
  ggplot(aes(x = strike, y = impvol, colour = type, size = volume)) +
  geom_point() +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  facet_wrap(type ~ biz_days) +
  theme(legend.position = "bottom") +
  labs(
    x = "Strike", y = "Implied Volatility",
    title = str_glue("Equity Options Volatility - {symbol_} {format(refdate)}")
  )

op_vol |>
  filter(!is.na(impvol)) |>
  ggplot(aes(x = adj_delta, y = impvol, colour = type, size = volume)) +
  geom_point() +
  geom_vline(xintercept = 0.5, alpha = 0.5, size = 1) +
  facet_wrap(~biz_days) +
  theme(legend.position = "bottom") +
  labs(
    x = "Delta", y = "Implied Volatility",
    title = str_glue("Equity Options Volatility - {symbol_} {format(refdate)}")
  )

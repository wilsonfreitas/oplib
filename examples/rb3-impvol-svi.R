
library(rb3)
library(bizdays)
library(tidyverse)

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

res <- with(op_vol |> filter(!is.na(impvol)), {
  svi_fit(
    impvol ^ 2, close_underlying[1], strike, time_to_maturity[1], rate[1]
  )
})

svi_vol <- with(op_vol |> filter(!is.na(impvol)), {
  F <- close_underlying * exp(rate[1] * time_to_maturity[1])
  strike_range <- range(strike)
  K <- seq(strike_range[1], strike_range[2], length.out = 100)
  var <- svi_var(res$a, res$b, res$m, res$rho, F / K, res$sigma)
  tibble(strike = K, svi_vol = sqrt(var))
})

op_vol |>
  filter(!is.na(impvol)) |>
  ggplot() +
  geom_point(aes(x = strike, y = impvol, colour = type, size = volume)) +
  geom_line(data = svi_vol, aes(x = strike, y = svi_vol)) +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")

op_vol |>
  filter(!is.na(impvol)) |>
  ggplot(aes(x = adj_delta, y = impvol, colour = type, size = volume)) +
  geom_point() +
  geom_vline(xintercept = 0.5, alpha = 0.5, size = 1) +
  facet_wrap(~ biz_days) +
  theme(legend.position = "bottom") +
  labs(
    x = "Delta", y = "Implied Volatility",
    title = str_glue("Equity Options Volatility - {symbol_} {format(refdate)}")
  )
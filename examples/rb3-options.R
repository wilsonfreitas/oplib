library(rb3)
library(bizdays)
library(tidyverse)
devtools::load_all()

refdate <- preceding(Sys.Date() - 1, "Brazil/B3")
ch <- cotahist_get(refdate, "daily")
yc <- yc_get(refdate)
op <- cotahist_equity_options_superset(ch, yc)

symbol_ <- "B3SA3"
op1 <- op |>
  filter(
    symbol.underlying == symbol_
  )

op1 |>
  mutate(
    cost = close / close.underlying,
    moneyness = (close.underlying * (1 + r_252)
    ^ (bizdays(refdate, maturity_date, "Brazil/ANBIMA") / 252)) / strike
  ) |>
  select(symbol, type, moneyness, cost) |>
  arrange(type, moneyness) |>
  View()
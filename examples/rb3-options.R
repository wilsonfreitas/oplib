library(rb3)
library(bizdays)
library(tidyverse)
devtools::load_all()

refdate <- preceding(Sys.Date() - 1, "Brazil/B3")
ch <- cotahist_get(refdate, "daily")
yc <- yc_get(refdate)
op <- cotahist_equity_options_superset(ch, yc)

mats <- unique(op$maturity_date) |> sort()

symbol_ <- "ABEV3"
op |>
  filter(
    symbol.underlying == symbol_,
    maturity_date %in% mats[1:2]
  ) |>
  mutate(
    cost = close / close.underlying,
    moneyness = (close.underlying * (1 + r_252)
    ^ (bizdays(refdate, maturity_date, "Brazil/ANBIMA") / 252)) / strike
  ) |>
  select(symbol, maturity_date, type, moneyness, cost) |>
  arrange(type, moneyness) |>
  ggplot(aes(x = moneyness, y = cost)) +
  geom_point() +
  facet_wrap(maturity_date ~ type)
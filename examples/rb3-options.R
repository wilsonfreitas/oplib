library(rb3)
library(fixedincome)
library(bizdays)
library(tidyverse)
devtools::load_all()

refdate <- getdate("last bizday", Sys.Date(), "Brazil/B3")
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
    rate = spotrate(r_252, "discrete", "business/252", "Brazil/ANBIMA"),
    biz_days = bizdays(refdate, maturity_date, "Brazil/ANBIMA"),
    future = close.underlying * compound(rate, term(biz_days, "days")),
    moneyness = future / strike
  ) |>
  select(symbol, maturity_date, type, moneyness, cost) |>
  arrange(type, moneyness) |>
  ggplot(aes(x = moneyness, y = cost)) +
  geom_point() +
  facet_wrap(maturity_date ~ type)
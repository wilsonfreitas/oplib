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

with(op1, {
  foption(
    symbol = symbol,
    underlying_symbol = symbol.underlying,
    type = type,
    strike = strike,
    maturity_date = maturity_date,
    calendar = "Brazil/B3"
  )
})

theovalue(opts)
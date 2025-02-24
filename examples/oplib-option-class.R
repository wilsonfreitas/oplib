library(rb3)
library(bizdays)
library(tidyverse)
library(oplib)

refdate <- getdate("last bizday", Sys.Date(), "Brazil/B3")
ch <- cotahist_get(refdate, "daily")
yc <- yc_get(refdate)

symbol_ <- "PETR4"

op1 <- cotahist_options_by_symbol_superset(symbol_, ch, yc)

options <- select(op1, symbol, type, strike, maturity_date)

underlying <- function(type) {

}

european_option <- function(type, strike, maturity, calendar, underlying) {

}

# pricing option
# - refdate
#   - business days
#   - time to maturity in years
# - interest rates: risk free, dividends, cost of carry ...
#   - constant rate
#   - interpolate from yield curve
# - spot price

# type, spot, strike, time, rate, yield, sigma
option_price(options, spot)

select(op1, refdate, close.underlying, maturity_date, r_252) |> distinct()


library(rb3)
library(oplib)
library(bizdays)
library(tidyverse)

refdate <- getdate("last bizday", Sys.Date(), "Brazil/B3")
ch <- cotahist_get(refdate, "yearly")
yc <- yc_mget(first_date = as.Date("2022-01-01"), last_date = refdate)

symbol_ <- "PETR4"

eqs <- ch[["HistoricalPrices"]] |>
  filter(.data$cod_negociacao == symbol_) |>
  rb3:::format_equity(TRUE)
eqs_opts <- ch[["HistoricalPrices"]] |>
  filter(.data$tipo_mercado %in% c(70, 80)) |>
  rb3:::format_options(TRUE) |>
  filter(cod_isin == eqs$cod_isin[1])
op <- inner_join(eqs_opts, eqs,
  by = c("refdate", "cod_isin"),
  suffix = c("", ".underlying")
) |>
  select(-c(.data$cod_isin)) |>
  mutate(
    fixing_maturity_date = following(.data$maturity_date, "Brazil/ANBIMA")
  ) |>
  inner_join(yc |> select(.data$refdate, .data$forward_date, .data$r_252),
    by = c("refdate", "fixing_maturity_date" = "forward_date")
  )

op1 <- op |>
  split(op$refdate) |>
  map_dfr(function(df) {
    min_date <- which.min(df$maturity_date)
    df |>
      filter(maturity_date == df$maturity_date[min_date])
  })

mats <- op1$maturity_date |>
  unique() |>
  sort()

op_vol <- op1 |>
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
    )
  ) |>
  select(
    refdate, symbol, volume,
    type, close.underlying, strike, time_to_maturity, rate,
    biz_days, close, high, low, bsm_impvol, delta
  )

op_atm <- op_vol |>
  split(op_vol$refdate) |>
  map_dfr(function(df) {
    df_type <- df[df$type == "Call", ]
    min_idx <- which.min(abs(abs(df_type$delta) - 0.75))
    df[min_idx, ]
  })

op_atm |>
  ggplot(aes(x = refdate, y = bsm_impvol)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = mats, color = "red")

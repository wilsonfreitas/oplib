library(rb3)
library(bizdays)
library(tidyverse)
library(oplib)

# ----

refdate <- getdate("last bizday", Sys.Date(), "Brazil/B3")
ch <- cotahist_get(refdate, "daily")
yc <- yc_get(refdate)
op <- cotahist_equity_options_superset(ch, yc)

symbol_ <- "PETR4"
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
    par = c(0.1, 0, 0),
    type, close.underlying, strike, rate, 0, time_to_maturity, close, 1,
    control = list(trace = 3)
  )
})

print(params)

# search space ----

obj_space <- function(df, typo) {
  function(par) {
    par <- c(params[1], par)
    with(df |> filter(type %in% typo, !is.na(adj_delta)), {
      csm_obj_min_price(
        par, type, close.underlying, strike, rate, 0, time_to_maturity, close, 1
      )
    })
  }
}

ff <- obj_space(op_vol, "Call")

comb <- expand.grid(
  mu3 = seq(-5, 5, length.out = 50),
  mu4 = seq(-5, 5, length.out = 50)
)

z <- sapply(seq_len(dim(comb)[1]), function(i) {
  parms <- as.numeric(comb[i, ])
  ff(parms)
})

comb$ll <- z
h <- hist(z, 20, plot = FALSE)
comb$llgroups <- cut(z, h$breaks)

t_vars <- uncons_regionD(comb$mu3, comb$mu4)
comb$mu3_t <- t_vars[, "mu3"]
comb$mu4_t <- t_vars[, "mu4"]

theme_set(theme_bw(base_size = 9))
ggplot(comb, aes(x = mu4, y = mu3, fill = ll)) +
  geom_tile() +
  labs(x = "kurtosis", y = "skewness", fill = "llh")

print(comb[which.min(comb$ll), ])

ggplot(comb, aes(x = mu4_t, y = mu3_t, colour = llgroups)) +
  geom_point() +
  labs(x = "kurtosis", y = "skewness", fill = "llh")

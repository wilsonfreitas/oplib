library(rb3)
library(bizdays)
library(tidyverse)
library(oplib)

# ----

refdate <- getdate("last bizday", Sys.Date(), "Brazil/B3")
ch <- cotahist_get(refdate, "daily")
yc <- yc_get(refdate)
op <- cotahist_equity_options_superset(ch, yc)

symbol_ <- "ABEV3"
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

grad_obj_csm <- function(par, type, spot, strike, rate, time, y.data, sy.data) {
  par_ <- uncons_regionD(par[2], par[3])
  par <- c(par[1], par_[1], par_[2])

  #numeric vega for Corrado-Su mod
  dcs.dsig <- csmvega(
    type, spot, strike, time, rate, 0,
    par[1], par[2], par[3]
  )

  #other derivatives for the gradient
  csmw_ <- csmw(par[1], time, par[2], par[3])
  dmod <- csmd(spot, strike, time, rate, 0, par[1], csmw_)
  q3 <- csmq3(spot, par[1], time, dmod, csmw_)
  q4 <- csmq4(spot, par[1], time, dmod, csmw_)

  premium <- csmprice(type, spot, strike, time, rate, 0, par[1], par[2], par[3])
  grad1 <- sum(2 * dcs.dsig * (premium -  y.data) / (sy.data ^ 2))
  grad2 <- sum(2 * q3 * (premium -  y.data) / (sy.data ^ 2))
  grad3 <- sum(2 * q4 * (premium -  y.data) / (sy.data ^ 2))

  ret <- c(grad1, grad2, grad3)
  ret
}

grad2_obj_csm <- function(par,
                          type, spot, strike, rate, time, y.data, sy.data) {
  bump <- 1e-4

  p <- c(par[1] + bump, par[2], par[3])
  du <- f_obj_csm(p, type, spot, strike, rate, time, y.data, sy.data)
  p <- c(par[1] - bump, par[2], par[3])
  dd <- f_obj_csm(p, type, spot, strike, rate, time, y.data, sy.data)
  grad1 <- (du - dd) / (2 * bump)

  p <- c(par[1], par[2] + bump, par[3])
  du <- f_obj_csm(p, type, spot, strike, rate, time, y.data, sy.data)
  p <- c(par[1], par[2] - bump, par[3])
  dd <- f_obj_csm(p, type, spot, strike, rate, time, y.data, sy.data)
  grad2 <- (du - dd) / (2 * bump)

  p <- c(par[1], par[2], par[3] + bump)
  du <- f_obj_csm(p, type, spot, strike, rate, time, y.data, sy.data)
  p <- c(par[1], par[2], par[3] - bump)
  dd <- f_obj_csm(p, type, spot, strike, rate, time, y.data, sy.data)
  grad3 <- (du - dd) / (2 * bump)

  c(grad1, grad2, grad3)
}

f_obj_csm <- function(par, type, spot, strike, rate, time, y.data, w.data) {
  sigma <- par[1]
  params <- uncons_regionD(par[2], par[3])
  par_ <- c(par[1], params)
  yf <- csmprice(type, spot, strike, time, rate, 0, par_[1], par_[2], par_[3])
  ret <- sum(((yf - y.data) * w.data)^2)
  ret
}

# ----

typo <- c("Put")
# res1 <- with(op_vol |> filter(type %in% typo, !is.na(adj_delta)), {
#   optim(
#     par = c(0.1, 0, 3), fn = f_obj_csm, gr = grad_obj_csm,
#     method = "L-BFGS-B",
#     type, close.underlying, strike, rate, time_to_maturity, close, vega,
#     control = list(trace = 3)
#   )
# })

# params <- uncons_regionD(res$par[2], res$par[3])

params <- with(op_vol |> filter(type %in% typo, !is.na(adj_delta)), {
  csm_fit_min_price(
    par = c(0.1, 0, 3),
    type, close.underlying, strike, rate, 0, time_to_maturity, close, vega,
    control = list(trace = 3)
  )
})

op_vol_f <- op_vol |>
  filter(type %in% typo) |>
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
  mutate(error = (close - theo_price) / close) |>
  ggplot(aes(x = strike, y = error, size = volume, colour = type)) +
  geom_point(alpha = 0.5)

op_vol_f |>
  ggplot(aes(x = strike, y = bsm_impvol, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = csm_impvol), colour = "black") +
  geom_vline(xintercept = close_underlying, alpha = 0.5, size = 1) +
  theme(legend.position = "none")

op_vol_f |>
  ggplot(aes(x = strike, y = close, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = theo_price), colour = "black") +
  theme(legend.position = "none") +
  facet_wrap(. ~ type)

# ----

obj_space <- function(df, typo) {
  function(par) {
    par <- c(0.40, par)
    with(df |> filter(type %in% typo, !is.na(adj_delta)), {
      f_obj_csm(
        par, type, close.underlying, strike, rate, time_to_maturity, close, 1
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
h <- hist(z, 10, plot = FALSE)
comb$llgroups <- cut(z, h$breaks)

comb[which.min(z), ]

theme_set(theme_bw(base_size = 9))
ggplot(comb, aes(x = mu4, y = mu3, fill = z)) +
  geom_tile() +
  labs(x = "kurtosis", y = "skewness", fill = "llh")

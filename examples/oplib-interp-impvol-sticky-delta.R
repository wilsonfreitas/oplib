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

op_data <- op1 |>
  filter(maturity_date %in% maturities[2]) |>
  mutate(
    biz_days = bizdays(
      refdate, following(maturity_date, "Brazil/ANBIMA"), "Brazil/ANBIMA"
    ),
    time_to_maturity = biz_days / 252,
    rate = log(1 + r_252),
    forward = close.underlying * exp(rate * time_to_maturity),
    impvol = bsmimpvol(
      close, type, close.underlying, strike, time_to_maturity, rate, 0
    ),
    delta = bsmdelta(
      type, close.underlying, strike, time_to_maturity, rate, 0, impvol
    )
  )

with(op_data, close >= forward)

# smile with all options ----
smile <- with(
  op_data,
  {
    sigma <- bsmimpvol(
      close, type, close.underlying, strike, time_to_maturity, rate, 0
    )
    delta <- bsmdelta(
      type, close.underlying, strike, time_to_maturity, rate, 0, sigma
    )

    delta <- ifelse(type == "Put", 1 + delta, delta)

    ix <- order(delta)
    delta <- delta[ix]
    sigma <- sigma[ix]

    xy.coords(
      x = delta[!is.na(sigma)],
      y = sigma[!is.na(sigma)],
      xlab = "Delta", ylab = "Volatility"
    )
  }
)

plot(smile, type = "b")


# smile with call for delta <= ½ and put for delta >= ½ ----

smile <- with(
  op_data,
  {
    sigma <- bsmimpvol(
      close, type, close.underlying, strike, time_to_maturity, rate, 0
    )
    delta <- bsmdelta(
      type, close.underlying, strike, time_to_maturity, rate, 0, sigma
    )

    delta <- ifelse(type == "Put", 1 + delta, delta)

    ix_call <- (type == "Call" & delta <= 0.5)
    ix_put <- (type == "Put" & delta >= 0.5)
    delta <- c(delta[ix_call], delta[ix_put])
    sigma <- c(sigma[ix_call], sigma[ix_put])
    ix <- order(delta)
    delta <- delta[ix]
    sigma <- sigma[ix]

    xy.coords(
      x = delta[!is.na(sigma)],
      y = sigma[!is.na(sigma)],
      xlab = "Delta", ylab = "Volatility"
    )
  }
)

nat <- splinefun(smile$x, smile$y, "f")
monoh <- splinefun(smile$x, smile$y, "m")
plot(smile, type = "b")
x_ <- setdiff(sort(c(seq(0, 1, 0.005), smile$x)), NA)
lines(x_, nat(x_), col = "blue")
lines(x_, monoh(x_), col = "darkgreen")

# delta interpolation ----

smilefun <- monoh

new_close_underlying <- op_data$close.underlying[1] * 1.01

df <- with(op_data, {
  delta <- 0.5
  new_sigma <- smilefun(delta)
  new_delta <- bsmdelta(
    type, new_close_underlying, strike, time_to_maturity, rate, 0, new_sigma
  )
  new_delta <- ifelse(type == "Call", new_delta, 1 + new_delta)
  err <- abs(new_delta - delta)
  count <- 0
  err_change <- 1
  while (any(err > 1e-4) & count < 100 & any(err_change > 1e-6)) {
    delta <- new_delta
    new_sigma <- smilefun(delta)
    new_delta <- bsmdelta(
      type, new_close_underlying, strike, time_to_maturity, rate, 0, new_sigma
    )
    new_delta <- ifelse(type == "Call", new_delta, 1 + new_delta)
    err_ <- abs(new_delta - delta)
    err_change <- abs(err - err_)
    err <- err_
    count <- count + 1
  }
  data.frame(
    close = close,
    price = bsmprice(
      type, new_close_underlying, strike, time_to_maturity, rate, 0, new_sigma
    ),
    delta = bsmdelta(
      type, new_close_underlying, strike, time_to_maturity, rate, 0, new_sigma
    ),
    sigma = new_sigma
  )
})

df

# rb3 create volatility surface ----

build_options <- function(type, strike, time) {
  ops <- data.frame(type, strike, time)
  ops["type"] <- factor(tolower(sapply(ops["type"], as.character)),
    levels = c("call", "put"), labels = c("call", "put")
  )
  ops
}

create_volatility_surface <- function(price, type, spot, strike, time, rate, yield) {
  op1 <- build_options(type, strike, time)
  op1$price <- price
  op1$spot <- spot
  op1$rate <- rate
  op1$yield <- yield
  split(op1, op1$time) |>
    lapply(function(df) {
      df |>
        with({
          sigma <- try(bsmimpvol(
            price, type, spot, strike, time, rate, yield
          ), silent = TRUE)
          if (is(sigma, "try-error")) {
            return(NULL)
          }
          delta <- bsmdelta(
            type, spot, strike, time, rate, yield, sigma
          )

          delta <- ifelse(type == "put", 1 + delta, delta)

          ix_call <- (type == "call" & delta <= 0.5)
          ix_put <- (type == "put" & delta >= 0.5)
          strike <- c(strike[ix_call], strike[ix_put])
          delta <- c(delta[ix_call], delta[ix_put])
          sigma <- c(sigma[ix_call], sigma[ix_put])
          ix <- order(delta)
          delta <- delta[ix]
          sigma <- sigma[ix]
          strike <- strike[ix]

          ix <- !is.na(sigma)
          data.frame(
            time = time[1],
            strike = strike[ix],
            delta = delta[ix],
            volatility = sigma[ix]
          )
        })
    }) |>
    bind_rows() |>
    arrange(time, strike)
}

volsurf <- with(op1, {
  biz_days <- bizdays(
    refdate, following(maturity_date, "Brazil/ANBIMA"), "Brazil/ANBIMA"
  )
  time_to_maturity <- biz_days / 252
  rate <- log(1 + r_252)
  create_volatility_surface(close, type, close.underlying, strike, time_to_maturity, rate, 0)
})

# plotting volatility surface ----

library(plotly)

plot_ly(volsurf,
  x = ~strike, y = ~time, z = ~volatility,
  type = "scatter3d", size = 10
)

plot_ly(volsurf,
  x = ~delta, y = ~time, z = ~volatility,
  type = "scatter3d", size = 10
)

# creating function to interpolate ----

interp_smile_sticky_strike <- function(type, new_spot, strike, time, rate, yield, stickystrikefun) {
  df <- data.frame(type, new_spot, strike, time, rate, yield)

  with(df, {
    delta <- 0.5
    new_sigma <- stickystrikefun(delta)
    new_delta <- bsmdelta(
      type, new_spot, strike, time, rate, yield, new_sigma
    )
    new_delta <- ifelse(type == "Call", new_delta, 1 + new_delta)
    err <- abs(new_delta - delta)
    count <- 0
    # err_change <- 1
    while (any(err > 1e-4) & count < 100) { #  & any(err_change > 1e-6)
      delta <- new_delta
      new_sigma <- stickystrikefun(delta)
      new_delta <- bsmdelta(
        type, new_spot, strike, time, rate, yield, new_sigma
      )
      new_delta <- ifelse(type == "Call", new_delta, 1 + new_delta)
      err_ <- abs(new_delta - delta)
      # err_change <- abs(err - err_)
      err <- err_
      count <- count + 1
    }
    new_sigma
  })
}

new_close_underlying <- op_data$close.underlying[1] * 1.01
op_data$sigma <- with(op_data, {
  interp_smile_sticky_strike(type, new_close_underlying, strike, time_to_maturity, rate, 0, smilefun)
})

op_data$new_delta <- with(op_data, {
  bsmdelta(type, new_close_underlying, strike, time_to_maturity, rate, 0, sigma)
})

op_data |>
  filter(!is.na(impvol)) |>
  ggplot(aes(x = strike)) +
  geom_point(aes(y = impvol, colour = type)) +
  geom_point(aes(y = sigma)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Strike", y = "Implied Volatility",
    title = str_glue("Equity Options Volatility - {symbol_} {format(refdate)}")
  )

op_data |>
  filter(!is.na(impvol)) |>
  ggplot(aes(x = strike)) +
  geom_point(aes(y = delta, colour = type)) +
  geom_point(aes(y = new_delta)) +
  theme(legend.position = "bottom") +
  labs(
    x = "Strike", y = "Delta",
    title = str_glue("Equity Options Volatility - {symbol_} {format(refdate)}")
  )

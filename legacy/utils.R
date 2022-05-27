# calculate historical volatility for the period in
# ohlc_data using log returns
calc_hist_vol <- function(ohlc_data, year_biz_days) {
  # ln(s_i+1 / s_i)
  close_vals <- ohlc_data[, 1]
  x <- sapply(2:nrow(ohlc_data), function(idx) {
    log((close_vals[[idx]]) / (close_vals[[idx - 1]]))
  })
  x_mean <- mean(x)
  x_sum <- sum((x - x_mean)^2)
  sqrt(x_sum * year_biz_days / (nrow(ohlc_data) - 2)) # no daily return for the first value
}

# Build volatility table for a given ohlc_data.
# Periods will be separated according to the number of
# months available in the data.
# Each row of the volatility table indicates when a period
# starts, and each column indicates the period length for
# which historical volatilities will be calculated.
build_vol_table <- function(ohlc_data) {
  start_date <- as.Date(min(zoo::index(ohlc_data)))
  end_date <- as.Date(max(zoo::index(ohlc_data)))
  sample_size <- length(zoo::index(ohlc_data))
  suppressWarnings(nmonth <- abs(interval(start_date, end_date)) %/% months(1))
  days_by_period <- sample_size %/% nmonth
  ymd_periods <- zoo::as.Date(sapply(0:(nmonth - 1), function(m) {
    zoo::index(ohlc_data[1 + m * days_by_period])
  }))
  year_days <- bizdays(
    from = paste(year(start_date), "01", "01", sep = "-"),
    to = paste(year(start_date), "12", "31", sep = "-"),
    cal = "Brazil/ANBIMA"
  )
  if (length(ymd_periods) > nmonth) ymd_periods <- ymd_periods[-length(ymd_periods)]

  vol_table <- data.frame(row.names = ymd_periods)
  for (m in 1:nmonth) {
    vol_table[, as.character(m * days_by_period)] <- sapply(ymd_periods, function(st_date) {
      st_date_index <- match(as.Date(st_date), zoo::index(ohlc_data))
      ed_date_index <- st_date_index + days_by_period * m
      if (ed_date_index <= (nmonth) * days_by_period + 1 && ed_date_index - st_date_index > 2) {
        if (ed_date_index > sample_size) ed_date_index <- sample_size
        tryCatch(
          {
            calc_hist_vol(ohlc_data[st_date_index:ed_date_index], year_days)
          },
          error = function(e) {
            browser()
          }
        )
      } else {
        NA
      }
    })
  }

  vol_table
}

# Calculate a confidence interval for a given alpha
confidence_interval <- function(alpha, n) {
  lower <- sqrt((n - 1) / qchisq(p = alpha / 2, df = n - 1))
  higher <- sqrt((n - 1) / qchisq(p = (1 - alpha / 2), df = n - 1))
  c(lower, higher)
}

# Calculate implied volatilities for all options traded on ref_d
# Risk free rate is given by PRE curve at ref_d
# Implied volatility is calculated using bsm
calc_impl_vol <- function(options_data, ref_d) {
  dividendYield <- 0
  options <- options_data[!is.na(options_data$close), ]
  options <- options[options$close != 0, ]
  options <- options[options$ref_date == ref_d, ]
  options <- transform(options, bznsdtexp = bizdays(ref_date, maturity_date, cal = "Brazil/ANBIMA"))
  options <- transform(options, dtexp = as.numeric(as.Date(maturity_date) - as.Date(ref_date)))
  rate_curve <- get_curve("PRE", ref_d)

  if (nrow(options) >= 1) {
    impl_vol <- apply(options, 1, function(option) {
      close <- as.numeric(option[["close"]])
      type <- option[["type"]]
      spot <- as.numeric(option[["spot"]])
      strike <- as.numeric(option[["strike"]])
      dtexp <- as.numeric(option[["dtexp"]])
      bzdtexp <- as.numeric(option[["bznsdtexp"]])
      riskFreeRate <- rate_curve[rate_curve$terms == dtexp, ]$value
      i <- dtexp - 5
      while (length(riskFreeRate) == 0) {
        riskFreeRate <- rate_curve[rate_curve$terms == i, ]$value
        i <- i + 1
      }
      ivol <- bsmimpvol(close, type, spot, strike, dtexp / 365, riskFreeRate / 100, dividendYield)
      data.frame(bz_days_to_exp = bzdtexp, days_to_exp = dtexp, value = ivol, type = type)
    })
    impl_vol <- do.call(rbind, impl_vol)
  } else {
    impl_vol <- NULL
  }
  impl_vol
}


#' Volatility cone graph generator
#'
#' Generates a graph showing a volatility cone for the stock given
#' by ticker, calculated with data starting from historical_sample_start_date
#' to the most recent available date. It also calculates implied volatilities
#' for that stock's options that were traded on ref_date and plots it on
#' the same graph.
#' The volatility estimator for the cone is annualized standard deviation
#' of log returns. Its confidence interval is for a 5% alpha.
#' The implied volatilities are calculated using the bissection method
#' and the Black-Scholes model.
#'
#' @param ticker Stock to be graphed
#' @param ref_date Reference date for the options and impl. vol calculation
#' @param historical_sample_start_date Start date for the historical data to be used
#' to calculate the volatility cone
#'
#' @return
#' ggplot object containing the graph with the estimated volatility cone, its
#' confidence intervals and implied volatilities for options traded on ref_date
#'
#'
#' @export
gen_vol_cone <- function(ticker, ref_date, historical_sample_start_date, historical_sample_end_date = today()) {
  options <- get_options(ticker, src = "stock", from = ref_date, to = ref_date)
  vals <- get_stocks(ticker, src = "bloomberg", from = historical_sample_start_date, to = historical_sample_end_date)

  vt <- build_vol_table(vals)

  impl_vols <- calc_impl_vol(options, ref_date)
  if (!is.null(impl_vols)) {
    impl_vols$value <- impl_vols$value * 100 # adjust to percentage values
    # remove entries for which implied volatilities couldn't be found
    # with the bissection method
    impl_vols <- impl_vols[impl_vols$value != 0.000001, ]
  }

  # Adjust historical volatility data and create confidence interval
  # with a 5% alpha
  max_values <- apply(vt[, names(vt)], 2, max, na.rm = TRUE) * 100
  mean_values <- apply(vt[, names(vt)], 2, mean, na.rm = TRUE) * 100
  min_values <- apply(vt[, names(vt)], 2, min, na.rm = TRUE) * 100
  ci <- confidence_interval(0.05, length(vals))
  maxcilow <- c(max_values * ci[1])
  maxcihi <- c(max_values * ci[2])
  meancilow <- c(mean_values * ci[1])
  meancihi <- c(mean_values * ci[2])
  mincilow <- c(min_values * ci[1])
  mincihi <- c(min_values * ci[2])

  cone <- melt(data.frame(
    days_to_exp = as.numeric(names(vt)),
    max_values = max_values,
    maxcilow = maxcilow,
    maxcihi = maxcihi,
    mean_values = mean_values,
    meancilow = meancilow,
    meancihi = meancihi,
    min_values = min_values,
    mincilow = mincilow,
    mincihi = mincihi
  ),
  id.vars = "days_to_exp"
  )

  len <- length(max_values)
  cone$cgroup <- c(rep("Maximum", times = 3 * len), rep("Mean", times = 3 * len), rep("Minimum", times = 3 * len))
  cone$type <- rep(c(rep("Estimate", times = len), rep("Confidence Interval", times = 2 * len)), times = 3)
  cone$type_detail <- rep(c(rep("Estimate", times = len), rep("Confidence Interval lo", times = len), rep("Confidence Interval hi", times = len)), times = 3)

  title <- paste(ticker, "estimated vol. cone with impl. vols for", ref_date)
  subtitle <- paste("Historical volatility calculated between", historical_sample_start_date, "and", historical_sample_end_date)

  plot <- ggplot(cone, aes(
    x = days_to_exp, y = value,
    group = interaction(cgroup, type_detail),
    colour = cgroup, linetype = type
  )) +
    geom_line() +
    geom_point(data = subset(cone, type == "Estimate")) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    scale_x_continuous("Days to Expiry", breaks = c(0, as.numeric(names(vt)))) +
    scale_y_continuous("Vol. (%)") +
    coord_cartesian(ylim = c(max(c(min(mincihi - 5), 0)), max(maxcilow) + 5)) +
    labs(shape = "Option type", linetype = "Line type", colour = "Estimated historical vols.", title = title, subtitle = subtitle) +
    theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(impl_vols)) {
    plot <- plot + geom_point(data = impl_vols, aes(x = bz_days_to_exp, y = value, shape = type), inherit.aes = FALSE)
  }
  plot
}
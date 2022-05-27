#' Computes Q_3 from Corrado-Su model.
#'
#' Implementation of Q_3 formula from Corrado-Su model. Q_3 measures the effects
#' of nonnormal skewness on the option price given by the model.
#'
#' @param spot current stock price.
#' @param sigma volatility of the stock price.
#' @param time time to option expiration in years.
#' @param d1 auxiliary value, same from Black-Scholes-Merton model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with Q_3.
#'
#' @seealso
#' \code{\link{csq4}} for Q_4 formula.\\
#' \code{\link{csprice}} for Corrado-Su option pricing formula.
#'
#' @examples
#' csq3(10, .3, 1.5, 0:100 / 100)
#' spot <- 25
#' sigma <- 0.2
#' time <- 0.5
#' strike <- 26
#' rate <- 0.12
#' yield <- 0
#' d1 <- bsmd1(spot, strike, time, rate, yield, sigma)
#' csq3(spot, sigma, time, d1)
#'
#' @export
csq3 <- function(spot, sigma, time, d1) {
  q3 <- spot * sigma * sqrt(time) / 6
  q3 * ((2 * sigma * sqrt(time) - d1) * dnorm(d1) + sigma * sigma * time *
    pnorm(d1))
}

#' Computes Q_4 from Corrado-Su model.
#'
#' Implementation of Q_4 formula from Corrado-Su model. Q_4 measures the effects
#' of nonnormal kurtosis on the option price given by the model.
#'
#' @param spot current stock price.
#' @param sigma volatility of the stock price.
#' @param time time to option expiration in years.
#' @param d1 auxiliary value, same from Black-Scholes-Merton model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with Q_4.
#'
#' @seealso
#' \code{\link{csq3}} for Q_3 formula.\\
#' \code{\link{csprice}} for Corrado-Su option pricing formula.
#'
#' @examples
#' csq4(10, .3, 1.5, 0:100 / 100)
#' spot <- 25
#' sigma <- 0.2
#' time <- 0.5
#' strike <- 26
#' rate <- 0.12
#' yield <- 0
#' d1 <- bsmd1(spot, strike, time, rate, yield, sigma)
#' csq4(spot, sigma, time, d1)
#'
#' @export
csq4 <- function(spot, sigma, time, d1) {
  q4 <- (d1^2 - 1 - 3 * sigma * sqrt(time) * (d1 - sigma * sqrt(time))) *
    dnorm(d1)
  q4 <- q4 + (sigma * sigma * sigma * time^(1.5) * pnorm(d1))
  q4 * spot * sigma * sqrt(time) / 24
}

#' Computes option price with Corrado-Su model.
#'
#' Implementation of Corrado-Su formula for option pricing.
#'
#' @param type 'call' for call option, any other value for put.
#' @param spot current stock price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param yield the dividends that are expected to be paid.
#' @param sigma volatility of the stock price.
#' @param mu3 skewness.
#' @param mu4 kurtosis.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with option prices.
#'
#' @seealso
#' \code{\link{csq3}} for Q_3 formula from Corrado-Su model.\\
#' \code{\link{csq4}} for Q_4 formula from Corrado-Su model.
#'
#' @examples
#' spot <- 25
#' sigma <- 0.2
#' time <- 0.5
#' strike <- 26
#' rate <- 0.12
#' yield <- 0
#' csprice("call", spot, strike, time, rate, yield, sigma, 0, 3) ==
#'   bsmprice("call", spot, strike, time, rate, yield, sigma)
#' csprice("call", spot, strike, time, rate, yield, sigma, 0:100 / 100, 3)
#' csprice("call", spot, strike, time, rate, yield, sigma, 0.5, 0:1000 / 100)
#'
#' @export
csprice <- function(type, spot, strike, time, rate, yield, sigma, mu3, mu4) {
  c_bs <- bsmprice("call", spot, strike, time, rate, yield, sigma)
  d1 <- bsmd1(spot, strike, time, rate, yield, sigma)
  q3 <- csq3(spot, sigma, time, d1)
  q4 <- csq4(spot, sigma, time, d1)
  callp <- c_bs + mu3 * q3 + (mu4 - 3) * q4
  df <- data.frame(callp, type, spot, yield, time, rate, strike)
  with(
    df,
    ifelse(
      tolower(as.character(type)) == "call",
      callp,
      callp - spot * exp(-yield * time) + strike * exp(-rate * time)
    )
  )
}

#' Computes Q_3 from Corrado-Su model options on futures.
#'
#' Implementation of Q_3 formula from Corrado-Su model for futures. Q_3 measures
#' the effects of nonnormal skewness on the option price given by the model.
#'
#' @param future current future price.
#' @param sigma volatility of the future price.
#' @param time time to option expiration in years.
#' @param d1 auxiliary value, same from Black-76 model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with Q_3.
#'
#' @seealso
#' \code{\link{blackcsq4}} for Q_4 formula for futures.\\
#' \code{\link{blackcsprice}} for Corrado-Su option pricing formula for futures.
#'
#' @examples
#' blackcsq3(10, .3, 1.5, 0:100 / 100)
#' future <- 25
#' sigma <- 0.2
#' time <- 0.5
#' strike <- 26
#' rate <- 0.12
#' d1 <- blackd1(future, strike, time, rate, sigma)
#' blackcsq3(future, sigma, time, d1)
#'
#' @export
blackcsq3 <- function(future, sigma, time, d1) {
  q3 <- future * sigma * sqrt(time) / 6
  q3 * ((2 * sigma * sqrt(time) - d1) * dnorm(d1) + sigma * sigma * time *
    pnorm(d1))
}


#' Computes Q_4 from Corrado-Su model options on futures.
#'
#' Implementation of Q_4 formula from Corrado-Su model for futures. Q_4 measures
#' the effects of nonnormal skewness on the option price given by the model.
#'
#' @param future current future price.
#' @param sigma volatility of the future price.
#' @param time time to option expiration in years.
#' @param d1 auxiliary value, same from Black-76 model.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with Q_4.
#'
#' @seealso
#' \code{\link{blackcsq3}} for Q_3 formula for futures.\\
#' \code{\link{blackcsprice}} for Corrado-Su option pricing formula for futures.
#'
#' @examples
#' blackcsq4(10, .3, 1.5, 0:100 / 100)
#' future <- 25
#' sigma <- 0.2
#' time <- 0.5
#' strike <- 26
#' rate <- 0.12
#' d1 <- blackd1(future, strike, time, rate, sigma)
#' blackcsq4(future, sigma, time, d1)
#'
#' @export
blackcsq4 <- function(future, sigma, time, d1) {
  q4 <- (d1^2 - 1 - 3 * sigma * sqrt(time) * (d1 - sigma * sqrt(time))) *
    dnorm(d1)
  q4 <- q4 + (sigma * sigma * sigma * time^(1.5) * pnorm(d1))
  q4 * future * sigma * sqrt(time) / 24
}

#' Computes call option price with Corrado-Su model for futures.
#'
#' Implementation of Corrado-Su formula for option pricing for futures.
#'
#' @param type 'call' for call option, any other value for put.
#' @param future current future price.
#' @param strike the strike price.
#' @param time time to option expiration in years.
#' @param rate the risk-free interest rate.
#' @param sigma volatility of the future price.
#' @param mu3 skewness.
#' @param mu4 kurtosis.
#'
#' @section Recycle rule:
#'
#' These arguments handle the recycle rule so vectors can be provided
#' and once those vectors differs in length the recycle rule is applied.
#'
#' @return
#' Vector with option prices.
#'
#' @seealso
#' \code{\link{blackcsq3}} for Q_3 formula from Corrado-Su model for futures.\\
#' \code{\link{blackcsq4}} for Q_4 formula from Corrado-Su model for futures.
#'
#' @examples
#' spot <- 25
#' sigma <- 0.2
#' time <- 0.5
#' strike <- 26
#' rate <- 0.12
#' blackcsprice("call", spot, strike, time, rate, sigma, 0, 3) ==
#'   blackprice("call", spot, strike, time, rate, sigma)
#' blackcsprice("call", spot, strike, time, rate, sigma, 0:100 / 100, 3)
#' blackcsprice("call", spot, strike, time, rate, sigma, 0.5, 0:1000 / 100)
#'
#' @export
blackcsprice <- function(type, future, strike, time, rate, sigma, mu3, mu4) {
  c_bs <- blackprice(type, future, strike, time, rate, sigma)
  d1 <- (log(future / strike) + (sigma * sigma / 2) * time) /
    (sigma * sqrt(time))
  q3 <- csq3(future, sigma, time, d1)
  q4 <- csq4(future, sigma, time, d1)
  callp <- c_bs + mu3 * q3 + (mu4 - 3) * q4
  ifelse(
    tolower(type) == "call",
    callp,
    callp - future + strike * exp(-rate * time)
  )
}

#' Corrado-Su model fit for futures.
#'
#' Finds and returns Corrado-Su model parameters sigma, mu3 (skewness) and
#'  mu4 (kurtosis)
#' for a given set of options data.
#'
#' @param prices market prices for options
#' @param types 'call' for call option, any other value for put.
#' @param futures current future price.
#' @param strikes the strike price.
#' @param terms time to option expiration in years.
#' @param rates the risk-free interest rate.
#' @param initial_guess optional initial guess for optimization
#'
#' @return
#' List with model parameters sigma, mu3 and mu4
#'
#' @export
blackcs_fit <- function(prices, types, futures, strikes, terms, rates,
                        initial_guess = NULL) {
  lo <- c(0, -1, 1)
  hi <- c(10, 1, 5)

  if (is.null(initial_guess)) {
    initial_guess <- c(CS.SIGMA, CS.SKEW, CS.KURT)

    fit <- nloptr::nloptr(
      x0 = initial_guess,
      eval_f = blackcs_minsq(prices, types, futures, strikes, terms, rates),
      lb = lo,
      ub = hi,
      opts = list(
        algorithm = "NLOPT_GN_ISRES",
        print_level = 0,
        xtol_rel = 1e-4,
        maxeval = 1000
      )
    )

    initial_guess <- fit$solution
  }

  fit <- nloptr::nloptr(
    x0 = initial_guess,
    eval_f = blackcs_minsq(prices, types, futures, strikes, terms, rates),
    lb = lo,
    ub = hi,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      print_level = 0,
      xtol_rel = 1e-4,
      maxeval = 1000
    )
  )

  curr_param <- fit$solution

  list(sigma = curr_param[1], mu3 = curr_param[2], mu4 = curr_param[3])
}

blackcs_minsq <- function(prices, types, futures, strikes, terms, rates) {
  function(x) {
    sigma <- x[1]
    mu3 <- x[2]
    mu4 <- x[3]
    csp <- blackcsprice(types, futures, strikes, terms, rates, sigma, mu3, mu4)
    if (any(is.nan(csp))) csp <- Inf
    if (any(csp < 0)) csp[csp < 0] <- 0
    sum((prices - csp)^2)
  }
}


# used to estimate de implied volatility
# TODO: adaptar para corrado-su
cs.sigma.obj.func <- function(sigma, fut, strike, rate, time, mu3, mu4, price) {
  csprice(fut, strike, sigma, rate, time, mu3, mu4) - price
}

# TODO: adaptar para corrado-su -- por enquanto aplica apenas para calls
calc.sigma.cs <- function(batdf) {
  n <- length(batdf$Date)
  sigma <- vector(mode = "numeric", length = n)
  for (i in seq_len(batdf$Date)) {
    root <- try(uniroot(
      f = cs.sigma.obj.func, interval = c(-0.1, 1), tol = 1e-15,
      fut = batdf$Fut[i], strike = batdf$Strike[i],
      rate = batdf$Rate[i], time = batdf$TBD[i],
      mu3 = batdf$mu3[i], mu4 = batdf$mu4[i],
      price = batdf$Opt.Price[i]
    ))
    if (inherits(root, "try-error")) {
      print(root)
      print(paste("Inputs:", batdf$Fut[i],
        strike = batdf$Strike[i],
        batdf$Rate[i], batdf$TBD[i], batdf$mu3[i], batdf$mu4[i],
        batdf$Opt.Price[i], sep = " "
      ))
      print(paste("i=", i))
      root <- list(root = -1)
    }
    sigma[i] <- root$root
  }
  return(sigma)
}

# Q14  56256
# V14  57151
# Z14  58056
# G15  58897
# M15  60635
# Z15  63594




#' Black Corrado-Su std. Deltas
#'
#' Finds and returns a dataframe with the strikes and implied
#' volatilites for the provided set of deltas or for the
#' default deltas.
#'
#' @param type Option type, "call" for calls or anything else for puts
#' @param future underlying future price at a given date
#' @param time time in years from a given date to a maturity date
#' @param rate interest free rate
#' @param sigma volatility of the stock price
#' @param mu3 skewness
#' @param mu4 kurtosis
#' @param deltas deltas for which strikes and vols should be calculated, 1-99. Defaults to usual
#' delta values.
#'
#' @return
#' Dataframe with columns delta, strike and impvol
#' @export
blackcs_deltas <- function(type, future, time, rate, sigma, mu3, mu4, deltas = DELTAS) {
  eps <- 1e-8

  strk_cs <- sapply(deltas / 100, function(delta) {
    uniroot(
      f = blackcs_strikes, interval = c(1, future * 10), tol = eps,
      delta = delta, type = type, future = future, time = time,
      rate = rate, sigma = sigma, mu3 = mu3, mu4 = mu4
    )$root
  })

  impvol_cs <- sapply(strk_cs, function(strike) {
    blackimpvol(
      blackcsprice(type, future, strike, time, rate, sigma, mu3, mu4),
      type, future, strike, time, rate
    )
  })


  data.frame(delta = deltas, strike = strk_cs, impvol = impvol_cs)
}

blackcs_strikes <- function(x, delta, type, future, time, rate, sigma, mu3, mu4) {
  csp <- blackcsprice(type, future, x, time, rate, sigma, mu3, mu4)
  if (csp < 0) {
    return(999)
  }
  varx <- tryCatch(
    {
      blackimpvol(csp, type, future, x, time, rate)
    },
    error = function(e) {
      1e-8
    }
  )
  x - blackstrike(future, delta, varx, time)
}



#' Corrado-Su std. Deltas
#'
#' Finds and returns a dataframe with the strikes and implied
#' volatilites for the provided set of deltas or for the
#' default deltas.
#'
#' @param type Option type, "call" for calls or anything else for puts
#' @param spot underlying spot price at a given date
#' @param time time in years from a given date to a maturity date
#' @param rate interest free rate
#' @param yield dividend yield
#' @param sigma volatility of the stock price
#' @param mu3 skewness
#' @param mu4 kurtosis
#' @param deltas deltas for which strikes and vols should be calculated, 1-99. Defaults to usual
#' delta values.
#'
#' @return
#' Dataframe with columns delta, strike and impvol
#' @export
cs_deltas <- function(type, spot, time, rate, yield, sigma, mu3, mu4, deltas = DELTAS) {
  ft <- spot * exp((rate - yield) * time)
  eps <- 1e-8

  strk_cs <- sapply(deltas / 100, function(delta) {
    uniroot(
      f = cs_strikes, interval = c(1, spot * 10), tol = eps,
      delta = delta, type = type, spot = spot, time = time,
      rate = rate, yield = yield, sigma = sigma, mu3 = mu3, mu4 = mu4, ft = ft
    )$root
  })

  impvol_cs <- blackimpvol(
    csprice(type, spot, strk_cs, time, rate, yield, sigma, mu3, mu4),
    type, ft, strk_cs, time, rate
  )


  data.frame(delta = deltas, strike = strk_cs, impvol = impvol_cs)
}

cs_strikes <- function(x, delta, type, spot, strike, time, rate, yield, sigma, mu3, mu4, ft) {
  csp <- csprice(type, spot, x, time, rate, yield, sigma, mu3, mu4)
  if (csp < 0) {
    return(999)
  }
  varx <- blackimpvol(csp, type, ft, x, time, rate)
  x - blackstrike(ft, delta, varx, time)
}
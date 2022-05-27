#' Computes european call option price via Heston Model
#'
#' Heston formula implementation, pricing european call options.
#'
#' @param S spot price
#' @param k the strike price.
#' @param tau annualized period between option reference date and maturity date
#' @param r interest free rate
#' @param q dividend yield rate
#' @param lambda speed of the reversion mean of the volatility
#' @param eta volatility of the volatility
#' @param rho correlation between the Brownian motion of St and vt.
#' @param vbar long term mean of volatility
#' @param v0 initial variance
#'
#' @return
#' Vector with price of european options, according to
#' Heston model.
#'
#' @export
heston_price <- function(S, K, tau, r, q, lambda, eta, rho, vbar, v0) {
  f <- S * exp((r - q) * tau)
  x <- log(f / K)
  integral <- heston_phi_transform(tau, x, lambda, eta, rho, vbar, v0)

  f - (sqrt(K * f) / (2 * pi)) * integral
}

#' Computes european call option price via Heston Model
#'
#' Heston formula implementation, pricing european call futures options.
#'
#' @param F future price
#' @param k the strike price.
#' @param tau annualized period between option reference date and maturity date
#' @param lambda speed of the reversion mean of the volatility
#' @param eta volatility of the volatility
#' @param rho correlation between the Brownian motion of St and vt.
#' @param vbar long term mean of volatility
#' @param v0 initial variance
#'
#' @return
#' Vector with price of european options, according to
#' Heston model.
#'
#' @export
heston_pricef <- function(f, K, tau, lambda, eta, rho, vbar, v0) {
  x <- log(f / K)

  integral <- heston_phi_transform(tau, x, lambda, eta, rho, vbar, v0)

  f - (sqrt(K * f) / (2 * pi)) * integral
}

#' Heston model fit
#'
#' Fits a heston model to the option data provided and returns
#' the model parameters.
#'
#' @param forwards forward prices of the underlying
#' @param strikes strike prices of the options
#' @param prices options closing prices for the reference dates
#' @param terms annualized periods between option reference dates and maturity dates
#' @param initial_guess vector with the initial guess for the parameters lambda, eta, rho, vbar and v0
#' @param uncertainties price values uncertainties in \% for weighting
#' @param do_imp_vol if TRUE will minimize implied vol error (market imp vols must be passed in prices param)
#' @param do_global if TRUE does global optimization to find starting point for local optimization
#'
#' @return
#' List with heston model parameters lambda, eta, rho, vbar and v0
#'
#' @export
heston_fitf <- function(forwards, strikes, prices, terms,
                        initial_guess = NULL, uncertainties = NULL,
                        do_imp_vol = FALSE, do_global = TRUE) {
  min_function <- heston_minsqf(forwards, strikes, prices, terms, uncertainties, do_imp_vol)
  heston_optimize(min_function, initial_guess, do_global)
}

heston_minsqf <- function(forwards, strikes, prices, terms, uncertainties = NULL, do_imp_vol = FALSE) {
  function(...) {
    x <- c(...)
    lambda <- x[1]
    eta <- x[2]
    rho <- x[3]
    vbar <- x[4]
    v0 <- x[5]
    heston_prices <- numeric()
    # Integration doesn't work well for vectors
    for (i in 1:length(forwards)) {
      heston_prices[i] <- heston_pricef(
        forwards[i], strikes[i], terms[i],
        lambda, eta, rho, vbar, v0
      )
    }

    if (do_imp_vol) {
      heston_prices <- tryCatch(
        {
          gmoa.pricing::bsmimpvol(heston_prices, "call", spots, strikes, terms, rates, yields)
        },
        error = function(e) {
          0
        }
      )
    }

    if (is.null(uncertainties)) {
      sum(((heston_prices - prices) / prices)^2)
    } else {
      sum(((((heston_prices - prices) / prices)) * 1 / ((uncertainties / 100)))^2)
    }
  }
}



#' Heston model fit
#'
#' Fits a heston model to the stock option data provided and returns
#' the model parameters.
#'
#' @param spots forward prices of the underlying
#' @param strikes strike prices of the options
#' @param prices options closing prices for the reference dates
#' @param rates interest free rates
#' @param yields dividend yields
#' @param terms annualized periods between option reference dates and maturity dates
#' @param initial_guess vector with the initial guess for the parameters lambda, eta, rho, vbar and v0
#' @param uncertainties price values uncertainties in \% for weighting
#' @param do_imp_vol if TRUE will minimize implied vol error (market imp vols must be passed in prices param)
#' @param do_global if TRUE does global optimization to find starting point for local optimization
#'
#' @return
#' List with heston model parameters lambda, eta, rho, vbar and v0
#'
#' @export
heston_fit <- function(spots, strikes, prices, rates, yields, terms,
                       initial_guess = NULL, uncertainties = NULL,
                       do_imp_vol = FALSE, do_global = TRUE) {
  min_function <- heston_minsq(spots, strikes, prices, rates, yields, terms, uncertainties, do_imp_vol)
  heston_optimize(min_function, initial_guess, do_global)
}

heston_minsq <- function(spots, strikes, prices, rates, yields, terms, uncertainties = NULL, do_imp_vol = FALSE) {
  function(...) {
    x <- c(...)
    lambda <- x[1]
    eta <- x[2]
    rho <- x[3]
    vbar <- x[4]
    v0 <- x[5]
    heston_prices <- numeric()
    # Integration doesn't work well for vectors
    for (i in 1:length(spots)) {
      heston_prices[i] <- heston_price(
        spots[i], strikes[i], terms[i], rates[i],
        yields[i], lambda, eta, rho, vbar, v0
      )
    }

    if (do_imp_vol) {
      heston_prices <- tryCatch(
        {
          gmoa.pricing::bsmimpvol(heston_prices, "call", spots, strikes, terms, rates, yields)
        },
        error = function(e) {
          0
        }
      )
    }

    if (is.null(uncertainties)) {
      sum(((heston_prices - prices) / prices)^2)
    } else {
      sum(((((heston_prices - prices) / prices)) * 1 / ((uncertainties / 100)))^2)
    }
  }
}

heston_phi <- function(k, tau, lambda, eta, rho, vbar, v0) {
  b <- lambda + complex(imaginary = rho * eta) * k
  d <- sqrt(b^2 + (eta^2) * (k) * (k + complex(imaginary = -1)))
  g <- (b - d) / (b + d)
  T_m <- (b - d) / (eta^2)
  t <- T_m * (1 - exp(-d * tau)) / (1 - g * exp(-d * tau))
  W <- lambda * vbar * (tau * T_m - 2 * log((1 - g * exp(-d * tau)) / (1 - g)) / (eta^2))

  exp(W + v0 * t)
}

heston_phi_transform <- function(tau, x, lambda, eta, rho, vbar, v0) {
  integrand <- function(k) {
    result <- 2 * Re(exp(complex(imaginary = -k * x)) * heston_phi(
      complex(real = k, imaginary = 0.5), tau,
      lambda, eta, rho, vbar, v0
    )) / (k^2 + 1 / 4)
    if (!all(is.numeric(result))) browser()
    if (any(is.nan(result))) browser()
    if (any(is.na(result))) browser()
    result
  }
  integrate(integrand, 0, 100, stop.on.error = F)$value
}

heston_optimize <- function(objective_function, initial_guess = NULL, do_global = FALSE) {
  if (is.null(initial_guess)) {
    x0 <- c(1e-6, 1e-6, -1, 1e-6, 1e-6)
  } else {
    x0 <- initial_guess
  }
  lo <- c(1e-12, 1e-12, -1, 1e-12, 1e-12)
  hi <- c(20, 5, 0, 1, 1)

  # Feller's condition
  # 2*lambda*vbar - eta^2 > 0
  # lambda, eta, rho, vbar, v0
  eval_g0 <- function(x) {
    c(x[2] * x[2] - (2 * x[1] * x[4]))
  }

  if (do_global) {
    gl_fit <- nloptr::nloptr(
      x0 = x0,
      eval_f = objective_function,
      lb = lo,
      ub = hi,
      eval_g_ineq = eval_g0,
      opts = list(
        algorithm = "NLOPT_GN_ISRES",
        print_level = 0,
        xtol_rel = 1.0e-6, maxeval = 500
      )
    )

    x0 <- gl_fit$solution
  }

  fit <- nloptr::nloptr(
    x0 = x0,
    eval_f = objective_function,
    lb = lo,
    ub = hi,
    eval_g_ineq = eval_g0,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      print_level = 0,
      xtol_rel = 1.0e-6, maxeval = 500
    )
  )


  params <- fit$solution

  list(
    lambda = params[1], eta = params[2], rho = params[3],
    vbar = params[4], v0 = params[5], objective = fit$objective
  )
}


#' Heston model std. Deltas
#'
#' Finds and returns a dataframe with the strikes and implied
#' volatilites for the provided set of deltas or for the
#' default deltas.
#'
#' @param fit a Heston model fit as returned by heston_fit
#' @param S underlying spot price at a given date
#' @param tau time in years from a given date to a maturity date
#' @param r interest free rate
#' @param q dividend yield
#' @param future underlying future price, can be provided instead of spot and yield
#' @param deltas deltas for which strikes and vols should be calculated. Defaults to usual
#' delta values.
#'
#' @return
#' Dataframe with columns delta, strike and impvol
#'
#' @export
heston_deltas <- function(fit, S = NULL, tau, r, q = NULL, future = NULL, deltas = DELTAS) {
  if (is.null(future)) {
    ft <- S * exp((r - q) * tau)
  } else {
    ft <- future
  }

  if (is.null(ft)) stop("Future price or spot, rate and yield must be provided.")

  eps <- 1e-8

  # heston price não funciona vetorialmente devido à integral,
  # então as chamadas precisam ficar em apply

  strk_heston <- sapply(deltas / 100, function(delta) {
    uniroot(
      f = heston_strikes, interval = c(1, ft * 10), tol = eps,
      t = tau, delta = delta, ft = ft, fit = fit, r = r
    )$root
  })

  impvol_heston <- sapply(strk_heston, function(strike) {
    blackimpvol(heston_pricef(ft, strike, tau, fit$lambda, fit$eta, fit$rho, fit$vbar, fit$v0), "call", ft, strike, tau, r)
  })

  data.frame(delta = deltas, strike = strk_heston, impvol = impvol_heston)
}

heston_strikes <- function(x, delta, ft, t, fit, r) {
  hp <- heston_pricef(ft, x, t, fit$lambda, fit$eta, fit$rho, fit$vbar, fit$v0)
  if (hp < 0) {
    return(999)
  }
  varx <- tryCatch(
    {
      blackimpvol(hp, "call", ft, x, t, r)
    },
    error = function(e) {
      1e-8
    }
  )
  x - blackstrike(ft, delta, varx, t)
}
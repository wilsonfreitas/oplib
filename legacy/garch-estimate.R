var.garch <- function(rets, omega, alpha, beta) {
  mu <- mean(rets)
  e2 <- (rets - mu)^2
  e2t <- omega + alpha * e2[c(-1, -length(e2))]
  s2 <- stats::filter(e2t, beta, "recursive", init = mean(e2))
  c(mean(e2), s2)
}

# source for student and ged formulas:
# https://hal.inria.fr/file/index/docid/270719/filename/B08027.pdf
.garch_likelihood <- function(rets, garchEnv, dist, scale = c(1, 1, 1)) {
  .rets <- rets
  .garchEnv <- garchEnv
  llh_function <- switch(dist,
    normal = {
      function(...) {
        x <- c(...)
        omega <- x[1] / scale[1] # 1e-5
        alpha <- x[2] / scale[2] # 1e-1
        beta <- x[3] / scale[3]
        garch <- var.garch(.rets, omega, alpha, beta)
        # likelihood = 0.5 sum(-log(2 pi) - log(garch) - rets^2/garch)
        # the implementation below is a simplification
        if (any(garch < 0)) {
          .llh <- get("llh", env = .garchEnv)
          llh <- .llh + 0.1 * (abs(.llh))
        } else {
          llh <- 0.5 * sum(log(garch) + .rets[-1]^2 / garch + log(2 * pi))
        }
        assign("llh", llh, env = .garchEnv)
        llh
      }
    },
    student = {
      function(...) {
        x <- c(...)
        omega <- x[1] / scale[1] # 1e-5
        alpha <- x[2] / scale[2] # 1e-1
        beta <- x[3] / scale[3]
        v <- x[4]
        n <- length(.rets)
        garch <- var.garch(.rets, omega, alpha, beta)
        s1 <- n * (lgamma((v + 1) / 2) - lgamma(v / 2) -
          (1 / 2) * log(pi * (v - 2)))
        ll <- 1 + ((.rets[-1])^2 / (garch * (v - 2)))
        if (is.na(v) || any(garch < 0) || any(ll < 0)) {
          .llh <- get("llh", env = .garchEnv)
          llh <- .llh + 0.1 * (abs(.llh))
        } else {
          llh <- -(s1 - 0.5 * (sum(log(garch) + (v + 1) * log(ll))))
        }
        assign("llh", llh, env = .garchEnv)
        llh
      }
    },
    ged = {
      function(...) {
        x <- c(...)
        omega <- x[1] / scale[1] # 1e-5
        alpha <- x[2] / scale[2] # 1e-1
        beta <- x[3] / scale[3]
        v <- x[4]
        n <- length(.rets)
        yv <- sqrt((2^(-2 / v)) * gamma(1 / v) / gamma(3 / v))
        garch <- var.garch(.rets, omega, alpha, beta)
        s1 <- n * (log(v / yv) - (1 + 1 / v) * log(2) - lgamma(1 / v))
        if (is.na(v) || any(garch < 0)) {
          .llh <- get("llh", env = .garchEnv)
          llh <- .llh + 0.1 * (abs(.llh))
        } else {
          llh <- -(s1 - 0.5 * (sum(log(garch) + ((sqrt(garch))^(-v)) * ((abs(.rets[-1] / yv))^(v)))))
        }
        assign("llh", llh, env = .garchEnv)
        llh
      }
    }
  )
  llh_function
}

garch_fit_1 <- function(rets, garchEnv, dist, bind_params) {
  tiny <- 1e-3
  alpha <- 0.1
  beta <- 0.8
  omega <- (1 - alpha - beta) * var(rets)
  mu <- mean(rets)
  v <- 1

  x0 <- c(omega, alpha, beta)

  # scaling and optimization bounds
  small <- 1e-12

  lo <- c(small, small, small)
  hi <- c(2 * var(rets), 1 - small, 1 - small)

  # setting up likelihood function
  garch_likelihood <- .garch_likelihood(rets, garchEnv, dist)

  # optimization

  eval_g0 <- function(x) {
    c(-x[1], -x[2], -x[3], sum(x[2:3]) - (1 - tiny))
  }


  # student needs the degrees of freedom to be optimized too
  if (dist == "student") {
    v <- 3
    x0 <- c(x0, v)
    lo <- c(lo, 3)
    hi <- c(hi, 10)
    eval_g0 <- function(x) {
      c(-x[1], -x[2], -x[3], sum(x[2:3]) - (1 - tiny), -x[4])
    }
  } else if (dist == "ged") { # ged needs nu to be estimated
    x0 <- c(x0, v)
    lo <- c(lo, 1)
    hi <- c(hi, 10)
    eval_g0 <- function(x) {
      c(-x[1], -x[2], -x[3], sum(x[2:3]) - (1 - tiny), -x[4])
    }
  }

  # optimization
  if (bind_params) {
    res <- nloptr::nloptr(
      x0 = x0,
      eval_f = garch_likelihood,
      lb = lo,
      ub = hi,
      eval_g_ineq = eval_g0,
      opts = list(
        algorithm = "NLOPT_LN_COBYLA",
        print_level = 0,
        xtol_rel = 1.0e-12, maxeval = 2000
      )
    )
  } else {
    res <- nloptr::nloptr(
      x0 = x0,
      eval_f = garch_likelihood,
      lb = lo,
      ub = hi,
      opts = list(
        algorithm = "NLOPT_LN_COBYLA",
        print_level = 0,
        xtol_rel = 1.0e-12, maxeval = 2000
      )
    )
  }

  par <- res$solution
  parms <- list(omega = par[1], alpha = par[2], beta = par[3])
  var.g <- var.garch(rets, parms$omega, parms$alpha, parms$beta)
  current_var <- sum(unlist(parms) * c(1, tail(rets, 1)^2, tail(var.g, 1)))
  fit <- as.garchparms(c(parms, sigma_t = annualize(as.volatility(as.variance(current_var)))))
  fit$shape <- par[4]
  fit
}

.garch_likelihood_2 <- function(rets, garchEnv) {
  .rets <- rets
  .var <- as.numeric(var(rets))
  .garchEnv <- garchEnv

  function(...) {
    x <- c(...)
    omega <- .var * (1 - sum(x))
    alpha <- x[1]
    beta <- x[2]
    garch <- var.garch(.rets, omega, alpha, beta)
    # likelihood = 0.5 sum(-log(2 pi) - log(garch) - rets^2/garch)
    # the implementation below is a simplification
    if (any(garch < 0)) {
      .llh <- get("llh", env = .garchEnv)
      llh <- .llh + 0.1 * (abs(.llh))
    } else {
      llh <- 0.5 * sum(log(garch) + .rets[-1]^2 / garch + log(2 * pi))
    }
    assign("llh", llh, env = .garchEnv)
    llh
  }
}

garch_fit_2 <- function(rets, garchEnv, bind_params) {
  tiny <- 1e-3
  alpha <- 0.2
  beta <- 0.7
  mu <- mean(rets)

  x0 <- c(alpha, beta)

  # scaling and optimization bounds
  small <- 1e-12

  lo <- c(small, small)
  hi <- c(1 - small, 1 - small)

  # setting up likelihood function
  garch_likelihood <- .garch_likelihood_2(rets, garchEnv)

  # optimization

  eval_g0 <- function(x) {
    c(-x[1], -x[2], sum(x[1:2]) - (1 - tiny))
  }

  if (bind_params) {
    res <- nloptr::nloptr(
      x0 = x0,
      eval_f = garch_likelihood,
      lb = lo,
      ub = hi,
      eval_g_ineq = eval_g0,
      opts = list(
        algorithm = "NLOPT_LN_COBYLA",
        print_level = 0,
        xtol_rel = 1.0e-12, maxeval = 2000
      )
    )
  } else {
    res <- nloptr::nloptr(
      x0 = x0,
      eval_f = garch_likelihood,
      lb = lo,
      ub = hi,
      opts = list(
        algorithm = "NLOPT_LN_COBYLA",
        print_level = 0,
        xtol_rel = 1.0e-12, maxeval = 2000
      )
    )
  }

  par <- res$solution
  parms <- list(omega = as.numeric(var(rets)) * (1 - sum(par)), alpha = par[1], beta = par[2])
  var.g <- var.garch(rets, parms$omega, parms$alpha, parms$beta)
  current_var <- sum(unlist(parms) * c(1, tail(rets, 1)^2, tail(var.g, 1)))
  fit <- as.garchparms(c(parms, sigma_t = annualize(as.volatility(as.variance(current_var)))))
  fit$shape <- NA
  fit
}


#' GARCH time series fitting
#'
#' Estimate parameters for a GARCH model, using Quasi-Maximum Likelihood Estimation
#' for optimization.
#'
#' @param rets log returns to be fitted
#' @param dist the conditional distribution to be used for estimation
#' @param bind_params if TRUE the fit will be bound by the alpha + beta < 1 restriction
#' @param force_asgard if TRUE will use asgard model instead of fair-garch
#'
#' @return
#' \code{"fit"} a list with the estimated omega, alpha and beta parameters,
#' and sigma, LTV, model, distribution and covariance matrix.
#'
#' @export
garch_fit <- function(rets,
                      dist = c("normal", "student", "ged"),
                      bind_params = TRUE,
                      force_asgard = FALSE) {
  dist <- match.arg(dist)
  garchEnv <- new.env(hash = TRUE)
  # Keep the calculated likelihood function values to handle NaN
  assign("llh", 1.0e99, env = garchEnv)
  fit <- garch_fit_1(rets, garchEnv, dist, bind_params)
  LTV <- ltv(fit)
  model <- "fair-garch"
  if ((LTV < 0 || as.numeric(annualize(as.volatility(LTV))) > 500) || force_asgard) {
    fit <- garch_fit_2(rets, garchEnv, bind_params)
    model <- "asgard"
  }
  LTV <- ltv(fit)
  LTV_annual <- as.numeric(annualize(as.volatility(LTV)))
  LTV <- as.numeric(LTV)

  SIGMA_T <- as.numeric(fit$sigma_t)

  fit$cvar <- cvar(rets, dist, fit$omega, fit$alpha, fit$beta, fit$shape)
  ci <- conf_interval(fit$omega, fit$alpha, fit$beta, fit$cvar, LTV, LTV_annual)


  list(
    omega = fit$omega, alpha = fit$alpha, beta = fit$beta, shape = fit$shape,
    sigma_t = SIGMA_T, LTV = LTV, LTV_annual = LTV_annual,
    model = model, dist = dist, cvar = fit$cvar, ci = ci
  )
}


cvar <- function(rets, dist, omega, alpha, beta, shape) {
  garchEnv <- new.env(hash = TRUE)
  assign("llh", 1.0e99, env = garchEnv) # Keep the calculated likelihood function values to handle NaN

  par <- c(omega, alpha, beta)
  if (!is.na(shape)) {
    par <- c(par, shape)
  }

  H <- hessian(.garch_likelihood(rets, garchEnv, dist), par,
    method.args = list(eps = 1e-12)
  )
  H <- 0.5 * (H + t(H))
  solve(H, tol = 1e-17)
}

conf_interval <- function(omega, alpha, beta, cvar, LTV, LTV_annual) {
  dv_dw <- 1 / (1 - alpha - beta)
  dv_da <- omega / ((1 - alpha - beta))^2
  dv_db <- dv_da
  var_vl <- (dv_dw^2) * cvar[1, 1] + (dv_da^2) * cvar[2, 2] +
    (dv_db^2) * cvar[3, 3] + 2 * dv_dw * dv_da * cvar[1, 2] +
    2 * dv_dw * dv_db * cvar[1, 3] + 2 * dv_da * dv_db * cvar[2, 3]

  .se_lt_var <- suppressWarnings(sqrt(var_vl))
  volp <- as.volatility(LTV + .se_lt_var)[[1]]
  volm <- as.volatility(LTV - .se_lt_var)[[1]]
  if (!is.nan(volm) && !is.null(volm) && !is.na(volm) && volm < 0) volm <- 0
  icmax <- sqrt(volp * 252) * 100
  icmin <- sqrt(volm * 252) * 100
  interv <- (icmax - icmin) / 2
  lo <- LTV_annual - interv
  hi <- LTV_annual + interv

  data.frame(
    Annualized_Volatility = LTV_annual,
    CI_low = lo,
    CI_high = hi,
    In_CI = (LTV_annual <= hi && LTV_annual >= lo),
    se_lt_var = .se_lt_var
  )
}
#' garchparms Class Creator
#'
#' Create list of class \code{"garchparms"}.
#'
#' @param omega Generalized AutoRegressive Conditional Heteroskedasticity
#' process paramater.
#' @param alpha GARCH process paramater.
#' @param beta GARCH process paramater.
#' @param sigma_t GARCH process paramater.
#'
#' @return
#' List with class \code{"garchparms"} from given arguments.
#'
#' @seealso
#' \code{\link{as.garchparms}}.
#'
#' @examples
#' garchparms(1, 1, 1, 1)
#' class(garchparms(1, 1, 1, 1))
#' typeof(garchparms(1, 1, 1, 1))
#'
#' @export
garchparms <- function(omega, alpha, beta, sigma_t) {
  as.garchparms(list(
    omega = omega,
    alpha = alpha,
    beta = beta,
    sigma_t = sigma_t
  ))
}

#' garchparms Class Coercing
#'
#' Generic function for coercing to \code{"garchparms"}. There is no default
#' method, so only classes \code{"data.frame"}, \code{"numeric"}, \code{"list"}
#' are supported.
#'
#' @param x an object whose class will determine the method to be dispatched. Must
#' have the values for omega, alpha, beta, sigma_t.
#' @param ... further arguments to be passed to the next method.
#'
#' @return
#' List with class \code{garchparms} from given argument. When x has
#' \code{"data.frame"} class, idx = FALSE or equivalent will pass numeric(0) values
#' and numbers in idx will produce NA.
#'
#' @examples
#' class(garchparms(1, 1, 1, 1)) == class(list(omega = 1, alpha = 1, beta = 1, sigma_t = 1))
#' x <- data.frame(omega = 1, alpha = 1, beta = 1, sigma_t = 1)
#' as.garchparms(x)
#' as.garchparms(x, 0)
#' as.garchparms(x, 2.579238)
#'
#' @name as.gachparms
#' @export
as.garchparms <- function(x, ...) {
  UseMethod("as.garchparms", x)
}

#' @rdname as.gachparms
#'
#' @param idx additional parameter.
#'
#' @export
as.garchparms.data.frame <- function(x, idx = 1, ...) {
  as.garchparms.list(unclass(x[idx, ]))
}

#' @rdname as.gachparms
#' @export
as.garchparms.numeric <- function(x, ...) {
  as.garchparms.list(as.list(x))
}

#' @rdname as.gachparms
#' @export
as.garchparms.list <- function(x, ...) {
  if (is.null(names(x))) {
    stop("names not defined")
  }
  aux <- c("omega", "alpha", "beta", "sigma_t") %in% tolower(names(x))
  if (!all(aux)) {
    stop("names mismatch")
  }
  x$sigma_t <- as.volatility(x$sigma_t, annualized = T, dib = 252)
  structure(x, class = "garchparms")
}

#' GARCH Information
#'
#' Gives relevant information about a GARCH model, such as parameters and
#' volatility. There is no default method, so only classes \code{"garch"},
#' \code{"fGARCH"}, \code{"uGARCHfit"} are supported.
#'
#' Methods for printing classes \code{"info.garch"}, \code{"info.fGARCH"},
#' \code{"info.uGARCHfit"} are also included.
#'
#' @param x a GARCH model.
#' @param ... further arguments to be passed to the next method.
#' @param dib days in one year, day count convention.
#'
#' @return
#' List with parameters, daily long term variance, annual long term variance, fit
#' of class "info.class(x)". Printing info classes will show daily and
#' annualized variance and volatility, plus parameters.
#'
#' @name info
#' @export
info <- function(x, ...) UseMethod("info")

#' @rdname info
#' @export
info.garch <- function(x, ..., dib = 252) {
  gamma <- ggamma(x)
  lt_var <- coef(x)["omega"] / gamma
  lt_vol <- sqrt(lt_var)
  lt_vol_year <- lt_vol * sqrt(dib)
  .vcov <- sqrt(diag(vcov(x)))
  names(.vcov) <- c("omega", "alpha", "beta")
  parms <- c(coef(x)["a0"], gamma, coef(x)["a1"], coef(x)["b1"])
  names(parms) <- c("omega", "gamma", "alpha", "beta")
  .se_lt_var <- sqrt(sum(diag(vcov(x)) * c(1, lt_var^2, lt_var^2)) / (gamma^2))
  longterm <- c(lt_var, lt_vol, .se_lt_var)
  names(longterm) <- c("variance", "volatility", "se")
  longterm.annual <- c(lt_vol_year^2, lt_vol_year)
  names(longterm.annual) <- c("variance", "volatility")
  ans <- list(
    garch.parms = parms, longterm = longterm,
    longterm.annual = longterm.annual,
    fit = x, se = .vcov
  )
  class(ans) <- "info.garch"
  ans
}

#' @rdname info
#' @export
print.info.garch <- function(x, ...) {
  a <- x
  class(a) <- NULL
  cat("\nLong Term Info (%):", "\n")
  print.default(a$longterm * rep(100, length(a$longterm)),
    digits = 3,
    print.gap = 4
  )
  cat("\nAnnualized Volatility (%):", "\n")
  print.default(a$longterm.annual * rep(1, length(a$longterm.annual)),
    digits = 3, print.gap = 4
  )
  cat("\nGARCH parameters:", "\n")
  print.default(a$garch.parms, digits = 3, print.gap = 4)
  invisible(x)
}

#' @rdname info
#' @export
info.fGARCH <- function(x, ..., dib = 252) {
  # long term info
  gamma <- ggamma(x)
  lt_var <- ltv(x)
  lt_vol <- as.volatility(lt_var)
  lt_vol_year <- annualize(lt_vol, dib)
  # se info
  .vcov <- sqrt(diag(x@fit$cvar))
  names(.vcov) <- c("omega", "alpha", "beta")
  # parameters
  parms <- c(coef(x)["omega"], gamma, coef(x)["alpha1"], coef(x)["beta1"])
  names(parms) <- c("omega", "gamma", "alpha", "beta")
  # long term variance se
  .se_lt_var <- sqrt(sum(diag(x@fit$cvar) * c(1, lt_var^2, lt_var^2)) / (gamma^2))
  longterm <- c(lt_var, lt_vol, .se_lt_var)
  names(longterm) <- c("variance", "volatility", "se")
  # annualized info
  longterm.annual <- c(lt_vol_year^2, lt_vol_year)
  names(longterm.annual) <- c("variance", "volatility")
  # info structure
  ans <- list(
    garch.parms = parms, longterm = longterm,
    longterm.annual = longterm.annual,
    fit = x, se = .vcov
  )
  class(ans) <- "info.garch"
  ans
}

#' @rdname info
#' @export
print.info.fGARCH <- function(x, ...) {
  print.info.garch(x, ...)
}

#' @rdname info
#' @export
info.uGARCHfit <- function(x, ..., dib = 252) {
  # long term info
  gamma <- ggamma(x)
  lt_var <- ltv(x)
  lt_vol <- as.volatility(lt_var)
  lt_vol_year <- annualize(lt_vol, dib)
  # se info
  .vcov <- sqrt(diag(x@fit$cvar))
  names(.vcov) <- c("omega", "alpha", "beta")
  # parameters
  parms <- c(coef(x)["omega"], gamma, coef(x)["alpha1"], coef(x)["beta1"])
  names(parms) <- c("omega", "gamma", "alpha", "beta")
  # long term variance se
  .se_lt_var <- sqrt(sum(diag(x@fit$cvar) * c(1, lt_var^2, lt_var^2)) / (gamma^2))
  longterm <- c(lt_var, lt_vol, .se_lt_var)
  names(longterm) <- c("variance", "volatility", "se")
  # annualized info
  longterm.annual <- c(lt_vol_year^2, lt_vol_year)
  names(longterm.annual) <- c("variance", "volatility")
  # info structure
  ans <- list(
    garch.parms = parms, longterm = longterm,
    longterm.annual = longterm.annual,
    fit = x, se = .vcov
  )
  class(ans) <- "info.garch"
  ans
}

#' @rdname info
#' @export
print.info.uGARCHfit <- function(x, ...) {
  print.info.garch(x, ...)
}

#' GARCH Model Long Term Variance
#'
#' Calculates long term variance for GARCH models. There is no default
#' method, so only classes \code{"garch"}, \code{"fGARCH"}, \code{"uGARCHfit"},
#' \code{"garchparms"} are supported.
#'
#' @param x object from where omega, alpha and beta parameters will be extracted.
#'
#' @return
#' A structure with daily long term variance and class \code{"variance"}.
#'
#' @seealso
#' .
#'
#' @examples
#' ltv(garchparms(1, 2:4 / 20, .7, 1))
#'
#' @name ltv
#'
#' @export
ltv <- function(x) UseMethod("ltv")

#' @rdname ltv
#' @export
ltv.fGARCH <- function(x) as.variance(coef(x)["omega"] / ggamma(x))

#' @rdname ltv
#' @export
ltv.garch <- function(x) as.variance(coef(x)["a0"] / ggamma(x))

#' @rdname ltv
#' @export
ltv.garchparms <- function(x) as.variance(x$omega / ggamma(x))

#' @rdname ltv
#' @export
ltv.uGARCHfit <- function(x) as.variance(coef(x)["omega"] / ggamma(x))


#' GARCH Model Long Term Variance Decay
#'
#' Calculates long term variance decay for given alpha and beta. There is
#' no default method, so only classes \code{"fGARCH"}, \code{"garch"},
#' \code{"garchparms"}, \code{"uGARCHfit"} are supported.
#'
#' @param x object from where alpha and beta parameters will be extracted.
#'
#' @return
#' Numeric with long term variance decay for given parameters.
#'
#' @examples
#' decay(garchparms(1, 0:10 / 10, 1, 1))
#'
#' @name decay
#' @export
decay <- function(x) UseMethod("decay")

#' @rdname decay
#' @export
decay.fGARCH <- function(x) {
  a <- log(1 / (coef(x)["alpha1"] + coef(x)["beta1"]))
  as.numeric(a)
}

#' @rdname decay
#' @export
decay.garch <- function(x) {
  a <- log(1 / (coef(x)["a1"] + coef(x)["b1"]))
  as.numeric(a)
}

#' @rdname decay
#' @export
decay.garchparms <- function(x) {
  a <- log(1 / (x$alpha + x$beta))
  as.numeric(a)
}

#' @rdname decay
#' @export
decay.uGARCHfit <- function(x) {
  a <- log(1 / (coef(x)["alpha1"] + coef(x)["beta1"]))
  as.numeric(a)
}

#' GARCH Model Gamma
#'
#' Computes gamma parameter for GARCH models. There is no default
#' method, so only classes \code{"garch"}, \code{"fGARCH"}, \code{"uGARCHfit"},
#' \code{"garchparms"} are supported.
#'
#' @param x an object whose class will determine the method to be dispatched.
#'
#' @return
#' Numeric with gamma calculated for given parameters.
#'
#' @seealso
#' .
#'
#' @examples
#' ggamma(garchparms(1, 2:4 / 20, .7, 1))
#'
#' @name ggamma
#' @export
ggamma <- function(x) UseMethod("ggamma")

#' @rdname ggamma
#' @export
ggamma.garch <- function(x) 1 - coef(x)["a1"] - coef(x)["b1"]

#' @rdname ggamma
#' @export
ggamma.fGARCH <- function(x) 1 - coef(x)["alpha1"] - coef(x)["beta1"]

#' @rdname ggamma
#' @export
ggamma.garchparms <- function(x) 1 - x$alpha - x$beta

#' @rdname ggamma
#' @export
ggamma.uGARCHfit <- function(x) 1 - coef(x)["alpha1"] - coef(x)["beta1"]

#' GARCh Model Variance Term Structure
#'
#' Computes variance term structure for GARCH models. There is a default method,
#' where parameters for formula may be passed directly, and methods for classes
#' \code{"fGARCH"}, \code{"garchparms"}, \code{"uGARCHfit"}, where parameters are
#' calcuated from the GARCH models.
#'
#' @param x a GARCH model, GARCH paramters or ignored for default method.
#' @param ... further arguments to be passed to the next method.
#' @param v0 for default method, sigma_t converted to daily variance.
#' @param a for default method, long term variance decay.
#' @param vl for default method, long term variance.
#'
#' @return
#' Term structure of variance as a function of time.
#'
#' @seealso
#' \code{\link{decay}}, \code{\link{ltv}}.
#'
#' @examples
#' vts(garchparms(1, 1, 1, 1))(1)
#' vts("ignored", .5, .5, .3)(1)
#'
#' @name vts
#' @export
vts <- function(x, ...) UseMethod("vts")

#' @rdname vts
#' @export
vts.default <- function(x, v0, a, vl, ...) {
  function(t) {
    q <- (1 - exp(-a * t)) / (a * t)
    as.variance(vl + q * (v0 - vl))
  }
}

#' @rdname vts
#' @export
vts.fGARCH <- function(x, v0 = NULL, ...) {
  v0 <- if (is.null(v0)) tail(fBasics::volatility(x, type = "h"), 1) else v0
  a <- decay(x)
  vl <- ltv(x)
  vts.default(x, v0 = v0, a = a, vl = vl)
}

#' @rdname vts
#' @export
vts.garchparms <- function(x, v0 = NULL, ...) {
  v0 <- if (is.null(v0)) daily(as.variance(x$sigma_t)) else v0
  a <- decay(x)
  vl <- ltv(x)
  vts.default(x, v0 = v0, a = a, vl = vl)
}

#' @rdname vts
#' @export
vts.uGARCHfit <- function(x, v0 = NULL, ...) {
  v0 <- if (is.null(v0)) daily(as.variance(x$sigma_t)) else v0
  a <- decay(x)
  vl <- ltv(x)
  vts.default(x, v0 = v0, a = a, vl = vl)
}

#' vts_update Volatility Term Structure Update
#'
#' Update volatility term structure given new return series data.
#'
#' @param garchdate current GARCH parameters estimation date.
#' @param refdate reference date for calculation.
#' @param returns return series to update vts.
#' @param omega omega GARCH parameter.
#' @param alpha alpha GARCH parameter.
#' @param beta beta GARCH parameter.
#' @param sigma_t sigma GARCH parameter.
#'
#' @return
#' Volatility Term Structure.
#'
#' @name vts_update
#' @export
vts_update <- function(garchdate, refdate, returns, omega, alpha, beta, sigma_t) {
  if (garchdate < refdate) {
    sigma_t <- sigma_t_update(returns, omega, alpha, beta, sigma_t)
    sigma_t <- as.numeric(sigma_t)
  }
  vts(garchparms(omega, alpha, beta, sigma_t))
}

# sigma(t) ----

#' sigma_t Calculation for GARCH Models
#'
#' Calculates sigma_t for GARCH models. There is no default method, so only
#' classes \code{"fGARCH"}, \code{"uGARCHfit"}, are supported.
#'
#' @param x a GARCH model with class \code{"fGARCH"} or \code{"uGARCHfit"}.
#' @param ... further arguments to be passed to the next method.
#'
#' @return
#' sigma_t for GARCH model in annualized volatility.
#'
#' @name sigma_t
#' @export
sigma_t <- function(x, ...) UseMethod("sigma_t")

#' @rdname sigma_t
#' @export
sigma_t.fGARCH <- function(x, ...) {
  last_ret <- tail(x@data, 1)
  last_var <- tail(x@sigma.t, 1)^2
  current_var <- sum(coef(x) * c(1, last_ret^2, last_var))
  annualize(as.volatility(as.variance(current_var)))
}

# #' @rdname sigma_t
# #' @export
# sigma_t.uGARCHfit <- function(fit) {
#   current_var <- rugarch::sigma(rugarch::ugarchforecast(fit, n.ahead = 1)) ^ 2
#   annualize(as.volatility(as.variance(current_var)))
# }

#' sigma_t_update Update GARCH paramater Sigma
#'
#' Given an additional return series, update sigma GARCH parameter, keeping the
#' others.
#'
#' @param returns return series to update sigma.
#' @param omega omega Garch parameter.
#' @param alpha alpha Garch parameter.
#' @param beta beta Garch parameter.
#' @param sigma_t sigma Garch parameter previously estimated.
#'
#' @return
#' Annualized sigma_t.
#'
#' @name sigma_t_update
#' @export
sigma_t_update <- function(returns, omega, alpha, beta, sigma_t) {
  sig <- daily(as.volatility(sigma_t, annualized = T, dib = 252))^2

  x <- omega + alpha * returns^2
  if (length(x) == 1) {
    sig <- x + beta * sig
  } else {
    sig <- filter(x, beta, method = "recursive", init = sig)
  }

  sig <- sig[length(sig)]
  annualize(as.volatility(as.variance(sig)))
}

# volatility ----

#' volatility_update Volatility Estimative Update
#'
#' Estimate new volatility value within initial date with known volatility
#' an desired final date.
#'
#' @param ticker underlying ticker identification.
#' @param refdate reference date for calculation.
#' @param matdate desired final date to witch update volatility.
#' @param garchdate current GARCH parameters estimation date.
#' @param omega omega GARCH parameter.
#' @param alpha alpha GARCH parameter.
#' @param beta beta GARCH parameter.
#' @param sigma_t sigma GARCH parameter.
#'
#' @return
#' Annualized volatility.
#'
#' @name volatility_update
#' @export
volatility_update <- function(ticker, refdate, matdate, garchdate, omega, alpha, beta, sigma_t) {
  vts <- vts_update(ticker, garchdate, refdate, omega, alpha, beta, sigma_t)
  bizdays <- bizdays::bizdays(refdate, matdate)

  volatility <- if (!is.null(vts)) vts(bizdays) else NA
  as.numeric(annualize(as.volatility(volatility))) / 100
}
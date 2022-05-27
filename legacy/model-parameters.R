#' Computes DI future option price via Black-Scholes formula.
#'
#' Black-Scholes formula implementation, pricing DI future options.
#'
#' @export
DELTAS <- c(1, 5, 10, 25, 35, 37, 50, 65, 67, 75, 90, 95, 99)

# Default SABR initial parameters
#' @export
ALPHA <- 0.001
#' @export
NU <- 0.001
#' @export
RHO <- 0.3
#' @export
BETA <- 0

# Default SVI initial parameters
#' @export
SVI.A <- 1e-12
#' @export
SVI.B <- 1e-8
#' @export
SVI.SIGMA <- 1e-8
#' @export
SVI.RHO <- -0.8
#' @export
SVI.M <- 1e-8
#' @export
SVI.ERR <- 0.02

# Default corrado-su initial parameters
#' @export
CS.SKEW <- 1e-8
#' @export
CS.KURT <- 1
#' @export
CS.SIGMA <- 1e-8
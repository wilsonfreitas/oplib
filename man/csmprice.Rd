% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csm.R
\name{csmprice}
\alias{csmprice}
\title{Computes option price with modified Corrado-Su model.}
\usage{
csmprice(type, spot, strike, time, rate, yield, sigma, mu3, mu4)
}
\arguments{
\item{type}{'call' for call option, any other value for put.}

\item{spot}{current stock price.}

\item{strike}{the strike price.}

\item{time}{time to option expiration in years.}

\item{rate}{the risk-free interest rate.}

\item{yield}{the dividends that are expected to be paid.}

\item{sigma}{volatility of the stock price.}

\item{mu3}{skewness.}

\item{mu4}{kurtosis.}
}
\value{
Vector with option prices.
}
\description{
Implementation of modified Corrado-Su formula for option pricing.
}
\section{Recycle rule}{


These arguments handle the recycle rule so vectors can be provided
and once those vectors differs in length the recycle rule is applied.
}

\examples{
spot <- 25
sigma <- 0.2
time <- 0.5
strike <- 26
rate <- 0.12
yield <- 0
skewness <- 1
kurtosis <- 4
csmprice(c("call", "put"), spot, strike, time, rate, yield, sigma, 0, 3) ==
  bsmprice(c("call", "put"), spot, strike, time, rate, yield, sigma)
csmprice(
  "put", spot, strike, time, rate, 0:100 / 100, sigma, skewness, kurtosis
)

}
\seealso{
\code{\link{csmw}} for w from modified Corrado-Su model.\\
\code{\link{csmd}} for d from modified Corrado-Su model.\\
\code{\link{csmq3}} for Q_3 from modified Corrado-Su model.\\
\code{\link{csmq4}} for Q_4 from modified Corrado-Su model.
}

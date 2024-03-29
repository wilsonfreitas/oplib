% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csm.R
\name{csmq3}
\alias{csmq3}
\title{Computes Q_3 from modified Corrado-Su model.}
\usage{
csmq3(spot, sigma, time, d1, w)
}
\arguments{
\item{spot}{current stock price.}

\item{sigma}{volatility of the stock price.}

\item{time}{time to option expiration in years.}

\item{d1}{d value from modified Corrado-Su model.}

\item{w}{constant from modified Corrado-Su model.}
}
\value{
Vector with Q_3.
}
\description{
Q_3 formula implementation from modified Corrado-Su model.
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
skewness <- 0
kurtosis <- 3
w <- csmw(sigma, time, skewness, kurtosis)
d <- csmd(spot, strike, time, rate, yield, sigma, w)
csmq3(spot, sigma, time, d, w)
csmq3(200:320 / 10, sigma, time, d, w)

}
\seealso{
\code{\link{csmw}} for w from modified Corrado-Su model.\\
\code{\link{csmd}} for d from modified Corrado-Su model.\\
\code{\link{csmq4}} for Q_4 from modified Corrado-Su model.\\
\code{\link{csmprice}} for modified Corrado-Su option pricing formula.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svi.R
\name{svi_fit}
\alias{svi_fit}
\title{SVI fit}
\usage{
svi_fit(variance, spot, strike, time, rate, yield = 0, initial_guess = NULL)
}
\arguments{
\item{variance}{variances for the given strikes}

\item{spot}{underlying spot price}

\item{strike}{options strikes}

\item{time}{time to maturity in years}

\item{rate}{risk free rate}

\item{yield}{yield rate}

\item{initial_guess}{Initial guess used in the optimization process}
}
\value{
A list with the SVI model parameters a, b, x, m, rho and sigma.
}
\description{
Fits a SVI model for the given smile as defined by strikes
and their respective variances. Time to maturity should also be provided
for the non-arbitrage restriction.
}
\seealso{
\code{\link{svi_var}} for SVI model variance calculation
}

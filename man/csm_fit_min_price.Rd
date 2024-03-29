% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csm.R
\name{csm_fit_min_price}
\alias{csm_fit_min_price}
\title{Corrado-Su Model Fit in Price}
\usage{
csm_fit_min_price(
  par,
  type,
  spot,
  strike,
  rate,
  yield,
  time,
  price,
  weights = 1,
  ...
)
}
\arguments{
\item{par}{a numeric value}

\item{type}{option type: \code{Call} or \code{Put}}

\item{spot}{a numeric value}

\item{strike}{a numeric value}

\item{rate}{a numeric value}

\item{yield}{a numeric value}

\item{time}{a numeric value}

\item{price}{a numeric value}

\item{weights}{a numeric value}

\item{...}{additional arguments passed to \code{optim}}
}
\description{
Function to fit Corrado-Su modelo that minimizes the error in price.
}

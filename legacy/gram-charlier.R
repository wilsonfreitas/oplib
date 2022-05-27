#' Check Gram Charlier limits from parameters
#'
#' This function verifies if kurtosis and skewness of a distribution satisfy
#' Gram-Charlier Kurtosis-Skewness limits.
#'
#' @param k kurtosis.
#' @param s skewness.
#' @param verbose TRUE or FALSE for whether to show or not warnings in console.
#'
#' @return
#' TRUE or FALSE, in or out of the limits, respectively.
#'
#' @seealso
#' \code{\link{check_gram_charlier_limits_4_series}}, \code{\link{adjmoments}}.
#'
#' @examples
#' check_gram_charlier_limits(4, .5)
#' check_gram_charlier_limits(6, -.5)
#' check_gram_charlier_limits(3, 0)
#'
#' @export
check_gram_charlier_limits <- function(k, s, verbose = FALSE) {
  data(gram_charlier_limits)
  rk <- range(gram_charlier_limits$points[, 1])
  if ((k < rk[1] || k > rk[2]) && verbose) {
    warning("kurtosis exceeded: ", k)
  }
  rs <- max(abs(gram_charlier_limits$points[, 2]))
  if (s > rs && verbose) {
    warning("skewness exceeded: ", s)
  }
  x <- gram_charlier_limits$fun(k) > s
  if (is.na(x)) FALSE else x
}

#' Check Gram Charlier limits from series
#'
#' This function verifies if kurtosis and skewness of a distribution satisfy
#' Gram-Charlier Kurtosis-Skewness limits.
#'
#' @param series return series to verify limits.
#' @param verbose TRUE or FALSE for whether to show or not warnings in console.
#'
#' @return
#' TRUE or FALSE, in or out of the limits, respectively.
#'
#' @seealso
#' \code{\link{check_gram_charlier_limits}}, \code{\link{adjmoments}}.
#'
#' @examples
#' check_gram_charlier_limits_4_series(rnorm(100000))
#'
#' @export
check_gram_charlier_limits_4_series <- function(series, verbose = FALSE) {
  k <- kurtosis(series, method = "moment")
  s <- abs(skewness(series, method = "moment"))
  check_gram_charlier_limits(k, s, verbose)
}

#' Out of the Boundaries Skewness and Kurtosis Adjustment
#'
#' In case skewness and kurtosis fall outside Gram-Charlier limits, use this
#' function to adjust them so they fall inside the borders of the region.
#'
#' @param s skewness
#' @param k kurtosis
#'
#' @return
#' Dataframe with adjusted skewness and kurtosis.
#'
#' @seealso
#' \code{\link{check_gram_charlier_limits}}, \code{\link{adjmoments}}.
#'
#' @examples
#' adjmoments(-200:200 / 100, 8)
#' adjmoments(0.5, 10:90 / 10)
#'
#' @export
adjmoments <- function(s, k) {
  data <- data.frame(s, k)
  if (nrow(data) == 0) {
    return()
  }

  data(gram_charlier_limits)
  edges <- as.data.frame(gram_charlier_limits$points)

  data <- within(data, {
    adj.skewness <- NA
    adj.kurtosis <- NA
  })
  na <- with(data, (is.na(s) | is.na(k)))

  kurtok <- with(data, !na & k >= 3 & k <= 7)
  data.kurtok <- data[kurtok, ]
  data.kurtok <- within(data.kurtok, {
    adj.kurtosis <- k
    skew.max <- approx(edges[, c(1, 3)], xout = k)[[2]]
    adj.skewness <- ifelse(abs(s) <= skew.max, s,
      sign(s) * skew.max
    )
  })
  data[kurtok, ] <- data.kurtok[, names(data)]

  sub.edges <- list(
    subset(edges[order(edges[, 3]), ], k <= 5),
    subset(edges[order(edges[, 3]), ], k > 5)
  )
  skewok <- with(data, !na & !kurtok & s >= -1 & s <= 1)
  data.skewok <- data[skewok, ]
  data.skewok <- within(data.skewok, {
    adj.skewness <- s
    adj.kurtosis <- ifelse(k < 3,
      approx(sub.edges[[1]][, c(3, 1)], xout = abs(s), ties = "ordered")[[2]],
      approx(sub.edges[[2]][, c(3, 1)], xout = abs(s), ties = "ordered")[[2]]
    )
  })
  data[skewok, ] <- data.skewok[, names(data)]

  nok <- with(data, !na & !kurtok & !skewok)
  data.nok <- data[nok, ]
  data.nok <- within(data.nok, {
    adj.kurtosis <- ifelse(k > 7, 7, 3)
    adj.skewness <- sign(s) * 0.5
  })
  data[nok, ] <- data.nok[, names(data)]

  return(data[, c("adj.skewness", "adj.kurtosis")])
}
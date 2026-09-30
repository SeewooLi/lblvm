ICOMP <- function(logL, Sigma) {
  ev <- eigen(Sigma, symmetric = TRUE)$values
  s <- length(ev)

  if (any(ev <= 0)) {
    penalty <- 0
    warning("Non-positive eigenvalues detected. ICOMP penalty set to 0.")
  } else {
    penalty <- (s / 2) * log(sum(ev) / s) - 0.5 * sum(log(ev))
  }

  return(-2 * logL + 2 * penalty)
}

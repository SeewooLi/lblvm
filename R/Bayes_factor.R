#' Bayes Factor from Laplace Approximation
#'
#' @param fit0 Model 0
#' @param fit1 Model 1
#'
#' @returns A list containing the following objects:
#' \describe{
#'   \item{log_marginal_M0}{Logarithm of the marginal probability of Model 0.}
#'   \item{log_marginal_M1}{Logarithm of the marginal probability of Model 1.}
#'   \item{logBF}{Logarithm of the Bayes factor}
#'   \item{BF}{The Bayes factor}
#' }
#' @export
#'
#' @details
#' Let \eqn{D} denote the observed data, \eqn{\tau} the parameter vector,
#' and \eqn{M_k} a model or hypothesis indexed by \eqn{k}. The marginal
#' probability of the data under \eqn{M_k} is
#' \deqn{
#'   I_k = P(D \mid M_k)
#'       = \int P(D \mid \tau, M_k)\pi(\tau \mid M_k)\,d\tau,
#' }
#' where \eqn{\pi(\tau \mid M_k)} is the prior density of \eqn{\tau}.
#' The marginal probability is approximated using a Laplace approximation:
#' \deqn{
#'   I_k \approx
#'   (2\pi)^{d/2}
#'   \lvert\hat{\Sigma}\rvert^{1/2}
#'   P(D \mid \hat{\tau}, M_k)
#'   \pi(\hat{\tau} \mid M_k),
#' }
#' where \eqn{d} is the dimension of the parameter vector, \eqn{\hat{\tau}}
#' is the mode of the posterior distribution, and \eqn{\hat{\Sigma}^{-1}} is the
#' Fisher information matrix. For example, `log_marginal_M0` denotes
#' \eqn{\log I_0}, the log marginal probability of the data under Model 0.
#'
#'
BF <- function(fit0, fit1) {
  cl <- match.call()
  name0 <- deparse(cl$fit0)
  name1 <- deparse(cl$fit1)

  logm0 <- laplace_log_marginal_invH(fit0$logL, fit0$l_prior_w, fit0$se$se_cov_raw)
  logm1 <- laplace_log_marginal_invH(fit1$logL, fit1$l_prior_w, fit1$se$se_cov_raw)

  logBF <- logm1 - logm0

  message("'", ifelse(logBF > 0, name1, name0), "' is more plausible.")

  return(list(
    log_marginal_M0 = logm0,
    log_marginal_M1 = logm1,
    logBF = logBF,
    BF = exp(logBF)
  ))
}

laplace_log_marginal_invH <- function(logLik_val, logPrior_val, Sigma) {
  k <- nrow(Sigma)

  logdetSigma <- as.numeric(
    determinant(Sigma, logarithm = TRUE)$modulus
  )

  return(logLik_val +
           logPrior_val +
           (k / 2) * log(2 * pi) + 0.5 * logdetSigma)
}


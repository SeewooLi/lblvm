#' Pseudo-Root-Mean-Square-Error-of-Approximation of Difference Test (\eqn{_pRMSEA_D})
#'
#' @param fit.ref A fitted object from the reference model.
#' @param fit.comp A fitted object from the comparison model.
#'
#' @returns A list containing \eqn{_pRMSEA_D} and its 90% confidence interval.
#'
#' @export
#'
rmsea_drm <- function(fit.ref, fit.comp){
  sample_size <- mean(colSums(!is.na(fit.ref$Options$data)))
  df <- fit.ref$eff_par - fit.comp$eff_par

  lrt <- -2*(fit.comp$logL - fit.ref$logL)
  lambda <- lrt - df # non-centrality parameter
  f0 <- max(lambda / (sample_size), 0)

  message("\r chi-square: ", sprintf("%.2f", lrt),",  df: ", sprintf("%.2f", df))

  if(lambda > 0){
    ci90 <- sqrt(
      (stats::qchisq(c(0.05, 0.95), df = df, ncp = lambda) - df)/ (df * sample_size)
    )
  } else {
    ci90 <- NULL
  }


  return(
    list(
      rmsea = sqrt(f0/df),
      ci90 = ci90
    )
  )
}

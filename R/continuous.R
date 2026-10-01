#' Run an exploratory factor analysis (EFA) for bounded-continuous responses
#'
#' @param data A matrix or data frame of continuous item responses.
#' @param dimension An integer specifying the number of dimensions of the
#'   latent variables.
#' @param range A numeric vector specifying the range of latent trait values
#'   for grid approximation in each dimension. Defaults to `c(-4, 4)`.
#' @param q An integer specifying the number of grid points per latent trait
#'   dimension for grid approximation. Defaults to `41`.
#' @param max_iter An integer specifying the maximum number of iterations for
#'   the expectation-maximization (EM) algorithm. Defaults to `200`.
#' @param threshold A numeric value specifying the convergence threshold on
#'   the maximum parameter change. Defaults to `0.0001`.
#'
#' @returns A list containing the following objects:
#' \describe{
#'   \item{par_est}{Estimated item parameters.}
#'   \item{se}{Standard errors of the item parameter estimates.}
#'   \item{fk}{Expected frequency (number of respondents) at each quadrature
#'     grid point.}
#'   \item{iter}{Number of EM algorithm iterations.}
#'   \item{quad}{Locations of the quadrature grid points.}
#'   \item{diff}{Value of the convergence criterion at the final iteration.}
#'   \item{prior}{Prior density at each quadrature grid point.}
#'   \item{posterior}{Posterior density at each quadrature grid point per respondent.}
#'   \item{Ak}{Posterior weight at each quadrature grid point.}
#'   \item{theta}{Expected a posteriori (EAP) estimates of the latent traits
#'     for each respondent.}
#'   \item{theta_se}{Standard errors of the EAP latent trait estimates.}
#'   \item{logL}{Approximated log-likelihood at the final iteration.}
#'   \item{f_means}{Estimated means of the latent factors.}
#'   \item{cov_mat}{Estimated covariance matrix of the latent factors.}
#'   \item{Options}{Input information used to call the function.}
#' }
#'
#' @export
#'
efa_cont <- function(data,dimension,range=c(-4,4),q=41,max_iter=200,
                    threshold=0.0001){

  # checking and reordering responses
  data <- check_data_cont(data)

  helper_matrices <- efa_helper_matrices_cont(dimension, data)

  fit <- fit_cont(data,
                  dimension = dimension,
                  contrast_m = helper_matrices$load_mat,
                  initialitem = helper_matrices$init_mat,
                  eq_constraint = helper_matrices$eq_constraint,
                  par_id = helper_matrices$par_id,
                  grouping = helper_matrices$grouping_index,
                  max_iter = max_iter,
                  threshold = threshold,
                  range = range,
                  q = q,
                  est_cov = FALSE
  )

  return(fit)
}



#' Run a confirmatory factor analysis (CFA) for bounded-continuous responses
#'
#' @param formula A symbolic description of the confirmatory factor analysis
#'   model. The formula specifies the relationships between latent factors and
#'   observed item responses, as well as any parameter constraints.
#' @param data A matrix or data frame of continuous item responses.
#' @param range A numeric vector specifying the range of latent trait values
#'   for grid approximation in each dimension. Defaults to `c(-4, 4)`.
#' @param q An integer specifying the number of grid points per latent trait
#'   dimension for grid approximation. Defaults to `41`.
#' @param max_iter An integer specifying the maximum number of iterations for
#'   the expectation-maximization (EM) algorithm. Defaults to `200`.
#' @param threshold A numeric value specifying the convergence threshold on
#'   the maximum parameter change. Defaults to `0.0001`.
#'
#' @returns A list containing the following objects:
#' \describe{
#'   \item{par_est}{Estimated item parameters.}
#'   \item{se}{Standard errors of the item parameter estimates.}
#'   \item{fk}{Expected frequency (number of respondents) at each quadrature
#'     grid point.}
#'   \item{iter}{Number of EM algorithm iterations.}
#'   \item{quad}{Locations of the quadrature grid points.}
#'   \item{diff}{Value of the convergence criterion at the final iteration.}
#'   \item{prior}{Prior density at each quadrature grid point.}
#'   \item{posterior}{Posterior density at each quadrature grid point per respondent.}
#'   \item{Ak}{Posterior weight at each quadrature grid point.}
#'   \item{theta}{Expected a posteriori (EAP) estimates of the latent traits
#'     for each respondent.}
#'   \item{theta_se}{Standard errors of the EAP latent trait estimates.}
#'   \item{logL}{Approximated log-likelihood at the final iteration.}
#'   \item{f_means}{Estimated means of the latent factors.}
#'   \item{cov_mat}{Estimated covariance matrix of the latent factors.}
#'   \item{Options}{Input information used to call the function.}
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' \dontrun{
#' # Items 1--5 on the first factor and the rest on the second factor
#' formula_string <- "
#' f1 ~ v1+v2+v3+v4+v5
#' f2 ~ v6+v7+v8+v9+v10
#'
#' v1.f1==v2.f1==v3.f2  # equality constraint
#'
#' v1.c <- 0  # fixed value
#' "
#'
#'
#' fit <- cfa_cont(formula_string, data)
#' }}
cfa_cont <- function(formula,data,range=c(-4,4),q=41,max_iter=200,
                     threshold=0.0001){

  # checking and reordering responses
  data <- check_data_cont(data)

  helper_matrices <- cfa_helper_matrices_cont(formula, data)

  fit <- fit_cont(data,
                  dimension = helper_matrices$d,
                  contrast_m = helper_matrices$load_mat,
                  initialitem = helper_matrices$init_mat,
                  eq_constraint = helper_matrices$eq_constraint,
                  par_id = helper_matrices$par_id,
                  grouping = helper_matrices$grouping_index,
                  max_iter = max_iter,
                  threshold = threshold,
                  range = range,
                  q = q,
                  est_cov = FALSE
  )

  return(fit)
}

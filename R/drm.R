#' Run an exploratory factor analysis (EFA) using the discretized response model (DRM)
#'
#' @param data A matrix or data frame of item responses.
#' @param dimension An integer specifying the number of dimensions of the
#'   latent variables.
#' @param t_prior A numeric vector specifying the standard deviations of the
#'   threshold priors.
#' @param range A numeric vector specifying the range of latent trait values
#'   for grid approximation in each dimension. Defaults to `c(-4, 4)`.
#' @param q An integer specifying the number of grid points per latent trait
#'   dimension for grid approximation. Defaults to `41`.
#' @param ngrid An integer specifying the number of grid points over the
#'   `0--1` range for grid approximation. Defaults to `1000`.
#' @param max_iter An integer specifying the maximum number of iterations for
#'   optimization. Defaults to `200` for `EM` and `400` for `MHRM`.
#' @param threshold A numeric value specifying the convergence threshold.
#'   Defaults to `0.000001` for `EM`, based on the gradient convergence
#'   criterion, and `0.001` for `MHRM`, based on the maximum parameter change.
#'   For `EM`, the gradient convergence criterion is defined as
#'   \eqn{\frac{\mathbf{g}'(\mathbf{H} + \mathbf{J})^{-1}\mathbf{g}}
#'   {|\log \mathcal{L}| + \delta}}.
#' @param eq_interval A logical value indicating whether category intervals
#'   are equally spaced. Defaults to `FALSE`.
#' @param estimation A character string specifying the estimation method:
#'   `EM` for the expectation-maximization algorithm or `MHRM` for the
#'   Metropolis-Hastings Robbins-Monro algorithm.
#'
#' @returns A list containing the following objects:
#' \describe{
#'   \item{par_est}{Item parameter estimates.}
#'   \item{iter}{Number of algorithm iterations.}
#'   \item{diff}{Value of the convergence criterion at the final iteration.}
#'   \item{logL}{Approximated log-likelihood.}
#'   \item{f_means}{Means of the latent factors.}
#'   \item{cov_mat}{Covariance matrix of the latent variables.}
#'   \item{Options}{Input information used to call the function.}
#' }
#'
#' When `estimation = "EM"`, the returned list also contains:
#' \describe{
#'   \item{se}{Standard errors of the parameter estimates.}
#'   \item{fk}{Expected frequency (number of observations) at each grid point.}
#'   \item{quad}{Locations of the quadrature grid points.}
#'   \item{prior}{Prior density at each grid point.}
#'   \item{posterior}{Posterior density at each grid point.}
#'   \item{Ak}{Posterior weight at each grid point.}
#'   \item{theta}{Expected a posteriori (EAP) estimate for each observation.}
#'   \item{theta_se}{Standard errors of the EAP estimates.}
#'   \item{l_prior_w}{Log prior density of the threshold parameters.}
#'   \item{AIC, BIC, ICOMP}{Akaike information criterion (AIC), Bayesian information criterion (BIC), and information complexity index, respectively.}
#'   \item{eff_par}{Effective number of parameters.}
#' }
#'
#' @export
#'
efa_drm <- function(data,dimension,t_prior=NULL,range=c(-4,4),q=41,ngrid=1000,max_iter=NULL,
                    threshold=NULL, eq_interval = FALSE, estimation="EM"){

  if(is.null(threshold)){
    if(estimation == "EM"){
      threshold <- 0.000001
    }else if(estimation == "MHRM"){
      threshold <- 0.001
    }
  }
  if(is.null(max_iter)){
    if(estimation == "EM"){
      max_iter <- 200
    }else if(estimation == "MHRM"){
      max_iter <- 400
    }
  }

  # checking and reordering responses
  data <- check_data_cat(data)

  if(is.null(t_prior)) t_prior <- t_prior_sd(N = nrow(data), category = apply(X = data, MARGIN = 2, FUN = max))

  helper_matrices <- efa_helper_matrices(dimension, data, eq_interval = eq_interval)

  fit <- dr_sim(data,
                dimension = dimension,
                contrast_m = helper_matrices$load_mat,
                initialitem = list(
                  helper_matrices$init_mat[,(1:(dimension+1))],
                  helper_matrices$init_mat[,-(1:(dimension+1))]
                ),
                eq_constraint = helper_matrices$eq_constraint,
                par_id = helper_matrices$par_id,
                grouping = helper_matrices$grouping_index,
                ngrid = ngrid,
                max_iter = max_iter,
                threshold = threshold,
                range = range,
                q = q,
                est_cov = FALSE,
                t_prior = t_prior,
                estimation=estimation
  )

  return(fit)
}

#' Run a confirmatory factor analysis (CFA) using the discretized response model (DRM)
#'
#' @param formula A character string specifying the predetermined model
#'   structure. See the examples for the supported specification.
#' @param data A matrix or data frame of item responses.
#' @param t_prior A numeric vector specifying the standard deviations of the
#'   threshold priors.
#' @param range A numeric vector specifying the range of latent trait values
#'   for grid approximation in each dimension. Defaults to `c(-4, 4)`.
#' @param q An integer specifying the number of grid points per latent trait
#'   dimension for grid approximation. Defaults to `41`.
#' @param ngrid An integer specifying the number of grid points over the
#'   `0--1` range for grid approximation. Defaults to `1000`.
#' @param max_iter An integer specifying the maximum number of iterations for
#'   optimization. Defaults to `200` for `EM` and `400` for `MHRM`.
#' @param threshold A numeric value specifying the convergence threshold.
#'   Defaults to `0.000001` for `EM`, based on the gradient convergence
#'   criterion, and `0.001` for `MHRM`, based on the maximum parameter change.
#'   For `EM`, the gradient convergence criterion is defined as
#'   \eqn{\frac{\mathbf{g}'(\mathbf{H} + \mathbf{J})^{-1}\mathbf{g}}
#'   {|\log \mathcal{L}| + \delta}}.
#' @param eq_interval A logical value indicating whether category intervals
#'   are equally spaced. Defaults to `FALSE`.
#' @param estimation A character string specifying the estimation method:
#'   `EM` for the expectation-maximization algorithm or `MHRM` for the
#'   Metropolis-Hastings Robbins-Monro algorithm.
#' @param est_cov A logical value indicating whether to estimate the covariance
#'   matrix of the latent variables \eqn{\theta}. Defaults to `TRUE`.
#'
#' @returns A list containing the following objects:
#' \describe{
#'   \item{par_est}{Item parameter estimates.}
#'   \item{iter}{Number of algorithm iterations.}
#'   \item{diff}{Value of the convergence criterion at the final iteration.}
#'   \item{logL}{Approximated log-likelihood.}
#'   \item{f_means}{Means of the latent factors.}
#'   \item{cov_mat}{Covariance matrix of the latent variables.}
#'   \item{Options}{Input information used to call the function.}
#' }
#'
#' When `estimation = "EM"`, the returned list also contains:
#' \describe{
#'   \item{se}{Standard errors of the parameter estimates.}
#'   \item{fk}{Expected frequency (number of observations) at each grid point.}
#'   \item{quad}{Locations of the quadrature grid points.}
#'   \item{prior}{Prior density at each grid point.}
#'   \item{posterior}{Posterior density at each grid point.}
#'   \item{Ak}{Posterior weight at each grid point.}
#'   \item{theta}{Expected a posteriori (EAP) estimate for each observation.}
#'   \item{theta_se}{Standard errors of the EAP estimates.}
#'   \item{l_prior_w}{Log prior density of the threshold parameters.}
#'   \item{AIC, BIC, ICOMP}{Akaike information criterion (AIC), Bayesian information criterion (BIC), and information complexity index, respectively.}
#'   \item{eff_par}{Effective number of parameters.}
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
#' v1.f1==v2.f1==v3.f1==v4.f1==v5.f1  # equality constraint
#'
#' v1.c <- 0  # fixed value
#' "
#'
#'
#' fit <- cfa_drm(formula_string, data)
#' }}
cfa_drm <- function(formula,data,t_prior=NULL,range=c(-4,4),q=41,ngrid=1000,max_iter=NULL,
                    threshold=NULL, eq_interval = FALSE, estimation="EM",est_cov = TRUE){

  if(is.null(threshold)){
    if(estimation == "EM"){
      threshold <- 0.000001
    }else if(estimation == "MHRM"){
      threshold <- 0.001
    }
  }
  if(is.null(max_iter)){
    if(estimation == "EM"){
      max_iter <- 500
    }else if(estimation == "MHRM"){
      max_iter <- 1000
    }
  }

  # checking and reordering responses
  data <- check_data_cat(data)

  if(is.null(t_prior)) t_prior <- t_prior_sd(N = nrow(data), category = apply(X = data, MARGIN = 2, FUN = max))

  helper_matrices <- cfa_helper_matrices(formula, data, eq_interval = eq_interval)

  fit <- dr_sim(data,
                dimension = helper_matrices$d,
                contrast_m = helper_matrices$load_mat,
                initialitem = list(
                  helper_matrices$init_mat[,(1:(helper_matrices$d+1))],
                  helper_matrices$init_mat[,-(1:(helper_matrices$d+1))]
                ),
                eq_constraint = helper_matrices$eq_constraint,
                par_id = helper_matrices$par_id,
                grouping = helper_matrices$grouping_index,
                ngrid = ngrid,
                max_iter = max_iter,
                threshold = threshold,
                range = range,
                q = q,
                est_cov = est_cov,
                t_prior = t_prior,
                estimation=estimation
  )

  return(fit)
}

#' Run a multiple-group (MG) confirmatory factor analysis (EFA) using the discretized response model (DRM)
#'
#' @param formula A character string specifying the predetermined model
#'   structure. See the examples for the supported specification.
#' @param data A matrix or data frame of item responses.
#' @param group A factor specifying group membership.
#' @param t_prior A numeric vector specifying the standard deviations of the
#'   threshold priors.
#' @param range A numeric vector specifying the range of latent trait values
#'   for grid approximation in each dimension. Defaults to `c(-4, 4)`.
#' @param q An integer specifying the number of grid points per latent trait
#'   dimension for grid approximation. Defaults to `41`.
#' @param ngrid An integer specifying the number of grid points over the
#'   `0--1` range for grid approximation. Defaults to `1000`.
#' @param max_iter An integer specifying the maximum number of iterations for
#'   optimization. Defaults to `200` for `EM` and `400` for `MHRM`.
#' @param threshold A numeric value specifying the convergence threshold.
#'   Defaults to `0.000001` for `EM`, based on the gradient convergence
#'   criterion, and `0.001` for `MHRM`, based on the maximum parameter change.
#'   For `EM`, the gradient convergence criterion is defined as
#'   \eqn{\frac{\mathbf{g}'(\mathbf{H} + \mathbf{J})^{-1}\mathbf{g}}
#'   {|\log \mathcal{L}| + \delta}}.
#' @param eq_interval A logical value indicating whether category intervals
#'   are equally spaced. Defaults to `FALSE`.
#' @param est_cov A logical value indicating whether to estimate the covariance
#'   matrix of the latent variables \eqn{\theta}. Defaults to `TRUE`.
#'
#' @returns A list containing the following objects:
#' \describe{
#'   \item{par_est}{Item parameter estimates.}
#'   \item{se}{Standard errors and related quantities for the parameter
#'     estimates, including raw and adjusted standard errors and inverse
#'     Hessian matrices.}
#'   \item{fk}{Expected frequency (number of observations) at each grid point.}
#'   \item{iter}{Number of algorithm iterations.}
#'   \item{quad}{Locations of the quadrature grid points.}
#'   \item{diff}{Value of the convergence criterion at the final iteration.}
#'   \item{prior}{Prior density at each grid point.}
#'   \item{posterior}{Posterior density at each grid point.}
#'   \item{Ak}{Posterior weight at each grid point.}
#'   \item{l_prior_w}{Log prior density of the threshold parameters.}
#'   \item{logL}{Approximated log-likelihood.}
#'   \item{AIC, BIC, ICOMP}{Akaike information criterion (AIC), Bayesian
#'     information criterion (BIC), and information complexity index,
#'     respectively.}
#'   \item{cov_mat}{Covariance matrix of the latent variables.}
#'   \item{f_mean}{Means of the latent variables.}
#'   \item{eff_par}{Effective number of parameters.}
#'   \item{Options}{Input information used to call the function.}
#' }
#'
#' @export
#'
cfa_drm_mg <- function(formula,data,group,t_prior=NULL,range=c(-4,4),q=41,ngrid=1000,max_iter=200,
                      threshold=0.000001, eq_interval = FALSE,est_cov = TRUE){

  # checking and reordering responses
  data <- check_data_cat(data)

  # grouping
  grouped <- transform_group(data, group)
  data <- grouped$data
  group <- grouped$group_int
  if(is.null(t_prior)) t_prior <- grouped$t_prior

  helper_matrices <- group_helper_matrices(formula,
                                           data,
                                           eq_interval = eq_interval,
                                           ngroup = length(unique(group)),
                                           group_label = grouped$group_char)

  fit <- dr_group(data,
                  group = group,
                  dimension = helper_matrices$d,
                  contrast_m = helper_matrices$load_mat,
                  initialitem = list(
                    helper_matrices$init_mat[,(1:(helper_matrices$d+1))],
                    helper_matrices$init_mat[,-(1:(helper_matrices$d+1))]
                  ),
                  eq_constraint = helper_matrices$eq_constraint,
                  par_id = helper_matrices$par_id,
                  grouping = helper_matrices$grouping_index,
                  ngrid = ngrid,
                  max_iter = max_iter,
                  threshold = threshold,
                  range = range,
                  q = q,
                  est_cov = est_cov,
                  t_prior = t_prior
  )

  names(fit$f_mean) <- names(fit$cov_mat) <- grouped$group_char
  fit$Options$group_label <- grouped$group_char
  return(fit)
}

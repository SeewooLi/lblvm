#' Title
#'
#' @param data
#' @param dimension
#' @param range
#' @param q
#' @param max_iter
#' @param threshold
#'
#' @returns
#' @export
#'
#' @examples
efa_cont <- function(data,dimension,range=c(-4,4),q=41,max_iter=200,
                    threshold=0.0001){

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



#' cfa
#'
#' @param formula
#' @param data
#' @param range
#' @param q
#' @param max_iter
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
cfa_cont <- function(formula,data,range=c(-4,4),q=41,max_iter=200,
                     threshold=0.0001){

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

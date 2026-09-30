################################################################################
# Exported Functions
################################################################################

#' Standard Deviations of the Threshold Priors
#'
#' @param N Sample size.
#' @param category The number of response options.
#' @param n The number of items that share the same prior via constraints. Defaults to `1`.
#'
#' @returns A vector resulting from the `1/sqrt(N/category*n)*sqrt(2*pi^2/3)` operation.
#' @export
#'
#' @examples
#' \donttest{
#' t_prior_sd(
#'   N = 100,
#'   category = 5,
#'   n = 1
#' )
#' }
t_prior_sd <- function(N, category, n=1){
  1/sqrt(N/category*n)*sqrt(2*pi^2/3)
}


#' Transformation of the Threshold Parameters to the Unit Scale
#'
#' @param x A vector of the raw threshold parameters.
#'
#' @returns Transformed thresholds.
#' @export
#'
#' @examples
#' \donttest{
#' cut_trans(c(0,0,0,0))
#' }
cut_trans <- function(x){
  x <- x[!is.na(x)]
  return(
    cumsum(
      c(1,exp(x[-length(x)]))/sum(c(1,exp(x)))
    )
  )
}

################################################################################
# Internal Functions
################################################################################


one_hot_3d <- function(M) {
  N <- nrow(M)            # number of persons
  I <- ncol(M)            # number of items
  K <- max(M, na.rm = TRUE) + 1         # number of categories (0-based)

  arr <- array(0, dim = c(I, N, K))

  for (i in 1:I) {
    mat <- matrix(0, N, K)
    mat[cbind(1:N, M[,i] + 1)] <- 1   # +1 because R indices start at 1
    arr[i,,] <- mat
  }

  return(arr)
}

jac_soft <- function(t) {
  npar <- length(t)
  e <- c(1, exp(t))
  D <- sum(e)

  N <- rep(0, npar)
  tmp <- 0
  for(i in 1:npar){
    tmp <- tmp + e[i]
    N[i] <- tmp
  }

  J <- matrix(0, npar, npar)

  for (i in 1:npar) {
    for (j in 1:npar) {
      if (j < i) {
        J[i,j] <- e[j] * (D - N[i]) / D^2
      } else {
        J[i,j] <- - N[i] * e[j] / D^2
      }
    }
  }
  return(J)
}

count_cat <- function(x) length(unique(x))
extract_cat <- function(x) sort(unique(x))

reorder_vec <- function(x){
  match(x, table = extract_cat(x))-1
}

reorder_mat <- function(x) apply(x, MARGIN = 2, FUN = reorder_vec)

P_2PL <- function(theta, a, b) 1/(1 + exp(-(as.vector(theta %*% a) + b)))

calculate_J <- function(psi,threshold) {
  nth <- length(psi)
  w_vars <- exp(c(-sum(psi),psi))
  P_vars <- w_vars / sum(w_vars)

  Diff_P <- P_vars[-1] - P_vars[1]
  J <- matrix(0, nth, nth)
  mask <- lower.tri(J)
  J <- mask * rep(Diff_P, each = nth) - (threshold %o% Diff_P) - ((!mask) * P_vars[1])
  return(J)
}


transform_group <- function(data, group) {

  if (nrow(data) != length(group)) {
    stop("Length of `group` must match number of rows in `data`.")
  }

  # Drop NA
  keep <- !is.na(group)
  data_filtered  <- data[keep, , drop = FALSE]
  group_filtered <- group[keep]
  if(sum(!keep) != 0) message(sum(!keep), " NA observations in `group` are excluded.")

  # INTEGER CASE
  if (is.integer(group_filtered)) {

    return(list(
      data = data_filtered,
      group_int = group_filtered,
      group_char = NULL
    ))
  }

  # CHARACTER CASE
  if (is.character(group_filtered)) {

    f <- factor(group_filtered)

    group_int  <- as.integer(f)          # observation-level
    group_char <- levels(f)              # unique labels

    return(list(
      data = data_filtered,
      group_int  = group_int,
      group_char = group_char
    ))
  }

  # FACTOR CASE
  if (is.factor(group_filtered)) {

    group_int  <- as.integer(group_filtered)   # observation-level
    group_char <- levels(group_filtered)       # unique labels

    return(list(
      data = data_filtered,
      group_int  = group_int,
      group_char = group_char
    ))
  }

  stop("`group` must be integer, character, or factor.")
}

P_DRM <- function(theta, a, b, nu, ncats=NULL, cut_score=NULL, return_mu = FALSE){
  p <- P_2PL(theta = theta, a = a, b = b)
  if(is.null(cut_score) & is.null(ncats)){
    stop("Specify either ncat or cut_score.")
  }else if(is.null(cut_score)){
    cut_score <- (1:(ncats-1))/ncats
  }else{
    ncats <- length(cut_score)+1
  }

  probs <- matrix(nrow = nrow(theta), ncol = ncats-1)

  for(i in 1:length(cut_score)){
    probs[,i] <- stats::pbeta(q = cut_score[i],
                       shape1 = p*nu,
                       shape2 = (1-p)*nu)
  }
  if(return_mu){
    return(
      list(prob = cbind(probs,1)-cbind(0,probs),
           mu = p)
    )
  }else {
    return(cbind(probs,1)-cbind(0,probs))
  }
}

diff_to_matrix <- function(par_id, vec) {
  mat <- as.matrix(par_id)
  replaced <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat),
                     dimnames = dimnames(mat))

  for (i in seq_along(vec)) {
    replaced[mat == i] <- vec[i]
  }

  replaced[is.na(replaced)] <- 0
  return(replaced)
}

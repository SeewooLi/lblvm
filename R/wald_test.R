vectorize_params <- function(param_mat, par_id) {

  values <- as.vector(param_mat)
  ids <- as.vector(par_id)

  keep <- !is.na(ids)
  values <- values[keep]
  ids <- ids[keep]

  split_vals <- split(values, ids)

  theta <- sapply(split_vals, function(v) {
    if (length(unique(v)) > 1) {
      warning("Inconsistent values found for same parameter ID")
    }
    v[1]
  })

  theta <- theta[order(as.numeric(names(theta)))]

  theta
}

build_wald_components <- function(constraints, par_id, theta, Sigma) {

  constraints <- trimws(constraints)
  lines <- unlist(strsplit(constraints, "\\s+"))
  lines <- lines[lines != ""]

  get_index <- function(token) {
    parts <- strsplit(token, "\\.")[[1]]
    item <- parts[1]
    par <- parts[2]
    par_id[item, par]
  }

  all_indices <- c()
  chains <- list()

  for (line in lines) {
    tokens <- unlist(strsplit(line, "=="))
    tokens <- trimws(tokens)

    idx <- sapply(tokens, get_index)

    chains[[length(chains) + 1]] <- idx
    all_indices <- c(all_indices, idx)
  }

  unique_idx <- sort(unique(all_indices))
  idx_map <- stats::setNames(seq_along(unique_idx), unique_idx)

  n_constraints <- sum(sapply(chains, function(x) length(x) - 1))
  k <- length(unique_idx)

  R <- matrix(0, nrow = n_constraints, ncol = k)

  row_counter <- 1

  for (idx in chains) {
    for (i in 1:(length(idx) - 1)) {
      c1 <- idx_map[as.character(idx[i])]
      c2 <- idx_map[as.character(idx[i + 1])]

      R[row_counter, c1] <- 1
      R[row_counter, c2] <- -1

      row_counter <- row_counter + 1
    }
  }

  theta_sub <- theta[unique_idx]
  Sigma_sub <- Sigma[unique_idx, unique_idx]

  list(
    R = R,
    theta = theta_sub,
    Sigma = Sigma_sub,
    idx = unique_idx
  )
}


#' Perform a Wald Test of Parameter Equality
#'
#' @param constraints A character string specifying equality constraints
#'   between model parameters.
#' @param fit A fitted model object returned by a compatible `lblvm`
#'   estimation function.
#'
#' @returns A list containing the following objects:
#' \describe{
#'   \item{W}{Wald test statistic.}
#'   \item{df}{Degrees of freedom for the Wald test.}
#'   \item{p_value}{P-value associated with the Wald test statistic.}
#' }
#'
#' @export
Wald_drm <- function(constraints, fit) {
  pars_mat <- cbind(
    fit$par_est[[1]],
    log(fit$par_est[[2]][, 1]),
    fit$par_est[[2]][, -1]
  )
  colnames(pars_mat) <- colnames(fit$Options$par_id)

  par_vec <- vectorize_params(
    pars_mat,
    fit$Options$par_id
  )

  out <- build_wald_components(
    constraints,
    fit$Options$par_id,
    par_vec,
    fit$se$se_cov_raw
  )

  R <- out$R
  theta <- out$theta
  Sigma <- out$Sigma

  W <- t(R %*% theta) %*%
    solve(R %*% Sigma %*% t(R)) %*%
    (R %*% theta)

  W <- as.vector(W)

  if (length(W) != 1) {
    stop("The statistic is not a scalar.")
  }

  df <- nrow(R)
  p <- 1 - stats::pchisq(W, df)

  list(
    W = W,
    df = df,
    p_value = p
  )
}


#' Calculate Local RMSEA
#'
#' @param constraints A character string specifying equality constraints
#'   between model parameters.
#' @param fit A fitted model object returned by a compatible `lblvm`
#'   estimation function.
#'
#' @returns A list containing the following objects:
#' \describe{
#'   \item{rmsea}{The local RMSEA estimate based on the Wald statistic.}
#'   \item{ci90}{A 90\% confidence interval for the local RMSEA. Returns
#'     `NULL` when the estimated noncentrality parameter is non-positive.}
#' }
#'
#' @export
local_rmsea <- function(constraints, fit) {
  sample_size <- mean(colSums(!is.na(fit$Options$data)))

  Wald_results <- Wald_drm(constraints, fit)
  df <- Wald_results$df
  W <- Wald_results$W

  lambda <- W - df
  f0 <- max(lambda / sample_size, 0)

  message(
    "\r chi-square: ",
    sprintf("%.2f", W),
    ",  df: ",
    sprintf("%.2f", df)
  )

  if (lambda > 0) {
    ci90 <- sqrt(
      (
        stats::qchisq(
          c(0.05, 0.95),
          df = df,
          ncp = lambda
        ) - df
      ) / (df * sample_size)
    )
  } else {
    ci90 <- NULL
  }

  list(
    rmsea = sqrt(f0 / df),
    ci90 = ci90
  )
}

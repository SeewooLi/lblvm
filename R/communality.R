#' Communality and Unique Dispersion
#'
#' @param fit A returned object from a model-fitting function.
#'
#' @returns A matrix with communality and unique dispersion indices.
#' @export
#'
communality <- function(fit){
  if(inherits(fit, "sg")){
    thresholds <- apply(fit$par_est[[2]][,-1], 1, cut_trans, simplify = FALSE)

    results <- matrix(nrow = ncol(fit$Options$data), ncol = 5)
    for(i in 1:nrow(results)){
      results[i, 1:2] <- null_dispersion(fit$Options$data[,i], thresholds[[i]])
    }

    results[, 3] <- fit$par_est[[2]][,1]

    results[, 5] <- results[, 2]/results[, 3]
    results[, 4] <- 1 - results[, 5]

    colnames(results) <- c("mu", "nu_null", "nu_model", "communality", "unique dispersion")
  } else if(inherits(fit, "mg")){
    group <- fit$Options$group
    n_group <- length(unique(group))
    n_item <- ncol(fit$Options$data)

    results <- list()

    for(g in 1:n_group){
      index <- (1:n_item) + n_item*(g-1)
      thresholds <- apply(fit$par_est[[2]][index,-1], 1, cut_trans, simplify = FALSE)

      tmp <- matrix(nrow = ncol(fit$Options$data), ncol = 5)
      for(i in 1:nrow(tmp)){
        tmp[i, 1:2] <- null_dispersion(fit$Options$data[group == g,i], thresholds[[i]])
      }

      tmp[, 3] <- fit$par_est[[2]][index,1]

      tmp[, 5] <- tmp[, 2]/tmp[, 3]
      tmp[, 4] <- 1 - tmp[, 5]

      colnames(tmp) <- c("mu", "nu_null", "nu_model", "communality", "unique dispersion")

      results[[g]] <- tmp
    }
    names(results) <- fit$Options$group_label
  }


  return(results)
}

null_dispersion <- function(data, cut_score=NULL){
  data <- reorder_vec(data)
  freq <- as.vector(table(data))

  ncats <- length(freq)

  if(is.null(cut_score)) {
    cut_score <- (1:ncats)/ncats
  } else {
    cut_score <- c(cut_score, 1)
  }

  pars <- c(0.5, 1)
  iter <- 0
  repeat{
    iter <- iter + 1

    mu <- pars[1]
    xi <- pars[2]
    nu <- exp(xi)

    p0 <- stats::pbeta(q = cut_score, shape1 = mu*nu,shape2 = (1-mu)*nu)
    p0 <- p0 - c(0, p0[-length(p0)])

    ngrid <- 1000

    beta_grid <- seq(1/2/ngrid, 1-1/2/ngrid, length=ngrid)

    ind_cat <- as.numeric(cut(beta_grid,breaks = c(0,cut_score),labels = 1:ncats))

    p_ <- stats::dbeta(beta_grid, shape1 = nu*mu, shape2 = nu*(1-mu))
    s1 <- log(beta_grid) - digamma(nu*mu)
    s2 <- log(1- beta_grid) - digamma(nu*(1-mu))
    num_mu <- p_*nu*(s1 - s2)
    num_xi <- nu*p_*(digamma(nu) + mu*s1 + (1-mu)*s2)


    l1m <- l1x <- c()
    for(ct in 1:ncats){
      l1m[ct] <- sum(num_mu[ind_cat==ct])/length(beta_grid)
      l1x[ct] <- sum(num_xi[ind_cat==ct])/length(beta_grid)
    }

    L1 <- c(
      freq %*% (l1m/p0),
      freq %*% (l1x/p0)
    )

    L2 <- sum(freq) * t(cbind(l1m, l1x)) %*% diag(1/p0) %*% cbind(l1m, l1x)

    diff <- L1 %*% solve(L2)

    pars <- pars + diff
    if(sum(abs(diff)) < 0.000001 | iter > 30) break
  }
  pars[2] <- exp(pars[2])

  return(as.vector(pars))
}

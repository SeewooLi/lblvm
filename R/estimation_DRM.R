sample_mhrm <- function(data, theta_c, item, cut_score, f_cov=NULL, sd=1){
  d <- ncol(theta_c)
  if(is.null(f_cov)){
    Sigma_inv <- diag(d)
  }else{
    Sigma_inv <- solve(f_cov)
  }

  AR <- 0

  theta_p <- theta_c
  for(i in 1:d){
    theta_p[,i] <- theta_p[,i] + stats::rnorm(nrow(theta_p), 0, sd)
    lr <- pmin(llik_ratio_mhrm(data, theta_c, theta_p, item, cut_score, Sigma_inv), 0)
    not_accepted <- (log(stats::runif(nrow(theta_p))) > lr)

    theta_p[not_accepted,i] <- theta_c[not_accepted,i]
    theta_c <- theta_p

    AR <- AR + (1 - mean(not_accepted)) / d
  }

  logL <- llik_mhrm(data, theta_c, item, cut_score, Sigma_inv)
  return(
    list(
      theta=theta_c,
      logL=logL,
      AR=AR
      )
  )
}

Mstep_sim <- function(E, item, par_id=NULL, grouping=NULL, ngrid = 1000,
                  max_iter=7, threshold=1e-7, prior, calculate_m=FALSE, data=NULL, group=NULL){
  grid <- E$grid
  d <- ncol(grid)

  n_item <- nrow(item[[1]])
  npar <-  max(par_id, na.rm = TRUE)
  se <- list()
  f_mat <- list()

  e.response <- E$e.response
  q <- nrow(grid)
  beta_grid <- seq(1/2/ngrid, 1-1/2/ngrid, length=ngrid)

  iter <- 0
  div <- 3

  if(!is.null(group)){
    n_group <- max(group)
    n_item_g <- n_item/n_group
  }

  repeat{
    iter <- iter + 1

    L1 <- rep(0, npar)
    L2 <- matrix(0, nrow = npar, ncol = npar)
    IMm <- matrix(0, nrow = npar, ncol = npar)
    se <- matrix(0, nrow = npar, ncol = npar)
    inv_L2 <- matrix(0, nrow = npar, ncol = npar)
    l_prior_w <- 0
    diff <- L1
    grad_norm <- 0
    eff_par <- 0

    # stacking up the gradients and information
    for(i in 1:n_item){
      if(!is.null(group)){
        g <- (i - 1) %/% n_item_g + 1
        item_col <- (i - 1) %% n_item_g + 1
      }

      L1L2 <- L1L2_sim(c(item[[1]][i,],item[[2]][i,]),
                       e.response[i,,],
                       grid, beta_grid, d,
                       posterior=if(calculate_m)E$posterior else NULL,
                       response=if(calculate_m){
                         if(!is.null(group)) data[group == g, item_col] else data[,i]
                         }else NULL,
                       calculate_m)

      ind <- par_id[i,]
      ind2 <- !is.na(ind)
      ind <- ind[ind2]

      L1[ind] <- L1[ind] + L1L2$Grad[ind2]
      L2[ind,ind] <- L2[ind,ind] + L1L2$IM[ind2,ind2]
      if(calculate_m){
        IMm[ind,ind] <- IMm[ind,ind] + L1L2$IMm[ind2,ind2]
      }
    }
    L1[!is.finite(L1)] <- sign(L1[!is.finite(L1)]) * 1
    original_L2 <- L2 # to later calculate the effective number of parameters

    # apply prior
    partial_id <- par_id[,(d+3):ncol(par_id)]
    item_t <- item[[2]][,-1]
    item_t <- sweep(item_t, 1, rowMeans(cbind(item_t, 0), na.rm = TRUE), FUN = "-") # centering around psi mean
    ind <- as.vector(partial_id)
    unique_ind <- unique(stats::na.omit(ind))
    current_t <- sapply(unique_ind, function(i){
      mean(item_t[ind==i], na.rm=TRUE)
    })
    item_number <- sapply(unique_ind, function(v) {
      which(partial_id == v, arr.ind = TRUE)[1, 1]
    })
    if(length(current_t) != 0){
      l_prior_w <- mvtnorm::dmvnorm(current_t, mean = rep(0, length(current_t)), sigma = diag(prior[item_number]^2), log = TRUE) # log-prior-weight
      L1[unique_ind] <- L1[unique_ind] - (current_t - 0)/(prior[item_number]^2)
      L2[unique_ind,unique_ind] <- L2[unique_ind,unique_ind] + diag(1/(prior[item_number]^2))
    }

    # calculate parameter changes
    for(i in 1:length(grouping)){
      inv_L2[grouping[[i]],grouping[[i]]] <- solve(L2[grouping[[i]],grouping[[i]]])
      if(calculate_m){
        missing_I <- IMm[grouping[[i]],grouping[[i]]]
        se[grouping[[i]],grouping[[i]]] <- solve(L2[grouping[[i]],grouping[[i]]] - missing_I)
        eff_par <- eff_par + sum(diag(
          (original_L2[grouping[[i]],grouping[[i]]] - missing_I) %*% inv_L2[grouping[[i]],grouping[[i]]]
        ))
      } else {
        eff_par <- eff_par + sum(diag(
          original_L2[grouping[[i]],grouping[[i]]] %*% inv_L2[grouping[[i]],grouping[[i]]]
        ))
      }

      diff[grouping[[i]]] <- as.vector(inv_L2[grouping[[i]],grouping[[i]]] %*% L1[grouping[[i]]])

      grad_norm <- grad_norm + diff[grouping[[i]]] %*% L1[grouping[[i]]]
    }

    # convert diff-vector to matrix
    diff <- diff_to_matrix(par_id, diff)

    # update parameters
    if( max(abs(diff)) > div){
      stepsize <- max(sum(abs(diff)),2)
      item[[1]] <- item[[1]] + diff[,1:(d+1)]/stepsize
      item[[2]] <- item[[2]] + diff[,(d+2):ncol(diff)]/stepsize
    } else {
      item[[1]] <- item[[1]] + diff[,1:(d+1)]/2
      item[[2]] <- item[[2]] + diff[,(d+2):ncol(diff)]/2
      div <- max(abs(diff))
    }

    if( div <= threshold | iter > max_iter) break
  }
  if(calculate_m){
    se_delta <- se
    if(length(current_t) != 0){
      for(i in 1:nrow(partial_id)){
        index <- partial_id[i,]
        if(all(!is.na(index))){
          J <- jac_soft(item[[2]][i,-1])
          se_delta[index,index] <- J %*% se[index,index] %*% t(J)
        }
      }
    }
    suppressWarnings({
      se_delta_reduced <- diff_to_matrix(par_id, sqrt(diag(se_delta)))
      se_raw <- diff_to_matrix(par_id, sqrt(diag(se)))
    })
    return(list(
      se_cov=se_delta,
      se_cov_raw=se,
      se=se_delta_reduced,
      se_raw=se_raw,
      eff_par=eff_par
      ))
  } else {
    return(list(
      item=item,
      grad_norm=grad_norm,
      l_prior_w=l_prior_w,
      H=original_L2,
      H_prior_added=L2,
      eff_par=eff_par))
  }
}

update_gH <- function(data_array, thetalist, item, par_id=NULL, grouping=NULL, ngrid = 1000,
                      prev_g, prev_H, prior, stepsize=1, f_cov=NULL){
  d <- ncol(thetalist[[1]])

  if(is.null(f_cov)){
    f_cov <- diag(d)
  }

  n_item <- nrow(item[[1]])
  npar <-  max(par_id, na.rm = TRUE)

  beta_grid <- seq(1/2/ngrid, 1-1/2/ngrid, length=ngrid)

  n_sample <- length(thetalist)

  # START the g, H calculation for each item
  L1 <- rep(0, npar)
  L2 <- matrix(0, nrow = npar, ncol = npar)
  diff <- L1
  grad_norm <- 0
  eff_par <- 0
  mu_hat <- score_Sigma <- H_Sigma <- 0

  # stacking up the gradients and information
  for(k in 1:n_sample){
    theta <- thetalist[[k]]
    for(i in 1:n_item){
      L1L2 <- L1L2_sim(c(item[[1]][i,],item[[2]][i,]),
                       data_array[i,,],
                       theta, beta_grid, d)

      ind <- par_id[i,]
      ind2 <- !is.na(ind)
      ind <- ind[ind2]

      L1[ind] <- L1[ind] + L1L2$Grad[ind2] / n_sample
      L2[ind,ind] <- L2[ind,ind] + L1L2$IM[ind2,ind2] / n_sample
    }
    # L1[!is.finite(L1)] <- sign(L1[!is.finite(L1)]) * 1

    n <- nrow(theta)
    mu_hat <- mu_hat + colMeans(theta) / n_sample
    Sigma_hat <- stats::cov(theta)
    Sigma_inv <- solve(f_cov)
    score_Sigma <- score_Sigma + (n/2) * as.vector(Sigma_inv - Sigma_inv %*% Sigma_hat %*% Sigma_inv) / n_sample
    H_Sigma <- H_Sigma -(n/2) * kronecker(Sigma_inv, Sigma_inv) / n_sample
  }

  # update gradient and Hessian
  prev_g[[1]] <- prev_g[[1]] + stepsize * (L1 - prev_g[[1]])
  L2 <- prev_H[[1]] <- prev_H[[1]] + stepsize * (L2 - prev_H[[1]])
  original_L2 <- L2 # to later calculate the effective number of parameters

  prev_g[[2]] <- prev_g[[2]] + stepsize * (score_Sigma - prev_g[[2]])
  prev_H[[2]] <- prev_H[[2]] + stepsize * (H_Sigma - prev_H[[2]])
  f_cov <- f_cov + stepsize * solve(prev_H[[2]], score_Sigma)

  # apply prior
  partial_id <- par_id[,(d+3):ncol(par_id)]
  item_t <- item[[2]][,-1]
  item_t <- sweep(item_t, 1, rowMeans(cbind(item_t, 0), na.rm = TRUE), FUN = "-") # centering around psi mean
  ind <- as.vector(partial_id)
  unique_ind <- unique(stats::na.omit(ind))
  current_t <- sapply(unique_ind, function(i){
    mean(item_t[ind==i], na.rm=TRUE)
  })

  item_number <- sapply(unique_ind, function(v) {
    which(partial_id == v, arr.ind = TRUE)[1, 1]
  })
  if(length(current_t) != 0){
    L1[unique_ind] <- L1[unique_ind] - (current_t - 0)/(prior[item_number]^2)
    L2[unique_ind,unique_ind] <- L2[unique_ind,unique_ind] + diag(1/(prior[item_number]^2))
  }

  # calculate parameter changes
  for(i in 1:length(grouping)){
    inv_L2 <- solve(L2[grouping[[i]],grouping[[i]]])
    diff[grouping[[i]]] <- as.vector(inv_L2 %*% L1[grouping[[i]]])

    grad_norm <- grad_norm + diff[grouping[[i]]] %*% L1[grouping[[i]]]
    eff_par <- eff_par + sum(diag(original_L2[grouping[[i]],grouping[[i]]] %*% inv_L2))
  }

  # convert diff-vector to matrix
  diff <- diff_to_matrix(par_id, diff)

  # update parameters
  item[[1]] <- item[[1]] + stepsize * diff[,1:(d+1)]
  item[[2]] <- item[[2]] + stepsize * diff[,(d+2):ncol(diff)]

  return(list(
    item=item,
    g_approx=prev_g,
    H_approx=prev_H,
    L1=L1,
    H=original_L2,
    H_prior_added=L2,
    grad_norm=grad_norm,
    eff_par=eff_par,
    diff=diff,
    factor_means=mu_hat,
    f_cov=f_cov))
}


L1L2_sim <- function(par, e.response, grid, beta_grid, d, posterior=NULL, response=NULL, calculate_m = FALSE){
  # par <- par[!is.na(par)]
  npar <- length(par)
  cut_score <- cut_trans(par[(d+3):npar])
  f <- rowSums(e.response)
  ncats <- length(cut_score)+1
  ncut <- ncats-1
  ind_cat <- as.numeric(cut(beta_grid,breaks = c(0,cut_score,1),labels = 1:ncats))

  # 1st derivatives
  l1s <- rep(list(0), npar)

  nu <- exp(par[d+2])

  p0 <- P_DRM_cpp(grid, par[1:d], par[d+1], nu, cut_score = cut_score, return_mu = TRUE)
  mu <- p0$mu
  p0 <- p0$prob

  l1s <- compute_l1_cpp(grid, nu, mu, beta_grid, ncats, ind_cat)

  # thresholds
  beta_vals <- matrix(nrow = nrow(grid), ncol=ncut)
  for(b in 1:ncut){
    beta_vals[,b] <- stats::dbeta(cut_score[b], shape1 = nu*mu, shape2 = nu*(1-mu))
  }

  for(b in 1:ncut){
    ind <- c(rep(0,b),rep(1,ncut-b))
    temp <- sweep(beta_vals,MARGIN = 2,STATS = (ind-cut_score)*cut_score[1]*exp(par[d+2+b]),FUN = "*")
    l1s[[b+d+2]] <- cbind(temp,0)-cbind(0,temp)
  }

  res_cpp <- Grad_IM_cpp(l1s,
                         e.response,
                         p0,
                         f,
                         posterior = posterior,
                         response  = response,
                         calculate_m)

  if(calculate_m){
    return(list(
      Grad = res_cpp$Grad,
      IM = res_cpp$IM,
      f = sum(f),
      IMm = res_cpp$IMm
    ))
  }else{
    return(list(
      Grad = res_cpp$Grad,
      IM = res_cpp$IM,
      f = sum(f)
    ))
  }
}

################################################################################
# THE MAIN DRIVER FUNCTION
################################################################################
dr_sim <- function(data, dimension=NULL, range=c(-4,4), q=11, t_prior=0.25, initialitem=NULL, ngrid=1000,
               max_iter=200, threshold=0.0000001,contrast_m=NULL, eq_constraint=NULL,
               par_id=NULL, grouping=NULL, est_cov=FALSE, f_cov=NULL, estimation="EM",threshold2=0.01){
  # data to matrix
  data <- as.matrix(data)
  if(is.null(colnames(data))) colnames(data) <- paste0("v", 1:ncol(data))


  # setting grid on the latent space
  x <- seq(range[1], range[2], length.out=q)
  grid_list <- replicate(dimension, x, simplify = FALSE)
  grid <- as.matrix(do.call(expand.grid, grid_list))

  # saving initial values
  init <- initialitem

  Options = list(initialitem=initialitem,dimension=dimension,data=data,range=range,
                 q=q, max_iter=max_iter, threshold=threshold,eq_constraint=eq_constraint,
                 par_id=par_id, grouping=grouping)

  # preparation for EM
  I <- initialitem
  if(is.null(f_cov)) f_cov <- diag(dimension)
  sds <- sqrt(diag(f_cov))
  prior <- mvtnorm::dmvnorm(grid,
                            mean = rep(0, dimension),
                            sigma = f_cov
  )
  prior <- prior/sum(prior)
  iter <- 0
  diff <- 1

  # preparation for MHRM
  if(estimation == "MHRM"){
    n_sample <- 3 # the number of samples in the stochastic imputation step

    thetalist <- list()
    theta <- matrix(0, nrow = nrow(data), ncol = dimension)
    data_array <- one_hot_3d(data)
    g_approx <- list(0,0)
    H_approx <- list(0,0)

    proposal_sd <- 2.38 # kernel for MH sampling
  }

  # prior for the thresholds
  if(length(t_prior) == 1) t_prior <- rep(t_prior, nrow(initialitem[[1]]))

  phase <- 1

  # EM
  repeat{
    iter <- iter + 1

    cut_score <- apply(initialitem[[2]][,-1], 1, cut_trans, simplify = FALSE)
    #--------------------------------
    # EM algorithm
    #--------------------------------
    if(estimation == "EM"){
      E <- Estep_cpp(data,
                     cbind(initialitem[[1]],initialitem[[2]][,1]),
                     grid,
                     prior,
                     cut_score)
      dim(E$e.response) <- c(nrow(initialitem[[1]]), nrow(grid), max(data, na.rm = TRUE) + 1)

      M <- Mstep_sim(E,
                     initialitem,
                     par_id,
                     grouping,
                     prior = if(phase==1) t_prior/10 else t_prior,
                     ngrid = if(phase==1) 100 else ngrid)

      initialitem <- M$item
      grad_norm <- M$grad_norm
      l_prior_w <- M$l_prior_w
      eff_par_c <- M$eff_par
      logL <- E$logL

      # estimating cov matrix
      factor_means <- as.vector(E$Ak%*%E$grid)
      cov_mat <- t(E$grid) %*% sweep(E$grid, 1, E$Ak, FUN = "*") - factor_means %*% t(factor_means)
      sds <- sqrt(diag(cov_mat))
      cov_mat <- cov_mat / (sds %o% sds)

      # update prior
      if(est_cov){
        f_cov <- cov_mat
        prior <- mvtnorm::dmvnorm(grid,
                                  mean = rep(0, dimension),
                                  sigma = f_cov)
        prior <- prior/sum(prior)
      }
      #--------------------------------
      # MHRM algorithm
      #--------------------------------
    } else if (estimation =="MHRM"){
      logL <- 0
      AR <- 0 # acceptance rate
      for(burnin in 1:n_sample){
        thetaUpdate <- sample_mhrm(data, theta, cbind(initialitem[[1]],initialitem[[2]][,1]), cut_score, f_cov, proposal_sd)
        theta <- thetalist[[burnin]] <- thetaUpdate$theta
        logL <- logL + thetaUpdate$logL / n_sample
        AR <- AR + thetaUpdate$AR / n_sample
      }
      if(AR < .3) proposal_sd <- proposal_sd * 0.8
      if(AR > .5) proposal_sd <- proposal_sd * 1.25

      updated_gH <- update_gH(data_array = data_array,
                              thetalist = thetalist,
                              item = initialitem,
                              par_id = par_id,
                              grouping = grouping,
                              ngrid = if(phase==1) 100 else ngrid,
                              prev_g = g_approx,
                              prev_H = H_approx,
                              prior = if(phase==1) t_prior/10 else t_prior,
                              stepsize = if(phase==1) 1 else 1/(iter-20),
                              f_cov = f_cov)
      initialitem <- updated_gH$item

      g_approx <- updated_gH$g_approx
      H_approx <- updated_gH$H_approx
      grad_norm <- updated_gH$grad_norm
      if(est_cov){
        f_cov <- updated_gH$f_cov
        sds <- sqrt(diag(f_cov))
        f_cov <- f_cov / (sds %o% sds)
      }
      eff_par <- updated_gH$eff_par
      factor_means <- rep(0, dimension) # updated_gH$factor_means
    }


    # adjust parameters
    if(dimension == 1){
      initialitem[[1]][,1] <- initialitem[[1]][,1] * sds
      initialitem[[1]][,2] <- initialitem[[1]][,2] + initialitem[[1]][,1] * factor_means
    }else{
      initialitem[[1]][,1:dimension] <- sweep(initialitem[[1]][,1:dimension], 2, sds, FUN = "*")
      initialitem[[1]][,dimension+1] <- initialitem[[1]][,dimension+1] + initialitem[[1]][,1:dimension] %*% factor_means
    }
    initialitem[[1]][is.na(par_id[,1:(dimension+1)])] <- init[[1]][is.na(par_id[,1:(dimension+1)])]
    if(!is.null(eq_constraint) & (max(eq_constraint)>0)){
      for(i in 1:max(eq_constraint)){
        initialitem[[1]][which(eq_constraint[,1:(dimension+1)]==i)] <- mean(initialitem[[1]][which(eq_constraint[,1:(dimension+1)]==i)])
        initialitem[[2]][which(eq_constraint[,(dimension+2):ncol(eq_constraint)]==i)] <- mean(initialitem[[2]][which(eq_constraint[,(dimension+2):ncol(eq_constraint)]==i)])
      }
    }

    diff <- c(
      max(abs(I[[1]]-initialitem[[1]]), na.rm = TRUE),
      max(abs(I[[2]]-initialitem[[2]]), na.rm = TRUE)
      )
    max_par <- par_id[which.max(cbind(abs(I[[1]]-initialitem[[1]]),abs(I[[2]]-initialitem[[2]])))]
    diff <- max(diff)
    grad_conv <- grad_norm/(abs(logL)+1e-06)
    if(is.na(grad_conv)) grad_conv <- 1
    I <- initialitem
    if(estimation == "EM"){
      message("\r","\r","Phase", phase, ", EM cycle = ",iter,", logL = ", sprintf("%.2f", logL), ", GradConv = ", sprintf("%.7f",grad_conv),", Max-Change = ", sprintf("%.5f", diff), " at (", paste(max_par, collapse = ","), ")  ", sep="",appendLF=FALSE)
      utils::flush.console()
      if((phase==2) & (iter >= max_iter | (grad_conv < threshold & diff < threshold2))) break
    } else if(estimation == "MHRM"){
      message("\r","\r","Phase", phase, ", MHRM cycle = ",iter, ", (AR|",sprintf("%.1f", AR),"), logL = ", sprintf("%.2f", logL),", Max-Change = ", sprintf("%.5f", diff), " at (", paste(max_par, collapse = ","), ")  ", sep="",appendLF=FALSE)
      utils::flush.console()
      if((phase==2) & (iter >= max_iter | diff < threshold)) break
    }
    if(diff<threshold2 & phase==1) phase <- 2
  }

  # preparation for outputs
  dimnames(initialitem[[1]]) <- list(colnames(data),c(paste0("f",1:dimension), "c"))
  dimnames(initialitem[[2]]) <- list(colnames(data),c("nu", paste("t", 1:(ncol(initialitem[[2]])-1), sep="")))
  initialitem[[2]][,1] <- exp(initialitem[[2]][,1])


  if (estimation =="EM"){
    theta <- E$posterior%*%E$grid
    theta_se <- sqrt(E$posterior%*%(E$grid^2)-theta^2)

    message("\n","\r","Calculating standard error...", sep="",appendLF=FALSE)
    utils::flush.console()
    M <- Mstep_sim(E,
                   I,
                   par_id,
                   grouping,
                   prior = t_prior,
                   ngrid = ngrid,
                   calculate_m = TRUE,
                   max_iter = 1,
                   data = data)
    eff_par <- M$eff_par
    return(structure(
      list(par_est=initialitem,
           se=M[c("se", "se_raw", "se_cov", "se_cov_raw")],
           fk=E$freq,
           iter=iter,
           quad=grid,
           diff=diff,
           prior=E$prior,
           posterior=E$posterior,
           Ak=E$Ak,
           theta = theta,
           theta_se = theta_se,
           l_prior_w=l_prior_w,
           logL= logL,
           AIC= -2*logL+2*eff_par_c,
           BIC= -2*logL+log(mean(colSums(!is.na(data))))*eff_par_c,
           ICOMP= ICOMP(logL, M$se_cov_raw),
           f_means=factor_means,
           cov_mat = f_cov,
           eff_par = eff_par_c,
           Options = Options # specified argument values
      ),
      class = c("sg", "dr", "list")
    )
    )
  } else if (estimation =="MHRM"){
    return(structure(
      list(par_est=initialitem,
           iter=iter,
           diff=diff,
           logL= logL,
           f_means=factor_means,
           cov_mat = f_cov,
           eff_par = eff_par,
           Options = Options # specified argument values
      ),
      class = c("sg", "dr", "list")
    )
    )
  }
}

dr_group <- function(
    data, group, dimension=NULL, range=c(-4,4), q=11, t_prior=0.25, initialitem=NULL, ngrid=1000,
    max_iter=200, threshold=0.000001, contrast_m=NULL, eq_constraint=NULL,
    par_id=NULL, grouping=NULL, est_cov=FALSE, f_cov=NULL, threshold2=0.01){

  if(is.null(colnames(data))) colnames(data) <- paste0("v", 1:ncol(data))

  # the number of groups
  ngroup <- length(unique(group))
  nitem_g <- nrow(initialitem[[1]])/ngroup

  # setting grid on the latent space
  x <- seq(range[1], range[2], length.out=q)
  grid_list <- replicate(dimension, x, simplify = FALSE)
  grid <- as.matrix(do.call(expand.grid, grid_list))

  # saving initial values
  init <- initialitem

  Options = list(initialitem=initialitem,dimension=dimension,data=data,range=range,
                 q=q, max_iter=max_iter, threshold=threshold,eq_constraint=eq_constraint,
                 par_id=par_id, grouping=grouping, group=group, threshold2-threshold2)

  # preparation for EM
  I <- initialitem
  prior0 <- mvtnorm::dmvnorm(grid,
                             mean = rep(0, dimension),
                             sigma = if(is.null(f_cov)) diag(dimension) else f_cov
  )
  prior0 <- prior0/sum(prior0)
  prior <- list()
  for(g in 1:ngroup){
    prior[[g]] <- prior0
  }

  iter <- 0
  diff <- 1

  if(length(t_prior) == 1) t_prior <- rep(t_prior, nrow(initialitem[[1]]))


  phase <- 1

  # EM
  repeat{
    iter <- iter + 1

    # E-step
    Ak <- list()
    freq <- list()
    posterior <- matrix(nrow = nrow(data), ncol = nrow(grid))
    logL <- 0
    e.response <- c()
    for(g in 1:ngroup){
      item_index <- ((g-1)*nitem_g + 1):(g*nitem_g)
      cut_score <- apply(initialitem[[2]][item_index,-1], 1, cut_trans, simplify = FALSE)

      E <- Estep_cpp(data[group==g,],
                      cbind(initialitem[[1]][item_index,],initialitem[[2]][item_index,1]),
                      grid,
                      prior[[g]],
                      cut_score)

      dim(E$e.response) <- c(nitem_g, nrow(grid), max(data, na.rm = TRUE) + 1)

      e.response <- abind::abind(e.response, E$e.response, along = 1)

      Ak[[g]] <- E$Ak
      freq[[g]] <- E$freq
      posterior[group == g, ] <- E$posterior
      logL <- logL + E$logL
    }

    # M-step
    M <- Mstep_sim(list(grid = grid, e.response = e.response),
                   initialitem,
                   par_id,
                   grouping,
                   prior = if(phase == 1) t_prior/10 else t_prior,
                   ngrid = if(phase == 1) 100 else ngrid)
    initialitem <- M$item
    l_prior_w <- M$l_prior_w
    eff_par_c <- M$eff_par

    # estimating cov matrix
    factor_means <- list()
    cov_mat <- list()
    sds <- list()
    for(g in 1:ngroup){
      factor_means[[g]] <- as.vector(Ak[[g]] %*% grid)
      cov_mat[[g]] <- t(grid) %*% sweep(grid, 1, Ak[[g]], FUN = "*") - factor_means[[g]] %*% t(factor_means[[g]])
      sds[[g]] <- sqrt(diag(cov_mat[[g]]))
    }


    # update prior
    if(est_cov){
      f_cov <- list()
      prior <- list()
      for(g in 1:ngroup){
        f_cov[[g]] <- cov_mat[[g]]
        prior[[g]] <- mvtnorm::dmvnorm(grid,
                                   mean = factor_means[[g]] - factor_means[[1]],
                                   sigma = f_cov[[g]] / (sds[[1]] %o% sds[[1]])
                                   )
        prior[[g]] <- prior[[g]]/sum(prior[[g]])
      }
    } else{
      f_cov <- list()
      prior <- list()
      for(g in 1:ngroup){
        f_cov[[g]] <- cov_mat[[g]]
        prior[[g]] <- mvtnorm::dmvnorm(grid,
                                       mean = factor_means[[g]] - factor_means[[1]],
                                       sigma = diag(diag(f_cov[[g]] / (sds[[1]] %o% sds[[1]])))
        )
        prior[[g]] <- prior[[g]]/sum(prior[[g]])
      }
    }

    # adjust parameters
    if(dimension == 1){
      initialitem[[1]][,1] <- initialitem[[1]][,1] * sds[[1]]
      initialitem[[1]][,2] <- initialitem[[1]][,2] + initialitem[[1]][,1] * factor_means[[1]]
    }else{
      initialitem[[1]][,1:dimension] <- sweep(initialitem[[1]][,1:dimension], 2, sds[[1]], FUN = "*")
      initialitem[[1]][,dimension+1] <- initialitem[[1]][,dimension+1] + initialitem[[1]][,1:dimension] %*% factor_means[[1]]
    }
    initialitem[[1]][is.na(par_id[,1:(dimension+1)])] <- init[[1]][is.na(par_id[,1:(dimension+1)])]
    if(!is.null(eq_constraint) & (max(eq_constraint)>0)){
      for(i in 1:max(eq_constraint)){
        initialitem[[1]][which(eq_constraint[,1:(dimension+1)]==i)] <- mean(initialitem[[1]][which(eq_constraint[,1:(dimension+1)]==i)])
        initialitem[[2]][which(eq_constraint[,(dimension+2):ncol(eq_constraint)]==i)] <- mean(initialitem[[2]][which(eq_constraint[,(dimension+2):ncol(eq_constraint)]==i)])
      }
    }

    diff <- c(
      max(abs(I[[1]]-initialitem[[1]]), na.rm = TRUE),
      max(abs(I[[2]]-initialitem[[2]]), na.rm = TRUE)
    )

    max_par <- par_id[which.max(cbind(abs(I[[1]]-initialitem[[1]]),abs(I[[2]]-initialitem[[2]])))]
    diff <- max(diff)
    grad_conv <- M$grad_norm/(abs(logL)+1e-06)
    if(is.na(grad_conv)) grad_conv <- 1
    I <- initialitem
    message("\r","\r","Phase", phase, ", EM cycle = ",iter,", logL = ", sprintf("%.2f", logL), ", GradConv = ", sprintf("%.7f",grad_conv),", Max-Change = ", sprintf("%.5f", diff), " at (", paste(max_par, collapse = ","), ")  ", sep="",appendLF=FALSE)
    utils::flush.console()
    if((phase==2) & (iter >= max_iter | (grad_conv < threshold & diff < threshold2))) break
    if(diff<threshold2 & phase==1) phase <- 2
  }

  # preparation for outputs
  dimnames(initialitem[[1]]) <- list(rownames(initialitem[[1]]),c(paste0("f",1:dimension), "c"))
  dimnames(initialitem[[2]]) <- list(rownames(initialitem[[1]]),c("nu", paste("t", 1:(ncol(initialitem[[2]])-1), sep="")))
  initialitem[[2]][,1] <- exp(initialitem[[2]][,1])

  message("\n","\r","Calculating standard error...", sep="",appendLF=FALSE)
  utils::flush.console()

  M <- Mstep_sim(list(grid = grid, e.response = e.response, posterior = posterior),
                 I,
                 par_id,
                 grouping,
                 prior = t_prior,
                 ngrid = ngrid,
                 calculate_m = TRUE,
                 max_iter = 1,
                 data = data,
                 group = group)
  eff_par <- M$eff_par

  return(structure(
    list(par_est=initialitem,
         se=M[c("se", "se_raw", "se_cov", "se_cov_raw")],
         fk=freq,
         iter=iter,
         quad=grid,
         diff=diff,
         prior=prior,
         posterior=posterior,
         Ak=Ak,
         l_prior_w=l_prior_w,
         logL= logL,
         AIC= -2*logL+2*eff_par_c,
         BIC= -2*logL+log(mean(colSums(!is.na(data))))*eff_par_c,
         ICOMP= ICOMP(logL, M$se_cov_raw),
         cov_mat = f_cov,
         f_mean = factor_means,
         eff_par = eff_par_c,
         Options = Options # specified argument values
    ),
    class = c("mg", "dr", "list")
  )
  )
}

Mstep_cont <- function(E, item, par_id=NULL, grouping=NULL,
                      max_iter=7, threshold=1e-7, data, calculate_m=FALSE){
  grid <- E$grid
  d <- ncol(grid)

  n_item <- nrow(item)
  npar <-  max(par_id, na.rm = TRUE)
  se <- list()

  posterior <- E$posterior
  q <- nrow(grid)

  iter <- 0
  div <- 3

  repeat{
    iter <- iter + 1

    L1 <- diff <- rep(0, npar)
    L2 <- matrix(0, nrow = npar, ncol = npar)
    IMm <- matrix(0, nrow = npar, ncol = npar)
    se <- matrix(0, nrow = npar, ncol = npar)
    inv_L2 <- matrix(0, nrow = npar, ncol = npar)

    # stacking up the gradients and information
    for(i in 1:n_item){
      L1L2 <- cont_L1L2_cpp(item = item[i,],
                            grid = grid,
                            data = data[,i],
                            Pk = posterior,
                            calculate_m = calculate_m)

      ind <- par_id[i,]
      ind2 <- !is.na(ind)
      ind <- ind[ind2]

      L1[ind] <- L1[ind] + L1L2$gradient[ind2]
      L2[ind,ind] <- L2[ind,ind] - L1L2$hessian[ind2,ind2]
      if(calculate_m){
        IMm[ind,ind] <- IMm[ind,ind] + L1L2$IMm[ind2,ind2]
      }
    }
    L1[!is.finite(L1)] <- sign(L1[!is.finite(L1)]) * 1

    # calculate parameter changes
    for(i in 1:length(grouping)){

      if(calculate_m){
        se[grouping[[i]],grouping[[i]]] <- solve(L2[grouping[[i]],grouping[[i]]] - IMm[grouping[[i]],grouping[[i]]])
      }

      diff[grouping[[i]]] <- as.vector(solve(L2[grouping[[i]],grouping[[i]]]) %*% L1[grouping[[i]]])
    }

    # convert diff-vector to matrix
    diff <- diff_to_matrix(par_id, diff)

    # update parameters
    if( max(abs(diff)) > div){
      stepsize <- max(sum(abs(diff)),2)
      item <- item + diff/stepsize
    } else {
      item <- item + diff/2
      div <- max(abs(diff))
    }

    if( div <= threshold | iter > max_iter) break
  }
  if(calculate_m){
    se_cov <- se
    se <- diff_to_matrix(par_id, sqrt(diag(se)))

    return(list(
      se_cov=se_cov,
      se=se
    ))
  } else {
    return(list(
      item=item,
      H=L2
      ))
  }
}

################################################################################
# THE MAIN DRIVER FUNCTION
################################################################################
fit_cont <- function(data, dimension=NULL, range=c(-4,4), q=11, t_prior=0.25, initialitem=NULL,
                   max_iter=200, threshold=0.0000001,contrast_m=NULL, eq_constraint=NULL,
                   par_id=NULL, grouping=NULL, est_cov=FALSE, f_cov=NULL){
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

  # EM
  repeat{
    iter <- iter + 1

    E <- Estep_Cont_cpp(data,
                        initialitem,
                        grid,
                        prior)

    M <- Mstep_cont(E,
                    initialitem,
                    par_id,
                    grouping,
                    data = data)

    initialitem <- M$item
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

    # adjust parameters
    if(dimension == 1){
      initialitem[,1] <- initialitem[,1] * sds
      initialitem[,2] <- initialitem[,2] + initialitem[,1] * factor_means
    }else{
      initialitem[,1:dimension] <- sweep(initialitem[,1:dimension], 2, sds, FUN = "*")
      initialitem[,dimension+1] <- initialitem[,dimension+1] + initialitem[,1:dimension] %*% factor_means
    }

    if(!is.null(eq_constraint) & (max(eq_constraint)>0)){
      for(i in 1:max(eq_constraint)){
        initialitem[which(eq_constraint[,1:(dimension+1)]==i)] <- mean(initialitem[which(eq_constraint[,1:(dimension+1)]==i)])
      }
    }

    diff <- max(abs(I-initialitem), na.rm = TRUE)

    max_par <- par_id[which.max(abs(I-initialitem))]
    diff <- max(diff)
    I <- initialitem

    message("\r","\r","EM cycle = ",iter,", logL = ", sprintf("%.2f", logL), ", Max-Change = ", sprintf("%.5f", diff), " at (", paste(max_par, collapse = ","), ")  ", sep="",appendLF=FALSE)
    utils::flush.console()
    if(iter >= max_iter | diff < threshold) break
  }

  # preparation for outputs
  dimnames(initialitem) <- list(colnames(data),c(paste0("f",1:dimension), "c", "nu"))
  initialitem[,dimension+2] <- exp(initialitem[,dimension+2])

  theta <- E$posterior%*%E$grid
  theta_se <- sqrt(E$posterior%*%(E$grid^2)-theta^2)

  message("\n","\r","Calculating standard error...", sep="",appendLF=FALSE)
  utils::flush.console()
  M <- Mstep_cont(E,
                  I,
                  par_id,
                  grouping,
                  data = data,
                  max_iter = 1,
                  calculate_m = TRUE)
  return(structure(
    list(par_est=initialitem,
         se=M[c("se_cov", "se")],
         fk=E$freq,
         iter=iter,
         quad=grid,
         diff=diff,
         prior=E$prior,
         posterior=E$posterior,
         Ak=E$Ak,
         theta = theta,
         theta_se = theta_se,
         logL= logL,
         f_means=factor_means,
         cov_mat = f_cov,
         Options = Options # specified argument values
    ),
    class = c("sg", "cont", "list")
  )
  )
}

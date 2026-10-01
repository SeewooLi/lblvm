# tests/testthat/test-model_fitting.R

# Helper function for data generation
P_2PL <- function(theta, a, b) 1 / (1 + exp(-(as.vector(theta %*% a) + b)))

P_DRM <- function(theta, a, b, nu, cut_score){
  p <- P_2PL(theta = theta, a = a, b = b)

  ncats <- length(cut_score)+1

  probs <- matrix(nrow = nrow(theta), ncol = ncats-1)

  for(i in 1:length(cut_score)){
    probs[,i] <- stats::pbeta(q = cut_score[i],
                              shape1 = p*nu,
                              shape2 = (1-p)*nu)
  }
  return(cbind(probs,1)-cbind(0,probs))
}

# Helper to generate test dataset
make_test_data_cont <- function() {
  set.seed(1)
  dimension <- 2

  N <- 300
  nitem <- 6

  item <- c(rep(1.5, nitem/2), rep(0, nitem/2))
  item <- cbind(item, rev(item), 0, log(10))

  theta <- mvtnorm::rmvnorm(
    N,
    mean = rep(0, dimension),
    sigma = diag(dimension)
  )

  data <- matrix(ncol = nitem, nrow = N)
  for (i in 1:nitem) {
    mu <- P_2PL(theta, item[i, 1:dimension], item[i, dimension + 1])
    data[, i] <- stats::rbeta(N, shape1 = mu * exp(item[i, dimension + 2]), shape2 = (1 - mu) * exp(item[i, dimension + 2]))
  }

  colnames(data) <- paste0("v", 1:nitem)
  return(data)
}

make_test_data_cat <- function() {
  set.seed(1)
  dimension <- 2

  N <- 300
  nitem <- 6

  item <- c(rep(1.5, nitem/2), rep(0, nitem/2))
  item <- cbind(item, rev(item), 0, log(10))

  theta <- mvtnorm::rmvnorm(
    N,
    mean = rep(0, dimension),
    sigma = diag(dimension)
  )

  data <- matrix(ncol = nitem, nrow = N)
  for(i in 1:nitem){
    p_mat <- P_DRM(theta, item[i,1:dimension], item[i,dimension+1], exp(item[i,dimension+2]), c(.2,.4,.6,.8))
    for(p in 1:N){
      data[p,i] <- sample(x = 1:5, size = 1, prob = p_mat[p,])
    }
  }

  colnames(data) <- paste0("v", 1:nitem)
  return(data)
}

test_that("efa_cont runs without error and returns valid par_est", {
  data <- make_test_data_cont()

  # Expect no errors during execution
  expect_error(
    fit <- efa_cont(data, 2, max_iter = 1),
    regexp = NA
  )

  # Check that par_est exists and has no NA or NaN values
  expect_true("par_est" %in% names(fit))
  expect_false(any(is.na(fit$par_est)))
})

test_that("cfa_cont runs without error and returns valid par_est", {
  data <- make_test_data_cont()

  formula_string <- "
  f1 ~ v1+v2+v3
  f2 ~ v4+v5+v6

  v1.f1==v2.f1==v3.f1

  v1.c <- 0
  "

  # Expect no errors during execution
  expect_error(
    fit <- cfa_cont(formula_string, data, max_iter = 1),
    regexp = NA
  )

  # Check that par_est exists and has no NA or NaN values
  expect_true("par_est" %in% names(fit))
  expect_false(any(is.na(fit$par_est)))
})

test_that("efa_cont runs without error and returns valid par_est", {
  data <- make_test_data_cat()

  # Expect no errors during execution
  expect_error(
    fit <- efa_drm(data, 2, max_iter = 1),
    regexp = NA
  )

  # Check that par_est exists and has no NA or NaN values
  expect_true("par_est" %in% names(fit))
  expect_false(any(is.na(fit$par_est)))
})

test_that("cfa_cont runs without error and returns valid par_est", {
  data <- make_test_data_cat()

  formula_string <- "
  f1 ~ v1+v2+v3
  f2 ~ v4+v5+v6

  v1.f1==v2.f1==v3.f1

  v1.c <- 0
  "

  # Expect no errors during execution
  expect_error(
    fit <- cfa_drm(formula_string, data, max_iter = 1),
    regexp = NA
  )

  # Check that par_est exists and has no NA or NaN values
  expect_true("par_est" %in% names(fit))
  expect_false(any(is.na(fit$par_est)))
})

# Run a multiple-group (MG) confirmatory factor analysis (EFA) using the discretized response model (DRM)

Run a multiple-group (MG) confirmatory factor analysis (EFA) using the
discretized response model (DRM)

## Usage

``` r
cfa_drm_mg(
  formula,
  data,
  group,
  t_prior = NULL,
  range = c(-4, 4),
  q = 41,
  ngrid = 1000,
  max_iter = 200,
  threshold = 1e-06,
  eq_interval = FALSE,
  est_cov = TRUE
)
```

## Arguments

- formula:

  A character string specifying the predetermined model structure. See
  the examples for the supported specification.

- data:

  A matrix or data frame of item responses.

- group:

  A factor specifying group membership.

- t_prior:

  A numeric vector specifying the standard deviations of the threshold
  priors.

- range:

  A numeric vector specifying the range of latent trait values for grid
  approximation in each dimension. Defaults to `c(-4, 4)`.

- q:

  An integer specifying the number of grid points per latent trait
  dimension for grid approximation. Defaults to `41`.

- ngrid:

  An integer specifying the number of grid points over the `0--1` range
  for grid approximation. Defaults to `1000`.

- max_iter:

  An integer specifying the maximum number of iterations for
  optimization. Defaults to `200` for `EM` and `400` for `MHRM`.

- threshold:

  A numeric value specifying the convergence threshold. Defaults to
  `0.000001` for `EM`, based on the gradient convergence criterion, and
  `0.001` for `MHRM`, based on the maximum parameter change. For `EM`,
  the gradient convergence criterion is defined as
  \\\frac{\mathbf{g}'(\mathbf{H} + \mathbf{J})^{-1}\mathbf{g}} {\|\log
  \mathcal{L}\| + \delta}\\.

- eq_interval:

  A logical value indicating whether category intervals are equally
  spaced. Defaults to `FALSE`.

- est_cov:

  A logical value indicating whether to estimate the covariance matrix
  of the latent variables \\\theta\\. Defaults to `TRUE`.

## Value

A list containing the following objects:

- par_est:

  Item parameter estimates.

- se:

  Standard errors and related quantities for the parameter estimates,
  including raw and adjusted standard errors and inverse Hessian
  matrices.

- fk:

  Expected frequency (number of observations) at each grid point.

- iter:

  Number of algorithm iterations.

- quad:

  Locations of the quadrature grid points.

- diff:

  Value of the convergence criterion at the final iteration.

- prior:

  Prior density at each grid point.

- posterior:

  Posterior density at each grid point.

- Ak:

  Posterior weight at each grid point.

- l_prior_w:

  Log prior density of the threshold parameters.

- logL:

  Approximated log-likelihood.

- AIC, BIC, ICOMP:

  Akaike information criterion (AIC), Bayesian information criterion
  (BIC), and information complexity index, respectively.

- cov_mat:

  Covariance matrix of the latent variables.

- f_mean:

  Means of the latent variables.

- eff_par:

  Effective number of parameters.

- Options:

  Input information used to call the function.

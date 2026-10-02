# Run an exploratory factor analysis (EFA) using the discretized response model (DRM)

Run an exploratory factor analysis (EFA) using the discretized response
model (DRM)

## Usage

``` r
efa_drm(
  data,
  dimension,
  t_prior = NULL,
  range = c(-4, 4),
  q = 41,
  ngrid = 1000,
  max_iter = NULL,
  threshold = NULL,
  eq_interval = FALSE,
  estimation = "EM"
)
```

## Arguments

- data:

  A matrix or data frame of item responses.

- dimension:

  An integer specifying the number of dimensions of the latent
  variables.

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

- estimation:

  A character string specifying the estimation method: `EM` for the
  expectation-maximization algorithm or `MHRM` for the
  Metropolis-Hastings Robbins-Monro algorithm.

## Value

A list containing the following objects:

- par_est:

  Item parameter estimates.

- iter:

  Number of algorithm iterations.

- diff:

  Value of the convergence criterion at the final iteration.

- logL:

  Approximated log-likelihood.

- f_means:

  Means of the latent factors.

- cov_mat:

  Covariance matrix of the latent variables.

- Options:

  Input information used to call the function.

When `estimation = "EM"`, the returned list also contains:

- se:

  Standard errors of the parameter estimates.

- fk:

  Expected frequency (number of observations) at each grid point.

- quad:

  Locations of the quadrature grid points.

- prior:

  Prior density at each grid point.

- posterior:

  Posterior density at each grid point.

- Ak:

  Posterior weight at each grid point.

- theta:

  Expected a posteriori (EAP) estimate for each observation.

- theta_se:

  Standard errors of the EAP estimates.

- l_prior_w:

  Log prior density of the threshold parameters.

- AIC, BIC, ICOMP:

  Akaike information criterion (AIC), Bayesian information criterion
  (BIC), and information complexity index, respectively.

- eff_par:

  Effective number of parameters.

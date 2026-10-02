# Run an exploratory factor analysis (EFA) for bounded-continuous responses

Run an exploratory factor analysis (EFA) for bounded-continuous
responses

## Usage

``` r
efa_cont(
  data,
  dimension,
  range = c(-4, 4),
  q = 41,
  max_iter = 200,
  threshold = 1e-04
)
```

## Arguments

- data:

  A matrix or data frame of continuous item responses.

- dimension:

  An integer specifying the number of dimensions of the latent
  variables.

- range:

  A numeric vector specifying the range of latent trait values for grid
  approximation in each dimension. Defaults to `c(-4, 4)`.

- q:

  An integer specifying the number of grid points per latent trait
  dimension for grid approximation. Defaults to `41`.

- max_iter:

  An integer specifying the maximum number of iterations for the
  expectation-maximization (EM) algorithm. Defaults to `200`.

- threshold:

  A numeric value specifying the convergence threshold on the maximum
  parameter change. Defaults to `0.0001`.

## Value

A list containing the following objects:

- par_est:

  Estimated item parameters.

- se:

  Standard errors of the item parameter estimates.

- fk:

  Expected frequency (number of respondents) at each quadrature grid
  point.

- iter:

  Number of EM algorithm iterations.

- quad:

  Locations of the quadrature grid points.

- diff:

  Value of the convergence criterion at the final iteration.

- prior:

  Prior density at each quadrature grid point.

- posterior:

  Posterior density at each quadrature grid point per respondent.

- Ak:

  Posterior weight at each quadrature grid point.

- theta:

  Expected a posteriori (EAP) estimates of the latent traits for each
  respondent.

- theta_se:

  Standard errors of the EAP latent trait estimates.

- logL:

  Approximated log-likelihood at the final iteration.

- f_means:

  Estimated means of the latent factors.

- cov_mat:

  Estimated covariance matrix of the latent factors.

- Options:

  Input information used to call the function.

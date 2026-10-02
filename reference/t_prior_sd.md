# Standard Deviations of the Threshold Priors

Standard Deviations of the Threshold Priors

## Usage

``` r
t_prior_sd(N, category, n = 1)
```

## Arguments

- N:

  Sample size.

- category:

  The number of response options.

- n:

  The number of items that share the same prior via constraints.
  Defaults to `1`.

## Value

A vector resulting from the `1/sqrt(N/category*n)*sqrt(2*pi^2/3)`
operation.

## Examples

``` r
# \donttest{
t_prior_sd(
  N = 100,
  category = 5,
  n = 1
)
#> [1] 0.5735737
# }
```

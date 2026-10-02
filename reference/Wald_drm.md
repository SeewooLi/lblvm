# Perform a Wald Test of Parameter Equality

Perform a Wald Test of Parameter Equality

## Usage

``` r
Wald_drm(constraints, fit)
```

## Arguments

- constraints:

  A character string specifying equality constraints between model
  parameters.

- fit:

  A fitted model object returned by a compatible `lblvm` estimation
  function.

## Value

A list containing the following objects:

- W:

  Wald test statistic.

- df:

  Degrees of freedom for the Wald test.

- p_value:

  P-value associated with the Wald test statistic.

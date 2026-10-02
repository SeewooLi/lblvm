# Calculate Local RMSEA

Calculate Local RMSEA

## Usage

``` r
local_rmsea(constraints, fit)
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

- rmsea:

  The local RMSEA estimate based on the Wald statistic.

- ci90:

  A 90\\ `NULL` when the estimated noncentrality parameter is
  non-positive.

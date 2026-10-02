# Plot Response Category Intervals

Creates a horizontal stacked bar plot showing the relative widths of
response category intervals defined by a vector of thresholds.

## Usage

``` r
plot_DRM(thresholds, item_name = NULL, legend = FALSE)
```

## Arguments

- thresholds:

  A numeric vector of threshold parameters.

- item_name:

  An optional character string specifying the label for the item.
  Defaults to `NULL`.

- legend:

  A logical value indicating whether to display the legend. Defaults to
  `FALSE`.

## Value

A `ggplot` object displaying the relative widths of the response
category intervals.

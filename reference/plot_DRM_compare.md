# Compare Response Category Intervals

Creates a horizontal stacked bar plot comparing the relative widths of
response category intervals across multiple sets of threshold
parameters.

## Usage

``` r
plot_DRM_compare(..., label = NULL, item_name = NULL, legend = FALSE)
```

## Arguments

- ...:

  Numeric vectors of threshold parameters to compare.

- label:

  An optional character vector specifying labels for the threshold
  vectors. Defaults to `NULL`, in which case labels are generated as
  `"Group1"`, `"Group2"`, and so on.

- item_name:

  An optional character string specifying the label for the item.
  Defaults to `NULL`.

- legend:

  A logical value indicating whether to display the legend. Defaults to
  `FALSE`.

## Value

A `ggplot` object displaying the relative widths of the response
category intervals for each set of thresholds.

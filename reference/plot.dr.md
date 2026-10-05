# Plot Response Category Intervals for All Items

Extracts the estimated response category thresholds from a fitted DRM
object and creates a horizontal stacked bar plot of the response
category intervals for each item. The resulting item-specific plots are
combined into a single plot using
[`gridExtra::grid.arrange()`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html).

## Usage

``` r
# S3 method for class 'dr'
plot(x, ...)
```

## Arguments

- x:

  A fitted DRM model object containing the estimated response category
  thresholds in `x$par_est[[2]]` and, optionally, item names in
  `x$Options$data`.

- ...:

  Additional arguments passed to
  [`gridExtra::grid.arrange()`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html).
  These can be used to control the arrangement and appearance of the
  combined plots.

## Value

A combined plot displaying the response category intervals for all items
in `x`. Each item is represented by a horizontal stacked bar whose
segments correspond to the relative widths of its response categories.

## See also

[`plot_DRM()`](https://seewooli.github.io/lblvm/reference/plot_DRM.md)
for plotting the response category intervals for a single item.

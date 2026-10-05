#' Plot Response Category Intervals for All Items
#'
#' Extracts the estimated response category thresholds from a fitted DRM
#' object and creates a horizontal stacked bar plot of the response category
#' intervals for each item. The resulting item-specific plots are combined
#' into a single plot using [gridExtra::grid.arrange()].
#'
#' @param x A fitted DRM model object containing the estimated response
#'   category thresholds in `x$par_est[[2]]` and, optionally, item names
#'   in `x$Options$data`.
#' @param ... Additional arguments passed to [gridExtra::grid.arrange()].
#'   These can be used to control the arrangement and appearance of the
#'   combined plots.
#'
#' @returns A combined plot displaying the response category intervals for
#'   all items in `x`. Each item is represented by a horizontal stacked bar
#'   whose segments correspond to the relative widths of its response
#'   categories.
#'
#' @seealso
#' [plot_DRM()] for plotting the response category intervals for a single
#' item.
#'
#' @export
plot.dr <- function(x, ...) {
  item_names <- rownames(x$Options$data)

  thresholds <- apply(
    x$par_est[[2]][, -1],
    1,
    cut_trans,
    simplify = FALSE
  )

  plts <- list()

  for (i in seq_along(thresholds)) {
    plts[[i]] <- plot_DRM(
      thresholds[[i]],
      if (is.null(item_names)) {
        paste0("item", sprintf("%02d", i))
      } else {
        item_names[i]
      }
    )
  }

  gridExtra::grid.arrange(
    grobs = plts,
    ncol = 1,
    ...
  )
}



#' Plot Response Category Intervals
#'
#' Creates a horizontal stacked bar plot showing the relative widths of
#' response category intervals defined by a vector of thresholds.
#'
#' @param thresholds A numeric vector of threshold parameters.
#' @param item_name An optional character string specifying the label for the
#'   item. Defaults to `NULL`.
#' @param legend A logical value indicating whether to display the legend.
#'   Defaults to `FALSE`.
#'
#' @returns A `ggplot` object displaying the relative widths of the response
#'   category intervals.
#'
#' @export
plot_DRM <- function(thresholds, item_name = NULL, legend = FALSE) {
  thresholds <- c(thresholds, 1)
  thresholds <- c(thresholds[1], diff(thresholds))

  data <- data.frame(
    Category = "a",
    Subgroup = paste0("Cat", seq_along(thresholds)),
    Value = rev(thresholds)
  )

  ppp <- ggplot2::ggplot(
    mapping = ggplot2::aes(x = data$Category, y = data$Value, fill = data$Subgroup)
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = "fill",
      show.legend = legend
    ) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(scale = 100),
      breaks = seq(0, 1, by = 1 / length(thresholds))
    ) +
    ggplot2::labs(x = item_name, y = NULL, fill = "Category") +
    ggplot2::theme(
      plot.margin = grid::unit(c(0, 0.1, 0, 0.1), "inches"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(angle = 0),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  ppp
}


#' Compare Response Category Intervals
#'
#' Creates a horizontal stacked bar plot comparing the relative widths of
#' response category intervals across multiple sets of threshold parameters.
#'
#' @param ... Numeric vectors of threshold parameters to compare.
#' @param label An optional character vector specifying labels for the
#'   threshold vectors. Defaults to `NULL`, in which case labels are generated
#'   as `"Group1"`, `"Group2"`, and so on.
#' @param item_name An optional character string specifying the label for the
#'   item. Defaults to `NULL`.
#' @param legend A logical value indicating whether to display the legend.
#'   Defaults to `FALSE`.
#'
#' @returns A `ggplot` object displaying the relative widths of the response
#'   category intervals for each set of thresholds.
#'
#' @export
plot_DRM_compare <- function(
    ...,
    label = NULL,
    item_name = NULL,
    legend = FALSE
) {
  t_list <- list(...)
  n_cat <- length(t_list)

  if (is.null(label)) {
    label <- paste0("Group", seq_len(n_cat))
  }

  if (length(label) != n_cat) {
    stop("Length of `label` must match number of threshold vectors.")
  }

  process_threshold <- function(t) {
    t <- c(t, 1)
    t <- c(t[1], diff(t))
    rev(t)
  }

  t_processed <- lapply(t_list, process_threshold)

  data <- do.call(rbind, lapply(seq_along(t_processed), function(i) {
    data.frame(
      Category = label[i],
      Subgroup = paste0("Cat", seq_along(t_processed[[i]])),
      Value = t_processed[[i]]
    )
  }))

  p <- ggplot2::ggplot(
    mapping = ggplot2::aes(x = data$Category, y = data$Value, fill = data$Subgroup)
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = "fill",
      show.legend = legend
    ) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      breaks = seq(
        0,
        1,
        length.out = length(t_processed[[1]]) + 1
      )
    ) +
    ggplot2::labs(x = item_name, y = NULL, fill = "Subgroup") +
    ggplot2::theme(
      plot.margin = grid::unit(c(0, 0.1, 0, 0.1), "inches"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(angle = 0),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  p
}

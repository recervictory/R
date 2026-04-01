# plotting/stacked_bar.R
# Stacked bar chart of cell metadata composition per cluster.

library(ggplot2)
library(ggtext)
library(dplyr)
library(Seurat)

# Plot a stacked bar chart showing the composition of metadata column y
# broken down by metadata column x (typically seurat_clusters).
# value_type: "percentage" (default) or "count".
# orders: optional factor level order for x.
plot_metadata_stacked_bar <- function(seurat_obj,
                                      x = "seurat_clusters",
                                      y = "sample",
                                      flip_plot = FALSE,
                                      value_type = "percentage",
                                      orders = NULL) {
  if (!y %in% names(seurat_obj@meta.data)) {
    stop("Specified metadata column does not exist in the Seurat object.")
  }

  data <- FetchData(seurat_obj, c(x, y))

  if (nrow(data) == 0) {
    stop("No data was extracted. Check your metadata column name and Seurat object.")
  }

  colnames(data)[1] <- x
  data[[x]] <- as.factor(data[[x]])

  data <- data %>%
    dplyr::count(!!dplyr::sym(x), !!dplyr::sym(y))

  total_counts <- data %>%
    dplyr::group_by(!!dplyr::sym(x)) %>%
    dplyr::summarise(Total = sum(n)) %>%
    dplyr::arrange(dplyr::desc(Total))

  data <- merge(data, total_counts, by = x)

  if (!is.null(orders)) {
    data[[x]] <- factor(data[[x]], levels = orders)
  } else {
    data[[x]] <- factor(data[[x]], levels = total_counts[[x]])
  }

  data <- data %>%
    dplyr::group_by(!!dplyr::sym(x)) %>%
    dplyr::mutate(Percentage = n / sum(n) * 100) %>%
    dplyr::ungroup()

  if (value_type == "count") {
    value_label <- "Count"
    data$Value <- data$n
  } else {
    value_label <- "Percentage"
    data$Value <- data$Percentage
  }

  p <- ggplot(data, aes_string(x = x, y = "Value", fill = y)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    xlab(x) +
    ylab(value_label) +
    labs(fill = y) +
    theme_minimal() +
    theme(axis.text.x = element_markdown(angle = 45, hjust = 1)) +
    scale_x_discrete()

  if (value_type == "count") {
    p <- p +
      ylim(c(0, max(data$Total) * 1.1)) +
      geom_text(
        aes_string(label = "ifelse(Percentage > 30, as.character(n), '')", y = "Value"),
        position = position_stack(vjust = 0.5), size = 3, angle = 90, color = "black"
      ) +
      geom_text(
        aes_string(label = "Total", y = "Total"),
        hjust = -0.1, angle = 90, color = "black"
      )
  } else if (value_type == "percentage") {
    p <- p +
      geom_text(
        aes_string(label = "ifelse(Percentage > 10, paste0(round(Percentage, 1), '%'), '')", y = "Value"),
        position = position_stack(vjust = 0.5), size = 3, angle = 90, color = "black"
      )
  }

  if (flip_plot) {
    p <- p + coord_flip()
  }

  p
}

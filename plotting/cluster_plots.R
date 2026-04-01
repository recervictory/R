# plotting/cluster_plots.R
# Functions for visualizing cluster purity, composition, and UMAP layouts.

library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
library(Seurat)

# Build a summary data frame categorizing each cluster by purity.
# Returns a tibble with columns: cluster, purity, count, category.
purity_df <- function(df, col_order,
                      name = "Cluster",
                      pure = c(0.90, 1.0),
                      mixed = c(0.5, 0.90),
                      take_max = TRUE) {
  df <- as_tibble(df)
  df <- df %>% dplyr::mutate(dplyr::across(dplyr::everything(), ~ round(. * 1, 2)))

  if (take_max) {
    col_values <- df %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), max, na.rm = TRUE, .names = "{col}"))
  } else {
    col_values <- df %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(.) / sum(. != 0), .names = "{col}"))
  }

  non_zero_counts <- df %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(. != 0), .names = "{col}"))

  result <- tibble(
    cluster = as.numeric(names(col_values)),
    purity  = round(as.numeric(col_values), 2),
    count   = as.numeric(non_zero_counts)
  )

  result <- result %>%
    dplyr::mutate(category = dplyr::case_when(
      purity >= pure[1]  & purity <= pure[2]  ~ "Pure",
      purity >= mixed[1] & purity <= mixed[2] ~ "Mixed",
      TRUE ~ "Impure"
    ))

  result$cluster <- factor(result$cluster, levels = col_order)
  return(result)
}


# Compute a normalized (column-proportional) contingency table from two
# Seurat metadata columns.
normalized_table <- function(seuratObject, row_colname = "", col_colname = "") {
  if (!(row_colname %in% colnames(seuratObject@meta.data))) {
    stop(sprintf("Column '%s' not found in metadata.", row_colname))
  }
  if (!(col_colname %in% colnames(seuratObject@meta.data))) {
    stop(sprintf("Column '%s' not found in metadata.", col_colname))
  }

  contingency_table <- table(
    seuratObject@meta.data[[row_colname]],
    seuratObject@meta.data[[col_colname]],
    useNA = "ifany"
  )

  normalized <- apply(contingency_table, 2, function(x) x / sum(x))
  return(normalized)
}


# Combine a UMAP plot, a split UMAP by purity category, and a purity tile bar
# into a single patchwork figure.
cluster_purity_split_plots <- function(seuratObject, cluster_df,
                                       group.by = "seurat_clusters",
                                       split_by = "purity_category",
                                       layout = NA) {
  if (is.na(layout)) {
    layout <- "
    ABBB
    ABBB
    ABBB
    ABBB
    ABBB
    ABBB
    ABBB
    CCCC
    "
  }

  A_seurat_clusters <- DimPlot(seuratObject, reduction = "umap", alpha = 0.5,
                                group.by = "seurat_clusters", label = TRUE) + NoLegend()

  B_seurat_clusters_split_purity <- DimPlot(seuratObject, reduction = "umap",
                                             group.by = group.by, split.by = split_by,
                                             label = TRUE, repel = TRUE)

  C_clusters_plot <- ggplot(cluster_df, aes(x = cluster, y = 1, fill = category)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = count), vjust = 0, size = 5) +
    labs(x = "Cluster", y = "", title = "Cluster's Purity") +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

  patch_image <- A_seurat_clusters + B_seurat_clusters_split_purity + C_clusters_plot +
    plot_layout(design = layout)

  return(patch_image)
}


# Combined metadata visualization: UMAP plots + stacked bar charts + purity tile.
# Requires plotting/stacked_bar.R to be sourced for plot_metadata_stacked_bar.
clusters_meta_visualization <- function(seuratObject,
                                         x = "seurat_clusters",
                                         y = "sample",
                                         layout = NA,
                                         cluster_table,
                                         col_order) {
  cluster_df <- purity_df(cluster_table, col_order)

  if (is.na(layout)) {
    layout <- "
AABBB
AABBB
AABBB
AABBB
CCDDD
CCDDD
CCDDD
CCDDD
EEEEE
"
  }

  clusters_plot <- ggplot(cluster_df, aes(x = cluster, y = 1, fill = category)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = count), vjust = 0, size = 5) +
    labs(x = "cluster", y = "", title = "Cluster's Purity") +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

  plot_count   <- plot_metadata_stacked_bar(seuratObject, x = x, y = y, value_type = "count",   orders = col_order)
  plot_percent <- plot_metadata_stacked_bar(seuratObject, x = x, y = y,                          orders = col_order)

  plot_sample   <- DimPlot(seuratObject, reduction = "umap", group.by = y, label = TRUE) +
    theme(legend.position = "none")
  plot_clusters <- DimPlot(seuratObject, reduction = "umap", group.by = x, label = TRUE) + NoLegend()

  patch_image <- plot_clusters + plot_count + plot_sample + plot_percent + clusters_plot +
    plot_layout(design = layout, guides = "collect")

  return(patch_image)
}

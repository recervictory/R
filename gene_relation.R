plot_gene_trends <- function(seurat_obj, gene_list, sort_by_gene = NULL, bin_size = 50, 
                             colors = NULL, cluster_colors = NULL, 
                             annotation_name = "seurat_clusters",
                             sort_by_gene_cutoff = 0,
                             nfeatures = 100,
                             min_range = 0.5,
                             sort_by_annotation = TRUE,
                             log_transform = FALSE,
                             verbose = TRUE) {
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  
  if (verbose) message("[Step 1] Preparing sort_by_gene...")
  if (is.null(sort_by_gene)) {
    hvg <- VariableFeatures(seurat_obj)
    assay_data <- GetAssayData(seurat_obj, slot = "data")[hvg, , drop = FALSE]
    gene_ranges <- apply(assay_data, 1, function(x) diff(range(x)))
    filtered_hvg <- names(sort(gene_ranges, decreasing = TRUE))[gene_ranges > min_range]
    sort_by_gene <- head(filtered_hvg, nfeatures)
    if (verbose) message("Using top ", length(sort_by_gene), " HVGs (range-filtered) for sorting.")
  }
  
  all_genes <- unique(c(gene_list, sort_by_gene))
  if (verbose) message("[Step 2] Fetching expression for ", length(all_genes), " genes...")
  expr_df <- FetchData(seurat_obj, vars = all_genes)
  
  missing_genes <- setdiff(sort_by_gene, colnames(expr_df))
  if (length(missing_genes) > 0) {
    stop(paste0("Missing sort_by_gene(s): ", paste(missing_genes, collapse = ", ")))
  }

  if (verbose) message("[Step 3] Calculating average expression for sort_by_gene...")
  expr_df <- expr_df %>%
    mutate(sort_by_gene_avg = rowMeans(across(all_of(sort_by_gene)), na.rm = TRUE))
  
  if (log_transform) {
    if (verbose) message("[Step 4] Applying log1p transformation to gene_list...")
    expr_df <- expr_df %>%
      mutate(across(all_of(gene_list), log1p))
  }
  
  if (verbose) message("[Step 5] Fetching cluster info from ", annotation_name)
  cluster_df <- FetchData(seurat_obj, vars = annotation_name)
  colnames(cluster_df) <- "cluster"
  
  expr_df <- expr_df %>%
    mutate(cluster = cluster_df$cluster) %>%
    filter(sort_by_gene_avg >= sort_by_gene_cutoff)
  
  if (verbose) message("[Step 6] Sorting cells...")
  expr_df <- if (sort_by_annotation) {
    expr_df %>%
      mutate(cluster = as.factor(cluster)) %>%
      arrange(cluster, sort_by_gene_avg)
  } else {
    expr_df %>% arrange(sort_by_gene_avg)
  }

  expr_df <- expr_df %>%
    mutate(cell_order = row_number(),
           bin = ntile(sort_by_gene_avg, bin_size))
  
  if (verbose) message("[Step 7] Calculating binned means...")
  binned_means <- expr_df %>%
    group_by(bin) %>%
    summarise(across(all_of(gene_list), mean), .groups = "drop") %>%
    mutate(cell_order = bin * (max(expr_df$cell_order) / bin_size))
  
  plot_df <- binned_means %>%
    pivot_longer(cols = all_of(gene_list), names_to = "gene", values_to = "expression")

  if (is.null(colors)) {
    colors <- setNames(RColorBrewer::brewer.pal(length(gene_list), "Set1"), gene_list)
  }
  
  sort_label <- if (is.null(sort_by_gene)) {
    paste0("Top ", nfeatures, " HVGs (range > ", min_range, ")")
  } else {
    paste(sort_by_gene, collapse = " + ")
  }
  x_label <- if (sort_by_annotation) {
    paste("Cells grouped by", annotation_name, "and sorted by", sort_label)
  } else {
    paste("Cells sorted by", sort_label)
  }

  if (verbose) message("[Step 8] Building expression trend plot...")
  main_plot <- ggplot(plot_df, aes(x = cell_order, y = expression, color = gene)) +
    geom_line(size = 1.2) +
    geom_text(data = plot_df %>% group_by(gene) %>% slice_tail(n = 1),
              aes(label = gene), hjust = -0.1, vjust = 0.5, show.legend = FALSE, fontface = "bold") +
    scale_color_manual(values = colors) +
    labs(x = x_label, y = "Average expression", title = "Smoothed gene expression trends") +
    xlim(NA, max(plot_df$cell_order) * 1.1) +
    theme_minimal()
  
  if (verbose) message("[Step 9] Preparing annotation bar...")
  annotation_df <- expr_df %>% dplyr::select(cell_order, cluster)
  
  if (is.null(cluster_colors)) {
    unique_clusters <- sort(unique(annotation_df$cluster))
    n_clusters <- length(unique_clusters)
    cluster_palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_clusters)
    cluster_colors <- setNames(cluster_palette, unique_clusters)
  }
  
  annotation_plot <- ggplot(annotation_df, aes(x = cell_order, y = 1, fill = factor(cluster))) +
    geom_tile() +
    scale_fill_manual(values = cluster_colors, name = annotation_name) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(0, 10, 0, 10)
    )

  if (verbose) message("[Step 10] Combining plots...")
  combined_plot <- main_plot / annotation_plot + plot_layout(heights = c(4, 0.6))
  
  if (verbose) message("✅ Done.")
  return(combined_plot)
}

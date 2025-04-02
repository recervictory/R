plot_gene_trends <- function(seurat_obj, gene_list, sort_by_gene, bin_size = 50, 
                             colors = NULL, cluster_colors = NULL, 
                             annotation_name = "seurat_clusters",
                             sort_by_gene_cutoff = 0) {
  
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  
  # Step 1: Combine all gene names needed (expression + sorting)
  all_genes <- unique(c(gene_list, sort_by_gene))
  
  # Step 2: Fetch gene expression
  expr_df <- FetchData(seurat_obj, vars = all_genes)
  
  # Step 3: Validate all sort_by_gene genes are present
  missing_genes <- setdiff(sort_by_gene, colnames(expr_df))
  if (length(missing_genes) > 0) {
    stop(paste0("The following sort_by_gene(s) not found in Seurat object: ", paste(missing_genes, collapse = ", ")))
  }

  # Step 4: Compute average expression of sort_by_gene(s)
  expr_df <- expr_df %>%
    mutate(sort_by_gene_avg = rowMeans(across(all_of(sort_by_gene)), na.rm = TRUE))
  
  # Step 5: Fetch cluster info
  cluster_df <- FetchData(seurat_obj, vars = annotation_name)
  colnames(cluster_df) <- "cluster"
  
  # Step 6: Combine annotation
  expr_df <- expr_df %>%
    mutate(cluster = cluster_df$cluster)
  
  # Step 7: Filter cells below cutoff
  expr_df <- expr_df %>%
    filter(sort_by_gene_avg >= sort_by_gene_cutoff)
  
  # Step 8: Sort and bin
  expr_df <- expr_df %>%
    arrange(sort_by_gene_avg) %>%
    mutate(cell_order = row_number(),
           bin = ntile(sort_by_gene_avg, bin_size))
  
  # Step 9: Calculate average expression per bin
  binned_means <- expr_df %>%
    group_by(bin) %>%
    summarise(across(all_of(gene_list), mean)) %>%
    mutate(cell_order = bin * (max(expr_df$cell_order) / bin_size))
  
  # Step 10: Convert to long format for plotting
  plot_df <- binned_means %>%
    pivot_longer(cols = all_of(gene_list), names_to = "gene", values_to = "expression")
  
  # Step 11: Gene colors
  if (is.null(colors)) {
    colors <- setNames(RColorBrewer::brewer.pal(n = length(gene_list), name = "Set1"), gene_list)
  }

  # Step 12: Main expression plot
  main_plot <- ggplot(plot_df, aes(x = cell_order, y = expression, color = gene)) +
    geom_line(size = 1.2) +
    geom_text(data = plot_df %>% group_by(gene) %>% slice_tail(n = 1),
              aes(label = gene), hjust = -0.1, vjust = 0.5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = colors) +
    labs(x = paste("Cells sorted by", paste(sort_by_gene, collapse = " + "), "avg expression"),
         y = "Average expression",
         title = "Smoothed gene expression trends") +
    xlim(NA, max(plot_df$cell_order) * 1.1) +
    theme_minimal()

  # Step 13: Cluster annotation
  annotation_df <- expr_df %>% select(cell_order, cluster)
  
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

  # Step 14: Combine plots
  combined_plot <- main_plot / annotation_plot + plot_layout(heights = c(4, 0.6))
  
  return(combined_plot)
}

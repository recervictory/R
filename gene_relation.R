library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)



library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

plot_gene_trends <- function(seurat_obj, gene_list, sort_by_gene, bin_size = 50, colors = NULL, cluster_colors = NULL, annotation_name = "seurat_clusters") {

  
  # Step 1: Fetch gene expression
  expr_df <- FetchData(seurat_obj, vars = gene_list)
  
  # Step 2: Check sort gene is present
  if (!(sort_by_gene %in% colnames(expr_df))) {
    stop(paste0("Gene '", sort_by_gene, "' not found in Seurat object"))
  }
  
  # Step 3: Fetch cluster info
  cluster_df <- FetchData(seurat_obj, vars = annotation_name)
  colnames(cluster_df) <- "cluster"
  
  # Step 4: Combine and sort by sort_by_gene
  expr_df <- expr_df %>%
    mutate(cluster = cluster_df$cluster) %>%
    arrange(.data[[sort_by_gene]]) %>%
    mutate(cell_order = row_number(),
           bin = ntile(.data[[sort_by_gene]], bin_size))
  
  # Step 5: Calculate average expression per bin
  binned_means <- expr_df %>%
    group_by(bin) %>%
    summarise(across(all_of(gene_list), mean)) %>%
    mutate(cell_order = bin * (max(expr_df$cell_order) / bin_size))
  
  # Step 6: Convert to long format for ggplot
  plot_df <- binned_means %>%
    pivot_longer(cols = all_of(gene_list), names_to = "gene", values_to = "expression")
  
  # Step 7: Define gene line colors
  if (is.null(colors)) {
    colors <- setNames(RColorBrewer::brewer.pal(n = length(gene_list), name = "Set1"), gene_list)
  }

  # Step 8: Create main line plot
  main_plot <- ggplot(plot_df, aes(x = cell_order, y = expression, color = gene)) +
    geom_line(size = 1.2) +
    geom_text(data = plot_df %>% group_by(gene) %>% slice_tail(n = 1),
              aes(label = gene), hjust = -0.1, vjust = 0.5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = colors) +
    labs(x = paste("Cells sorted by", sort_by_gene, "expression"),
         y = "Average expression", 
         title = "Smoothed gene expression trends") +
    xlim(NA, max(plot_df$cell_order) * 1.1) +
    theme_minimal()

  # Step 9: Create annotation plot
  annotation_df <- expr_df %>% select(cell_order, cluster)
  
  # Default cluster colors if not provided
  if (is.null(cluster_colors)) {
  unique_clusters <- sort(unique(annotation_df$cluster))
  n_clusters <- length(unique_clusters)
  
  # Generate enough distinct colors
  cluster_palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_clusters)
  cluster_colors <- setNames(cluster_palette, unique_clusters)
}

  annotation_plot <- ggplot(annotation_df, aes(x = cell_order, y = 1, fill = factor(cluster))) +
    geom_tile() +
    scale_fill_manual(values = cluster_colors) +
    theme_void() +
          theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(0, 10, 0, 10)
  )

  # Step 10: Combine both plots
  combined_plot <- main_plot / annotation_plot + plot_layout(heights = c(4, 0.4))
  
  return(combined_plot)
}

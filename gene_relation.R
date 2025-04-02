library(Seurat)
library(ggplot2)
library(dplyr)




plot_gene_trends <- function(seurat_obj, gene_list, sort_by_gene, bin_size = 50, colors = NULL) {
  # Step 1: Fetch gene expression
  expr_df <- FetchData(seurat_obj, vars = gene_list)
  
  # Step 2: Check sort gene is present
  if (!(sort_by_gene %in% colnames(expr_df))) {
    stop(paste0("Gene '", sort_by_gene, "' not found in Seurat object"))
  }
  
  # Step 3: Sort by selected gene and bin
  expr_df <- expr_df %>%
    arrange(.data[[sort_by_gene]]) %>%
    mutate(cell_order = row_number(),
           bin = ntile(.data[[sort_by_gene]], bin_size))
  
  # Step 4: Calculate average expression for each gene in each bin
  binned_means <- expr_df %>%
    group_by(bin) %>%
    summarise(across(all_of(gene_list), mean)) %>%
    mutate(cell_order = bin * (max(expr_df$cell_order) / bin_size))
  
  # Step 5: Reshape to long format for ggplot
  plot_df <- binned_means %>%
    pivot_longer(cols = all_of(gene_list), names_to = "gene", values_to = "expression")
  
  # Step 6: Color mapping
  if (is.null(colors)) {
    colors <- setNames(RColorBrewer::brewer.pal(n = length(gene_list), name = "Set1"), gene_list)
  }

  # Step 7: Plot
  p <- ggplot(plot_df, aes(x = cell_order, y = expression, color = gene)) +
    geom_line(size = 1.2) +
    # Add labels at the end of each gene line
    geom_text(data = plot_df %>% group_by(gene) %>% slice_tail(n = 1),
              aes(label = gene), hjust = -0.1, vjust = 0.5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = colors) +
    labs(x = paste("Cells sorted by", sort_by_gene, "expression"),
         y = "Average expression", 
         title = "Smoothed expression trends across sorted cells") +
    xlim(NA, max(plot_df$cell_order) * 1.1) +  # Extra space for labels
    theme_minimal()

  return(p)
}

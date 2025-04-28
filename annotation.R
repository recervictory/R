#### Annotation ####

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R") 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")



save_plot <- function(plot_object, plot_dir, count = "00", file_name = "image", project_name = "", format = "png", batch_name = "", plot_type = "scatter", width = 800, height = 800, dpi = 72, ...) {
  
  # Ensure plot_dir exists
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  # Normalize file path and construct file name
  file_extension <- switch(format,
                           png = "png",
                           pdf = "pdf",
                           jpeg = "jpg",
                           tiff = "tiff",
                           "png") # default to png if format is not recognized
  
  file_name <- paste(count, project_name, batch_name, file_name, plot_type, sep = "_")
  file_name <- paste0(file_name, ".", file_extension)
  file_path <- file.path(plot_dir, file_name)
  file_path <- normalizePath(file_path, mustWork = FALSE)
  
  # Save the plot in the specified format
  switch(format,
         png = {
           png(file_path, width = width, height = height, res = dpi)
           print(plot_object)
           dev.off()
         },
         pdf = {
           width_in <- width / dpi
           height_in <- height / dpi
           pdf(file_path, width = width_in, height = height_in)
           print(plot_object)
           dev.off()
         },
         jpeg = {
           jpeg(file_path, width = width, height = height, res = dpi)
           print(plot_object)
           dev.off()
         },
         tiff = {
           tiff(file_path, width = width, height = height, res = dpi)
           print(plot_object)
           dev.off()
         },
         {
           # Default to PNG if format is not recognized
           png(file_path, width = width, height = height, res = dpi)
           print(plot_object)
           dev.off()
         })
  
  # Log file saved message
  if (exists("flog.info")) {
    flog.info("%s File Saved to %s", toupper(format), file_path)
  } else {
    message(toupper(format), " File Saved to ", file_path)
  }
}



### ClusterTree Object plots
run_clustree <- function(seuratObject, resolutions = c(0.1, 0.25, 0.5, 1.0, 1.5, 2), log_level = "INFO") {
  library(futile.logger)
  library(clustree)
  
  # Set logging level
  flog.threshold(log_level)
  flog.info("😎 Function Name: run_clustree")
  
  # Log the input resolutions
  flog.info("🔧 Resolutions provided: %s", paste(resolutions, collapse = ", "))
  
  # Initialize the Seurat object for clustering
  clustree_seurat_object <- seuratObject
  
  # Perform clustering for each resolution
  for (res in resolutions) {
    flog.info("🔍 Starting FindClusters with resolution: %0.2f", res)
    
    # Run clustering at the current resolution
    clustree_seurat_object <- FindClusters(clustree_seurat_object, graph.name = "RNA_snn", resolution = res, algorithm = 1, verbose = FALSE)
    
    # Log completion of clustering for the current resolution
    flog.info("✅ FindClusters completed successfully for resolution: %0.2f", res)
  }
  
  # Log the start of clustree plot generation
  flog.info("📊 Generating clustree plot for resolutions: %s", paste(resolutions, collapse = ", "))
  
  # Generate the clustree plot
  plot <- clustree(clustree_seurat_object@meta.data, prefix = "RNA_snn_res.")
  
  flog.info("🎉 Clustree plot generated successfully")
  
  # Return the plot object
  return(plot)
}


#### Noramlized TABLE #####
normalized_table <- function(seuratObject, row_colname ="", col_colname = "") {
  # Check if column names exist in the metadata
  if (!(row_colname %in% colnames(seuratObject@meta.data))) {
    stop(sprintf("Column '%s' not found in metadata.", row_colname))
  }
  if (!(col_colname %in% colnames(seuratObject@meta.data))) {
    stop(sprintf("Column '%s' not found in metadata.", col_colname))
  }
  
  # Create contingency table
  contingency_table <- table(
    seuratObject@meta.data[[row_colname]],
    seuratObject@meta.data[[col_colname]],
    useNA = "ifany"
  )
  
  # Normalize by the total number of cells in each column (clusters or other categories)
  normalized_table <- apply(
    contingency_table, 
    2, # Apply function to columns
    function(x) x / sum(x) # Normalize each column
  )
  
  return(normalized_table)
}


#### Find Purity of Clusters
purity_df <- function(df, col_order, name = "Cluster", pure = c(0.90, 1.0), mixed = c(0.5, 0.90), take_max = TRUE) {
  # Convert the dataframe to a tibble for easier manipulation
  df <- as_tibble(df)
  
  # Change values to percentage and round to 2 decimal points
  df <- df %>% mutate(across(everything(), ~ round(. * 1, 2)))
  
  # Calculate the average column purity or maximum column purity
  if (take_max) {
    col_values <- df %>%
      summarise(across(everything(), max, na.rm = TRUE, .names = "{col}"))
  } else {
    col_values <- df %>%
      summarise(across(everything(), ~ sum(.) / sum(. != 0), .names = "{col}"))
  }
  
  non_zero_counts <- df %>%
    summarise(across(everything(), ~ sum(. != 0), .names = "{col}"))
  
  # Create the result data.frame
  result <- tibble(
    cluster = as.numeric(names(col_values)),
    purity = round(as.numeric(col_values), 2),
    count = as.numeric(non_zero_counts)
  )
  
  # Categorize based on purity thresholds
  result <- result %>%
    mutate(category = case_when(
      purity >= pure[1] & purity <= pure[2] ~ "Pure",
      purity >= mixed[1] & purity <= mixed[2] ~ "Mixed",
      TRUE ~ "Impure"
    ))
  
  # Order clusters based on the input order
  result$cluster <- factor(result$cluster, levels = col_order)
  
  return(result)
}




######## Cluster Purity ########
cluster_purity_split_plots <- function(seuratObject, cluster_df,group.by = "seurat_clusters", split_by = "purity_category", layout = NA) {
  # Check if layout is NA and assign default layout if true
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
  
  
  # Generate Seurat cluster UMAP plot
  A_seurat_clusters <- DimPlot(seuratObject, reduction = "umap", alpha = 0.5, group.by = "seurat_clusters", label = TRUE) + NoLegend()
  
  # Generate split UMAP plot by purity category
  B_seurat_clusters_split_purity <- DimPlot(seuratObject, reduction = "umap", group.by = group.by, split.by = split_by, label = TRUE,repel = TRUE) 
  
  # Generate clusters plot
  C_clusters_plot <- ggplot(cluster_df, aes(x = cluster, y = 1, fill = category)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = count), vjust = 0, size = 5) +
    labs(x = "Cluster", y = "", title = "Cluster's Purity") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank())
  
  # Combine plots using the specified layout
  patch_image <- A_seurat_clusters +
    B_seurat_clusters_split_purity +
    C_clusters_plot +
    plot_layout(design = layout)
  
  return(patch_image)
}

####### Create Annotation ######
create_annotation <- function(seuratObject, col_names, by = "seurat_clusters") {
  
  # Extract meta.data from Seurat object
  meta.data <- seuratObject@meta.data
  
  # Select columns of interest from meta.data
  annotation <- meta.data[, col_names]
  
  # Group by 'seurat_clusters' and summarize counts
  annotation <- annotation %>%
    group_by(across(all_of(col_names))) %>%
    summarize(count = n(), .groups = 'drop') %>%
    group_by(across(all_of(by))) %>%
    slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Prepare final annotation dataframe
  annotation <- as.data.frame(annotation)
  
  # Remove unnecessary columns and format
  annotation[[by]] <- as.character(annotation[[by]])
  row.names(annotation) <- annotation[[by]]
  annotation[[by]] <- NULL
  annotation$count <- NULL
  
  return(annotation)
} 


plot_metadata_stacked_bar <- function(seurat_obj, x = "seurat_clusters", y = "sample", flip_plot = FALSE, value_type = "percentage", orders = NULL) {
  # Ensure the metadata column exists
  if (!y %in% names(seurat_obj@meta.data)) {
    stop("Specified metadata column does not exist in the Seurat object.")
  }
  
  # Extract data
  data <- FetchData(seurat_obj, c(x, y))
  head(data)
  
  # Check if data extraction was successful
  if (nrow(data) == 0) {
    stop("No data was extracted. Check your metadata column name and Seurat object.")
  }
  
  # Renaming for clarity, using the value of `x` for naming
  colnames(data)[1] <- x
  data[[x]] <- as.factor(data[[x]])
  
  
  # Calculating counts
  data <- data %>%
    dplyr::count(!!dplyr::sym(x), !!dplyr::sym(y))
  
  # Calculate total counts for sorting and for annotation
  total_counts <- data %>%
    dplyr::group_by(!!dplyr::sym(x)) %>%
    dplyr::summarise(Total = sum(n)) %>%
    dplyr::arrange(desc(Total))
  
  # Merge total counts back into data
  data <- merge(data, total_counts, by = x)
  
  # Order based on total count
  if (!is.null(orders)) {
    data[[x]] <- factor(data[[x]], levels = orders)
  } else {
    # Order based on total count if no order is provided
    data[[x]] <- factor(data[[x]], levels = total_counts[[x]])
  }
  
  # Calculate percentages for all cases because it's needed for conditional labeling
  data <- data %>%
    dplyr::group_by(!!dplyr::sym(x)) %>%
    dplyr::mutate(Percentage = n / sum(n) * 100) %>%
    dplyr::ungroup()
  
  # Decide the value and label based on the selected value_type
  if (value_type == "count") {
    value_label <- "Count"
    data$Value <- data$n
  } else {  # default to percentage
    value_label <- "Percentage"
    data$Value <- data$Percentage
  }
  
  # Generate default ggplot2 discrete colors
  x_levels <- levels(data[[x]])
  num_colors <- length(x_levels)
  
  
  
  # Plot
  p <- ggplot(data, aes_string(x = x, y = "Value", fill = y)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    xlab(x) +
    ylab(value_label) +
    labs(fill = y) +
    theme_minimal() +
    theme(
      axis.text.x = element_markdown(angle = 45, hjust = 1)
    ) +
    scale_x_discrete()
  
  # Conditional labels based on value_type and percentage threshold
  if (value_type == "count") {
    p <- p +
      ylim(c(0, max(data$Total) * 1.1)) +
      geom_text(aes_string(label = "ifelse(Percentage > 30, as.character(n), '')", y = "Value"), 
                position = position_stack(vjust = 0.5), 
                size = 3, 
                angle = 90, 
                color = "black") +
      geom_text(aes_string(label = "Total", y = "Total"), 
                hjust = -0.1,angle = 90, color = "black") 
  } else if (value_type == "percentage") {
    p <- p +
      geom_text(aes_string(label = "ifelse(Percentage > 10, paste0(round(Percentage, 1), '%'), '')", y = "Value"), 
                position = position_stack(vjust = 0.5), 
                size = 3, 
                angle = 90, 
                color = "black")
  }
  
  if (flip_plot) {
    p <- p + coord_flip()
  }
  
  p
}


###### Cluster Meta Vis #######
clusters_meta_visualization  <- function(seuratObject, x = "seurat_clusters", y = "sample", layout = NA ,cluster_table, col_order) {
  
  
  
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
    labs(x = "cluster", y = "", title = "Cluters's Purity") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank())
  
  # Generate the plots for metadata stacked bars
  plot_count <- plot_metadata_stacked_bar(seuratObject, x = x, y = y, value_type = "count", orders = col_order) 
  plot_percent <- plot_metadata_stacked_bar(seuratObject, x = x, y = y, orders = col_order)
  
  # Generate the UMAP plots
  plot_sample <- DimPlot(seuratObject, reduction = "umap", group.by = y, label = TRUE) + theme(legend.position = "none")
  plot_clusters <- DimPlot(seuratObject, reduction = "umap", group.by = x, label = TRUE) + NoLegend()
  
  # Combine the plots
  patch_image <- plot_clusters + plot_count + plot_sample + plot_percent + clusters_plot + plot_layout(design = layout, guides = "collect")
  
  # Return the combined plot
  return(patch_image)
}









########## Replace NA #####
replace_surat_NA <- function(seurat_obj, col_name, replacement_value = "Unknown") {
  # Check if the specified column exists in the metadata
  if (!(col_name %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", col_name, "does not exist in the Seurat object metadata."))
  }
  
  na_count <- sum(is.na(seurat_obj@meta.data[[col_name]]))
  
  flog.info("NA Found : %s",na_count)
  
  # Replace NA values in the specified column with the replacement value
  seurat_obj@meta.data[[col_name]][is.na(seurat_obj@meta.data[[col_name]])] <- replacement_value
  
  # Return the modified Seurat object
  return(seurat_obj)
}
  
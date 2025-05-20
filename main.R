#### Import Library #####
library(futile.logger)
library(Cairo)
flog.info("🤖 Function File: main.R")


#### Description save_plot  Function######
# The save_plot function is designed to save a given plot object to a specified directory in various formats such as PNG, PDF, JPEG, and TIFF.
# It handles the creation of the directory if it does not exist and constructs a standardized file name based on various parameters like project
# name, batch, plot type, and format. The function also allows specifying image dimensions and resolution.


save_plot <- function(plot_object, plot_dir,
                      count = "00",
                      file_name = "image",
                      project_name = "",
                      format = "png",
                      batch_name = "",
                      plot_type = "scatter",
                      width = 800, height = 800,
                      dpi = 72, ...) {
  flog.info("😎 Function Name: save_plot")

  # Ensure plot_dir exists
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  # Determine file extension
  file_extension <- switch(format,
    png = "png",
    pdf = "pdf",
    jpeg = "jpg",
    jpg = "jpg",
    tiff = "tiff",
    "png"
  )

  # Construct the file name and path
  file_name <- paste(count, project_name, batch_name, file_name, plot_type, sep = "_")
  file_name <- paste0(file_name, ".", file_extension)
  file_path <- file.path(plot_dir, file_name)
  file_path <- normalizePath(file_path, mustWork = FALSE)

  # Save the plot based on format
  tryCatch({
    if (format == "png") {
      if (capabilities("cairo")) {
        Cairo::CairoPNG(file_path, width = width, height = height, res = dpi)
      } else {
        png(file_path, width = width, height = height, res = dpi)
      }
      print(plot_object)
      dev.off()
    } else if (format %in% c("jpeg", "jpg")) {
      jpeg(file_path, width = width, height = height, res = dpi)
      print(plot_object)
      dev.off()
    } else if (format == "tiff") {
      tiff(file_path, width = width, height = height, res = dpi)
      print(plot_object)
      dev.off()
    } else if (format == "pdf") {
      width_in <- width / dpi
      height_in <- height / dpi
      pdf(file_path, width = width_in, height = height_in)
      print(plot_object)
      dev.off()
    } else {
      # Default to PNG
      png(file_path, width = width, height = height, res = dpi)
      print(plot_object)
      dev.off()
    }

    # Log success
    if (exists("flog.info")) {
      flog.info("%s file saved to %s", toupper(format), file_path)
    } else {
      message(toupper(format), " file saved to ", file_path)
    }
  }, error = function(e) {
    flog.error("Error saving plot: %s", e$message)
    stop(e)
  })
}


# Seurat Analysis Workflow
runSeuratWorkflow <- function(object) {
  flog.info("😎 Function Name: runSeuratWorkflow - Loaction ")

  if (!"umap" %in% names(Reductions(object))) {
    flog.info("UMAP not found, starting Seurat Analysis Workflow.")
    object <- NormalizeData(object = object)
    object <- FindVariableFeatures(object = object)
    object <- ScaleData(object = object)
    object <- RunPCA(object = object, verbose = FALSE)
    object <- FindNeighbors(object = object, dims = 1:30)
    object <- FindClusters(object = object, resolution = 0.5)
    object <- RunUMAP(object = object, dims = 1:30, verbose = FALSE)
  } else {
    flog.info("UMAP already exists, skipping workflow.")
  }
  flog.info("Completed Seurat Analysis Workflow.")
  return(object)
}



filter_qc <- function(seurat_obj, min_nFeature_RNA = 200, max_nFeature_RNA = 9000,
                      min_nCount_RNA = 1000, max_nCount_RNA = 45000, mito_high_cutoff = 25, ...) {
  # Add the QC_Filtered column with default value "Passed"
  seurat_obj@meta.data$QC_Consensus_Filtered <- "Passed"

  # Apply the conditions to update QC_Filtered column
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$nCount_RNA <= min_nCount_RNA] <- "Low nCOUNT"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$nFeature_RNA <= min_nFeature_RNA] <- "Low Genes per Cell"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$QCs_miQC_keep == FALSE] <- "MiQC+"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$nFeature_RNA >= max_nFeature_RNA] <- "High Genes per Cell"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$nCount_RNA >= max_nCount_RNA] <- "High nCOUNT"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$QCs_doubletFinder == "doublet"] <- "DoubletFinder+"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$QCs_scDblFinder_class == "doublet"] <- "DblFinder+"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$QCs_doubletFinder == "doublet" &
    seurat_obj@meta.data$QCs_scDblFinder_class == "doublet"] <- "Doublet"
  return(seurat_obj)
}


plot_metadata_stacked_bar <- function(seurat_obj, x = "ident", y = "sample", flip_plot = FALSE, value_type = "percentage", orders = NA, palette = NA) {
  # Ensure the metadata column exists
  if (!y %in% names(seurat_obj@meta.data)) {
    stop("Specified metadata column does not exist in the Seurat object.")
  }

  # Extract data
  data <- FetchData(seurat_obj, c(x, y))

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



  # Calculate percentages for all cases because it's needed for conditional labeling
  data <- data %>%
    dplyr::group_by(!!dplyr::sym(x)) %>%
    dplyr::mutate(Percentage = n / sum(n) * 100) %>%
    dplyr::ungroup()

  # Decide the value and label based on the selected value_type
  if (value_type == "count") {
    value_label <- "Count"
    data$Value <- data$n
  } else { # default to percentage
    value_label <- "Percentage"
    data$Value <- data$Percentage
  }

  data[[x]] <- factor(data[[x]], levels = total_counts[[x]])

  # Plot
  p <- ggplot(data, aes_string(x = x, y = "Value", fill = y)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    xlab(x) +
    ylab(value_label) +
    labs(fill = y) +
    scale_fill_manual(values = palette) + # Use manual scale_fill_manual to apply generated palette
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1) # Correct usage of element_text for rotating x-axis labels
    ) +
    scale_x_discrete()

  # Conditional labels based on value_type and percentage threshold
  if (value_type == "count") {
    p <- p +
      geom_text(aes_string(label = "ifelse(Percentage > 30, as.character(n), '')", y = "Value"),
        position = position_stack(vjust = 0.5),
        size = 3,
        angle = 90,
        color = "black"
      ) +
      geom_text(aes_string(label = "Total", y = "Total"),
        hjust = 0.01, size = 4, angle = 90, color = "black"
      )
  } else if (value_type == "percentage") {
    p <- p +
      geom_text(aes_string(label = "ifelse(Percentage > 10, paste0(round(Percentage, 1), '%'), '')", y = "Value"),
        position = position_stack(vjust = 0.5),
        size = 3,
        angle = 90,
        color = "black"
      )
  }

  if (flip_plot) {
    p <- p + coord_flip()
  }

  p
}

read_csv_and_create_variables <- function(file_name) {
  flog.info("😎 Function Name: read_csv_and_create_variables ")
  data <- fread(file_name)
  
  for (i in 1:nrow(data)) {
    assign(data$filter[i], data$value[i], envir = .GlobalEnv)
  }
  
  flog.info("Variables have been recreated from the CSV file '%s'.", file_name)
}


library(dplyr)
library(Seurat)

filterSeurat <- function(seurat_obj,
                         min_nFeature_RNA = 200,
                         max_nFeature_RNA = 9000,
                         min_nCount_RNA = 1000,
                         max_nCount_RNA = 45000,
                         mito_high_cutoff = 25,
                         samples = c(),
                         filter_dir = "../01_metadata/") {
  # Start function log
  flog.info("😎 Function Name: filterSeurat")
  flog.info("🛠️  Parameters - min_nFeature_RNA: %d, max_nFeature_RNA: %d, min_nCount_RNA: %d, max_nCount_RNA: %d, mito_high_cutoff: %d", 
            min_nFeature_RNA, max_nFeature_RNA, min_nCount_RNA, max_nCount_RNA, mito_high_cutoff)
  
  # Ensure the directory for filter results exists
  if (!dir.exists(filter_dir)) {
    dir.create(filter_dir)
    flog.info("📁 Directory created: %s", filter_dir)
  } else {
    flog.info("📁 Directory already exists: %s", filter_dir)
  }
  
  # Capture the number of cells before filtering
  pre_filter <- seurat_obj@meta.data
  flog.info("🔍 Initial number of cells before filtering: %d", nrow(pre_filter))
  
  # Apply the filters
  flog.info("⚙️  Applying filters on nFeature_RNA, nCount_RNA, percent_mito, QC_Consensus_Filtered, and samples")
  seurat_obj <- seurat_obj %>%
    subset(nFeature_RNA >= min_nFeature_RNA & nFeature_RNA <= max_nFeature_RNA &
             nCount_RNA > min_nCount_RNA & nCount_RNA < max_nCount_RNA &
             percent_mito <= mito_high_cutoff &
             QC_Consensus_Filtered != "Doublet" &
             !seurat_obj@meta.data$sample %in% samples)
  
  # Capture the number of cells after filtering
  post_filter <- seurat_obj@meta.data
  flog.info("✅ Filtering completed. Number of cells after filtering: %d", nrow(post_filter))
  
  # Determine removed cells
  removed_cells <- pre_filter[!rownames(pre_filter) %in% rownames(post_filter),]
  flog.info("📉 Number of cells removed during filtering: %d", nrow(removed_cells))
  
  # Create data frame for removed entries with reasons
  removed_UMI <- data.frame(
    cells = rownames(removed_cells),
    filterReason = ifelse(removed_cells$nFeature_RNA < min_nFeature_RNA | removed_cells$nFeature_RNA > max_nFeature_RNA, "nFeature_RNA", 
                          ifelse(removed_cells$nCount_RNA < min_nCount_RNA | removed_cells$nCount_RNA > max_nCount_RNA, "nCount_RNA", 
                                 ifelse(removed_cells$percent_mito > mito_high_cutoff, "percent_mito", "sample_removal")))
  )
  
  # Log details of removed cells
  flog.info("💾 Writing details of removed cells to file: %s", file.path(filter_dir, "04_removed_umis.csv"))
  write.csv(removed_UMI, file.path(filter_dir, "04_removed_umis.csv"), row.names = FALSE)
  flog.info("📁 File saved: %s", file.path(filter_dir, "04_removed_umis.csv"))
  
  # Return the filtered Seurat object
  flog.info("🚀 Filter process completed. Returning filtered Seurat object.")
  return(seurat_obj)
}


process_seurat_metadata <- function(seurat_obj, 
                                    remove_names = c("pANN_0.25_0.09_1834", "RNA_snn_res.0.5"), 
                                    rename_single = c("QC_Consensus_Filtered", "QCs_Consensus"), 
                                    rename_names = c("percent_mito", "log10GenesPerUMI", "percent_top50", "percent_oxphos", "percent_apop", 
                                                     "percent_dna_repair", "percent_ieg", "S.Score", "G2M.Score", "Phase"), 
                                    rename_prefix = "QCs") {
  
  flog.info("😎 Function Name: process_seurat_metadata")
  
  # Extract metadata
  metadata <- seurat_obj@meta.data
  flog.info("🛠️ Extracted metadata with %d columns and %d rows", ncol(metadata), nrow(metadata))
  
  # 1. Remove specified columns
  if (!is.null(remove_names)) {
    columns_to_remove <- intersect(colnames(metadata), remove_names)
    if (length(columns_to_remove) > 0) {
      flog.info("❌ Removing columns: %s", paste(columns_to_remove, collapse = ", "))
      metadata <- metadata[, !colnames(metadata) %in% remove_names]
    } else {
      flog.info("⚠️ No matching columns found to remove")
    }
  }
  
  # 2. Rename a specific column
  if (!is.null(rename_single) && length(rename_single) == 2) {
    oldname <- rename_single[1]
    newname <- rename_single[2]
    if (oldname %in% colnames(metadata)) {
      flog.info("🔄 Renaming column '%s' to '%s'", oldname, newname)
      colnames(metadata)[colnames(metadata) == oldname] <- newname
    } else {
      flog.info("⚠️ Column '%s' not found for renaming", oldname)
    }
  }
  
  # 3. Add prefix to specified columns
  if (!is.null(rename_names)) {
    renamed_columns <- list()
    for (name in rename_names) {
      if (name %in% colnames(metadata)) {
        new_colname <- paste0(rename_prefix, "_", name)
        flog.info("🔄 Renaming column '%s' to '%s'", name, new_colname)
        colnames(metadata)[colnames(metadata) == name] <- new_colname
        renamed_columns <- append(renamed_columns, new_colname)
      }
    }
    if (length(renamed_columns) > 0) {
      flog.info("📝 Prefixed columns: %s", paste(renamed_columns, collapse = ", "))
    } else {
      flog.info("⚠️ No columns were renamed with prefix '%s'", rename_prefix)
    }
  }
  
  # Update the Seurat object with modified metadata
  seurat_obj@meta.data <- metadata
  flog.info("✅ Metadata processing completed. Updated metadata has %d columns", ncol(metadata))
  
  return(seurat_obj)
}




flog.info("Function Load Completed...")




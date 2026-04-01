# utils/filter.R
# Functions for filtering and splitting Seurat objects prior to integration.

library(Seurat)
library(futile.logger)

# Load a Seurat RDS, subset to cells where meta_column == target_value,
# and remove batches with fewer than minimum_cells_per_group cells.
# Returns a list: list(seurat_object = ..., summary = data.frame(...)).
FilterSeuratForIntegration <- function(
  path,
  target_value,
  meta_column,
  filter_column_name = "batch",
  minimum_cells_per_group = 20
) {
  flog.info("📦 Loading %s", path)

  seu <- readRDS(path)

  sample_name        <- basename(path)
  total_cells_main   <- ncol(seu)

  summary_df <- data.frame(
    sample               = sample_name,
    total_cells_main     = total_cells_main,
    total_cells_target   = NA,
    batches_removed      = NA,
    final_cells_retained = NA,
    status               = NA,
    stringsAsFactors     = FALSE
  )

  # Verify meta_column exists
  if (!meta_column %in% colnames(seu@meta.data)) {
    flog.warn("⚠️ Column '%s' not found in %s", meta_column, sample_name)
    summary_df$status <- "meta_column_missing"
    return(list(seurat_object = NULL, summary = summary_df))
  }

  seu[[meta_column]] <- trimws(seu[[meta_column]][, 1])

  # Verify target_value exists
  if (!target_value %in% seu[[meta_column]][, 1]) {
    flog.warn("🚫 %s → No '%s' found", sample_name, target_value)
    summary_df$status <- "target_not_found"
    return(list(seurat_object = NULL, summary = summary_df))
  }

  # Subset to target cells
  seu_subset <- tryCatch(
    subset(seu, cells = colnames(seu)[seu[[meta_column]][, 1] == target_value]),
    error = function(e) NULL
  )

  if (is.null(seu_subset) || ncol(seu_subset) == 0) {
    summary_df$status <- "subset_failed_or_empty"
    return(list(seurat_object = NULL, summary = summary_df))
  }

  total_cells_target            <- ncol(seu_subset)
  summary_df$total_cells_target <- total_cells_target

  # Remove under-represented batch groups
  batches_removed <- 0

  if (filter_column_name %in% colnames(seu_subset@meta.data)) {
    group_counts  <- table(seu_subset[[filter_column_name]][, 1])
    keep_groups   <- names(group_counts[group_counts >= minimum_cells_per_group])
    remove_groups <- names(group_counts[group_counts <  minimum_cells_per_group])

    batches_removed <- length(remove_groups)

    if (length(keep_groups) == 0) {
      summary_df$batches_removed <- batches_removed
      summary_df$status          <- "all_groups_removed"
      return(list(seurat_object = NULL, summary = summary_df))
    }

    seu_subset <- subset(
      seu_subset,
      cells = colnames(seu_subset)[seu_subset[[filter_column_name]][, 1] %in% keep_groups]
    )
  }

  summary_df$batches_removed      <- batches_removed
  summary_df$final_cells_retained <- ncol(seu_subset)
  summary_df$status               <- "success"

  flog.info("🧬 %s → Final cells retained: %d", sample_name, ncol(seu_subset))

  return(list(seurat_object = seu_subset, summary = summary_df))
}


# Subset a Seurat object to cells belonging to subset_group_name,
# then split the result into a list of Seurat objects by batch_column_name.
splitSeuratObjectToList <- function(seuratObject, subset_group_name, batch_column_name = "batch") {
  flog.info("Starting the Seurat object processing function.")

  tryCatch({
    if (!"seurat_clusters_labeled" %in% colnames(seuratObject@meta.data)) {
      stop("'seurat_clusters_labeled' column is not present in the metadata.")
    }

    flog.info("Subsetting the Seurat object with subset_group_name: %s", subset_group_name)
    subsetSeuratObject <- subset(x = seuratObject,
                                  subset = seurat_clusters_labeled == subset_group_name,
                                  invert = FALSE)

    flog.info("Extracting the count matrix from the subset Seurat object.")
    subset.countMatrix <- LayerData(subsetSeuratObject, assay = "RNA", layer = "counts")

    flog.info("Extracting metadata from the subset Seurat object.")
    subset.meta.data <- subsetSeuratObject@meta.data

    flog.info("Creating a new Seurat object with the subset count matrix and metadata.")
    subsetSeuratObject <- CreateSeuratObject(
      counts   = subset.countMatrix,
      project  = subset_group_name,
      meta.data = subset.meta.data
    )

    if (!batch_column_name %in% colnames(subsetSeuratObject@meta.data)) {
      stop(sprintf("The batch column '%s' is not present in the metadata.", batch_column_name))
    }

    flog.info("Splitting the Seurat object by batch column: %s", batch_column_name)
    splitSeuratList <- split(subsetSeuratObject, f = subsetSeuratObject$batch)

    flog.info("Seurat object processing completed.")
    return(splitSeuratList)

  }, error = function(e) {
    flog.error("Error occurred: %s", e$message)
    stop("The splitSeuratObjectToList function failed. Please check the input parameters and try again.")
  })
}

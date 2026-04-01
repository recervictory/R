# utils/consensus.R
# Functions for deriving consensus cell-type and meta-program labels across clusters.

library(Seurat)
library(dplyr)
library(futile.logger)

# Smooth CT_consensus labels by assigning the cluster-level mode to every cell.
# Returns a data frame with CT_consensus and CT_consensus_cluster columns.
create_consensus_CT_clusters <- function(seuratObject, cutoff = 0, resolution = 20) {
  seuratObject <- FindClusters(seuratObject, resolution = resolution)

  meta_data <- data.frame(
    CT_consensus    = seuratObject$CT_consensus,
    seurat_clusters = seuratObject$seurat_clusters
  )

  # Calculate the mode of a vector
  get_mode <- function(v) {
    uniq_v <- unique(v)
    uniq_v[which.max(tabulate(match(v, uniq_v)))]
  }

  new_df <- meta_data %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::mutate(CT_consensus_cluster = get_mode(CT_consensus)) %>%
    dplyr::ungroup()

  return(new_df)
}


# Derive consensus meta-program (MP) labels for each cell by finding the top-3
# scoring MPs per cell, then smoothing by cluster-level mode.
# Returns a data frame with MPs_consensus, MPs_Second_consensus, MPs_Third_consensus columns.
create_consensus_Mps <- function(seuratObject, cutoff = 0.2, resolution = 20, cutoff_list = NULL) {
  seuratObject$MPs_consensus_Score        <- NULL
  seuratObject$MPs_Second_consensus_Score <- NULL
  seuratObject$MPs_Third_consensus_Score  <- NULL

  seuratObject <- FindClusters(seuratObject, resolution = resolution)

  prefixes  <- c("MPs")
  meta_data <- seuratObject@meta.data

  # Identify MP score columns
  mp.col.names <- colnames(meta_data)[sapply(colnames(meta_data), function(col) any(startsWith(col, prefixes)))]

  df <- meta_data[, mp.col.names, drop = FALSE]
  df <- df[sapply(df, is.numeric)]

  # Apply per-column cutoffs (set values below threshold to NA)
  if (!is.null(cutoff_list)) {
    for (col in names(cutoff_list)) {
      if (col %in% colnames(df)) {
        df[[col]][df[[col]] < cutoff_list[[col]]] <- NA
      }
    }
  }

  # Pre-allocate result vectors
  n_cells              <- nrow(df)
  max_scores           <- numeric(n_cells)
  max_cols             <- character(n_cells)
  second_max_scores    <- numeric(n_cells)
  second_max_cols      <- character(n_cells)
  third_max_scores     <- numeric(n_cells)
  third_max_cols       <- character(n_cells)
  confidence           <- character(n_cells)

  for (i in seq_len(n_cells)) {
    row            <- as.numeric(df[i, ])
    sorted_indices <- order(row, decreasing = TRUE, na.last = NA)
    sorted_values  <- row[sorted_indices]

    max_scores[i]        <- sorted_values[1]
    max_cols[i]          <- names(df)[sorted_indices[1]]
    second_max_scores[i] <- sorted_values[2]
    second_max_cols[i]   <- names(df)[sorted_indices[2]]
    third_max_scores[i]  <- sorted_values[3]
    third_max_cols[i]    <- names(df)[sorted_indices[3]]
    confidence[i]        <- ifelse(sorted_values[1] > cutoff, "Confident", "Not Confident")
  }

  # Calculate the mode of a vector
  get_mode <- function(v) {
    uniq_v <- unique(v)
    uniq_v[which.max(tabulate(match(v, uniq_v)))]
  }

  new_df <- data.frame(
    MPs_consensus_Score       = max_scores,
    MPs_consensus             = sub("^MPs_", "", max_cols),
    MPs_Second_consensus_Score = second_max_scores,
    MPs_Second_consensus      = sub("^MPs_", "", second_max_cols),
    MPs_Third_consensus_Score = third_max_scores,
    MPs_Third_consensus       = sub("^MPs_", "", third_max_cols),
    MPs_consensus_Confidence  = confidence,
    MPs_consensus_cluster     = seuratObject$seurat_clusters
  )

  new_df <- new_df %>%
    dplyr::group_by(MPs_consensus_cluster) %>%
    dplyr::mutate(
      MPs_consensus_cluster_program = get_mode(MPs_consensus),
      MPs_Second_cluster_program    = get_mode(MPs_Second_consensus),
      MPs_Third_cluster_program     = get_mode(MPs_Third_consensus)
    ) %>%
    dplyr::ungroup()

  return(new_df)
}

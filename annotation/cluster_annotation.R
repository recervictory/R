# annotation/cluster_annotation.R
# Functions for assigning human-readable labels to Seurat clusters.

library(Seurat)
library(futile.logger)

# Assign user-defined names to Seurat clusters and store them in a new
# metadata column. Clusters not listed in cluster_list receive other_clusters_name.
name_seurat_cluster <- function(seurat_object, cluster_list,
                                other_clusters_name = "Unknown",
                                meta.data_colname = "seurat_clusters_labeled") {
  cluster_assignments <- rep(other_clusters_name, ncol(seurat_object))
  unique_clusters <- unique(seurat_object$seurat_clusters)

  for (cluster_name in names(cluster_list)) {
    clusters <- cluster_list[[cluster_name]]
    cluster_assignments[seurat_object$seurat_clusters %in% clusters] <- cluster_name
    flog.info(paste("Assigning cluster:", cluster_name, "with cluster numbers:", paste(clusters, collapse = ", ")))
  }

  not_found_clusters <- setdiff(unique_clusters, unlist(cluster_list))
  if (length(not_found_clusters) > 0) {
    flog.info(paste("Assigning cluster: Unknown with cluster numbers:", paste(not_found_clusters, collapse = ", ")))
  } else {
    flog.info("All clusters present in cluster_list.")
  }

  seurat_object[[meta.data_colname]] <- cluster_assignments
  return(seurat_object)
}


# Assign user-defined names to Seurat clusters using annotate_seurat_cluster.
# Cells not matched by cluster_list receive other_clusters_name (default "Malignant").
annotate_seurat_cluster <- function(seurat_object, cluster_list,
                                    other_clusters_name = "Malignant",
                                    meta.data_colname = "seurat_clusters_labeled") {
  cluster_assignments <- rep(other_clusters_name, ncol(seurat_object))
  unique_clusters <- unique(seurat_object$seurat_clusters)

  for (cluster_name in names(cluster_list)) {
    clusters <- cluster_list[[cluster_name]]
    cluster_assignments[seurat_object$seurat_clusters %in% clusters] <- cluster_name
    flog.info(paste("Assigning cluster:", cluster_name, "with cluster numbers:", paste(clusters, collapse = ", ")))
  }

  not_found_clusters <- setdiff(unique_clusters, unlist(cluster_list))
  if (length(not_found_clusters) > 0) {
    flog.info(paste("Assigning cluster: Unknown with cluster numbers:", paste(not_found_clusters, collapse = ", ")))
  } else {
    flog.info("All clusters present in cluster_list.")
  }

  seurat_object[[meta.data_colname]] <- cluster_assignments
  return(seurat_object)
}


# Build an annotation data frame by taking the majority label in each cluster.
# col_names should include the grouping column (by) and at least one annotation column.
create_annotation <- function(seuratObject, col_names, by = "seurat_clusters") {
  meta.data <- seuratObject@meta.data
  annotation <- meta.data[, col_names]

  annotation <- annotation %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(col_names))) %>%
    dplyr::summarize(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  annotation <- as.data.frame(annotation)
  annotation[[by]] <- as.character(annotation[[by]])
  row.names(annotation) <- annotation[[by]]
  annotation[[by]] <- NULL
  annotation$count <- NULL

  return(annotation)
}

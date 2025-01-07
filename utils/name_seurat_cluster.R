name_seurat_cluster <- function(seurat_object, cluster_list, other_clusters_name = "Unknown", meta.data_colname = 'seurat_clusters_labeled') {
  
  # Initialize a vector for storing cluster assignments
  cluster_assignments <- rep(other_clusters_name, ncol(seurat_object))  # Default to "Unknown"
  
  # Get the unique cluster numbers from the Seurat object
  unique_clusters <- unique(seurat_object$seurat_clusters)
  
  # Assign clusters based on the cluster list
  for (cluster_name in names(cluster_list)) {
    # Get the cluster numbers for the current named cluster
    clusters <- cluster_list[[cluster_name]]
    
    # Assign the cluster name to the corresponding cells in the cluster_assignments vector
    cluster_assignments[seurat_object$seurat_clusters %in% clusters] <- cluster_name
    
    # Log information about the clusters being assigned
    flog.info(paste("Assigning cluster:", cluster_name, "with cluster numbers:", paste(clusters, collapse = ", ")))
  }
  
  # Log any cluster numbers in the Seurat object not present in the cluster_list
  not_found_clusters <- setdiff(unique_clusters, unlist(cluster_list))
  
  if (length(not_found_clusters) > 0) {
    flog.info(paste("Assigning cluster: Unkown with cluster numbers:", paste(not_found_clusters, collapse = ", ")))
  } else {
    flog.info("All clusters present in cluster_list.")
  }
  
  # Add the new cluster assignments as a metadata column in the Seurat object
  seurat_object[[meta.data_colname]] <- cluster_assignments
  
  # Return the modified Seurat object
  return(seurat_object)
}

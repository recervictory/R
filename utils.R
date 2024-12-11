library(Seurat)
library(futile.logger)


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

# -------------------------------

splitSeuratObjectToList <- function(seuratObject, subset_group_name, batch_column_name = "batch") {
  # Start logging
  flog.info("Starting the Seurat object processing function.")
  
  tryCatch({
    # Check if 'seurat_clusters_labeled' exists in metadata
    if (!"seurat_clusters_labeled" %in% colnames(seuratObject@meta.data)) {
      stop("'seurat_clusters_labeled' column is not present in the metadata.")
    }
    
    
    
    # Step 1: Subset Seurat object
    flog.info("Subsetting the Seurat object with subset_group_name: %s", subset_group_name)
    subsetSeuratObject <- subset(x = seuratObject, subset = seurat_clusters_labeled == subset_group_name, invert = FALSE)
    
    # Step 2: Extract count matrix
    flog.info("Extracting the count matrix from the subset Seurat object.")
    subset.countMatrix <- LayerData(subsetSeuratObject, assay = "RNA", layer = "counts")
    
    # Step 3: Extract metadata
    flog.info("Extracting metadata from the subset Seurat object.")
    subset.meta.data <- subsetSeuratObject@meta.data
    
    # Step 4: Create a new Seurat object
    flog.info("Creating a new Seurat object with the subset count matrix and metadata.")
    subsetSeuratObject <- CreateSeuratObject(
      counts = subset.countMatrix,
      project = subset_group_name,
      meta.data = subset.meta.data
    )
    
    # Step 5: Check if the batch column exists
    if (!batch_column_name %in% colnames(subsetSeuratObject@meta.data)) {
      stop(sprintf("The batch column '%s' is not present in the metadata.", batch_column_name))
    }
    
    # Step 6: Split the Seurat object by batch
    flog.info("Splitting the Seurat object by batch column: %s", batch_column_name)
    splitSeuratList <- split(subsetSeuratObject, f = subsetSeuratObject$batch)
    
    # End logging
    flog.info("Seurat object processing completed.")
    
    # Return the split list of Seurat objects
    return(splitSeuratList)
    
  }, error = function(e) {
    flog.error("Error occurred: %s", e$message)
    stop("The splitSeuratObjectToList function failed. Please check the input parameters and try again.")
  })
}


processSeuratForIntegrated <- function(
  seuratObject, 
  dims = 1:30, 
  resolution = 1, 
  cluster.name = "unintegrated_clusters", 
  reduction.name = "umap_unintegrated"
) {
  # Start logging
  flog.info("Starting the Seurat object processing and clustering function.")
  
  tryCatch({
    # Step 1: Normalize Data
    flog.info("Normalizing the data.")
    seuratObject <- NormalizeData(seuratObject)
    
    # Step 2: Find Variable Features
    flog.info("Finding variable features.")
    seuratObject <- FindVariableFeatures(seuratObject)
    
    # Step 3: Scale Data
    flog.info("Scaling the data.")
    seuratObject <- ScaleData(seuratObject)
    
    # Step 4: Run PCA
    flog.info("Running PCA with dimensions: %s", paste(dims, collapse = ", "))
    seuratObject <- RunPCA(seuratObject)
    
    # Step 5: Find Neighbors
    flog.info("Finding neighbors using PCA reduction and dimensions: %s", paste(dims, collapse = ", "))
    seuratObject <- FindNeighbors(seuratObject, dims = dims, reduction = "pca")
    
    # Step 6: Find Clusters
    flog.info("Finding clusters with resolution: %s and cluster name: %s", resolution, cluster.name)
    seuratObject <- FindClusters(seuratObject, resolution = resolution)
    seuratObject@meta.data[[cluster.name]] <- seuratObject$seurat_clusters
    
    # Step 7: Run UMAP
    flog.info("Running UMAP with dimensions: %s and reduction name: %s", paste(dims, collapse = ", "), reduction.name)
    seuratObject <- RunUMAP(seuratObject, dims = dims, reduction = "pca", reduction.name = reduction.name)
    
    # Final Logging
    flog.info("Seurat object processing and clustering completed.")
    
    # Return the processed Seurat object
    return(seuratObject)
    
  }, error = function(e) {
    flog.error("Error occurred: %s", e$message)
    stop("The processAndClusterSeurat function failed. Please check the input parameters and try again.")
  })
}



integrateAndClusterSeurat <- function(
  seuratObject, 
  integrationMethod, 
  orig.reduction = "pca", 
  new.reduction = "integrated", 
  dims = 1:30, 
  resolution = 0.5
) {
  # Start logging
  flog.info("Starting the Seurat object integration and clustering function.")
  
  tryCatch({
    # Step 1: Integrate Layers
    flog.info("Integrating layers using method: %s, original reduction: %s, new reduction: %s", 
              deparse(substitute(integrationMethod)), orig.reduction, new.reduction)
    seuratObject <- IntegrateLayers(
      object = seuratObject, 
      method = integrationMethod, 
      orig.reduction = orig.reduction, 
      new.reduction = new.reduction,
      verbose = TRUE
    )
    
    # Step 2: Re-join layers after integration
    flog.info("Re-joining layers after integration.")
    seuratObject[["RNA"]] <- JoinLayers(seuratObject[["RNA"]])
    
    # Step 3: Find Neighbors
    flog.info("Finding neighbors using integrated reduction and dimensions: %s", paste(dims, collapse = ", "))
    seuratObject <- FindNeighbors(seuratObject, reduction = new.reduction, dims = dims)
    
    # Step 4: Find Clusters
    flog.info("Finding clusters with resolution: %s", resolution)
    seuratObject <- FindClusters(seuratObject, resolution = resolution)
    
    # Step 5: Run UMAP
    flog.info("Running UMAP with dimensions: %s and reduction: %s", paste(dims, collapse = ", "), new.reduction)
    seuratObject <- RunUMAP(seuratObject, dims = dims, reduction = new.reduction)
    
    # Final Logging
    flog.info("Seurat object integration and clustering completed.")
    
    # Return the processed Seurat object
    return(seuratObject)
    
  }, error = function(e) {
    flog.error("Error occurred: %s", e$message)
    stop("The integrateAndClusterSeurat function failed. Please check the input parameters and try again.")
  })
}

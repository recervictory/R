# core/seurat_workflow.R
# Core Seurat analysis workflow functions: normalization, PCA, UMAP, clustering.

library(Seurat)
library(futile.logger)

# Run the full Seurat preprocessing pipeline on a Seurat object.
# Skips the workflow if UMAP already exists.
runSeuratWorkflow <- function(object) {
  flog.info("😎 Function Name: runSeuratWorkflow")

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


# Normalize, find variable features, scale, run PCA, and find neighbors.
# Returns the Seurat object ready for clustering or integration.
process_seurat_object <- function(seurat_object, dims = 1:30, nfeatures = 3000) {
  flog.info("NormalizeData Started...")
  seurat_object <- NormalizeData(seurat_object)
  flog.info("NormalizeData Completed...")

  flog.info("FindVariableFeatures Started...")
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  flog.info("FindVariableFeatures Completed...")

  flog.info("ScaleData Started...")
  seurat_object <- ScaleData(seurat_object)
  flog.info("ScaleData Completed...")

  flog.info("RunPCA Started...")
  seurat_object <- RunPCA(seurat_object)
  flog.info("RunPCA Completed...")

  flog.info("FindNeighbors Started...")
  seurat_object <- FindNeighbors(seurat_object, dims = dims, reduction = "pca")
  flog.info("FindNeighbors Completed...")

  return(seurat_object)
}


# Run the unintegrated Seurat pipeline: normalize, scale, PCA, neighbors, clusters, UMAP.
# Stores clusters under cluster.name and UMAP under reduction.name.
processSeuratForIntegrated <- function(
  seuratObject,
  dims = 1:30,
  resolution = 1,
  cluster.name = "unintegrated_clusters",
  reduction.name = "umap_unintegrated"
) {
  flog.info("Starting the Seurat object processing and clustering function.")

  tryCatch({
    flog.info("Normalizing the data.")
    seuratObject <- NormalizeData(seuratObject)

    flog.info("Finding variable features.")
    seuratObject <- FindVariableFeatures(seuratObject)

    flog.info("Scaling the data.")
    seuratObject <- ScaleData(seuratObject)

    flog.info("Running PCA with dimensions: %s", paste(dims, collapse = ", "))
    seuratObject <- RunPCA(seuratObject)

    flog.info("Finding neighbors using PCA reduction and dimensions: %s", paste(dims, collapse = ", "))
    seuratObject <- FindNeighbors(seuratObject, dims = dims, reduction = "pca")

    flog.info("Finding clusters with resolution: %s and cluster name: %s", resolution, cluster.name)
    seuratObject <- FindClusters(seuratObject, resolution = resolution)
    seuratObject@meta.data[[cluster.name]] <- seuratObject$seurat_clusters

    flog.info("Running UMAP with dimensions: %s and reduction name: %s", paste(dims, collapse = ", "), reduction.name)
    seuratObject <- RunUMAP(seuratObject, dims = dims, reduction = "pca", reduction.name = reduction.name)

    flog.info("Seurat object processing and clustering completed.")
    return(seuratObject)

  }, error = function(e) {
    flog.error("Error occurred: %s", e$message)
    stop("The processSeuratForIntegrated function failed. Please check the input parameters and try again.")
  })
}

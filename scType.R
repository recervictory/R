

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R") 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")



# Function definition
scType_cell_prediction <- function(seuratObject, 
                                   databasepat, 
                                   tissue, 
                                   metedatadata_col_name, 
                                   overwrite = TRUE,
                                   assay = "RNA",
                                   verbose = TRUE) {
  # Setup logging
  library(logging)
  
  if (verbose) {
    basicConfig(level = 'INFO')
  } else {
    basicConfig(level = 'ERROR')
  }
  
  flog.info("Starting the processing of Seurat object.")
  
  # Prepare gene sets
  flog.info("Preparing gene sets from database.")
  gs_list <- gene_sets_prepare(databasepat, tissue)
  flog.info("Gene sets prepared.")
  
  # Check for Seurat package version
  seurat_package_v5 <- isFALSE('scale.data' %in% names(attributes(seuratObject[["RNA"]])))
  flog.info(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))
  
  # Extract scaled data
  flog.info("Extracting scaled data.")
  scRNAseqData_scaled <- if (seurat_package_v5) {
    as.matrix(GetAssayData(seuratObject, layer = "scale.data")) #[["RNA"]]$scale.data
  } else {
    as.matrix(seuratObject[[assay]]@scale.data)
  }
  
  
  # Calculate scores
  flog.info("Calculating SCTYPE scores.")
  es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  flog.info("SCTYPE scores calculated.")
  
  # Process results
  flog.info("Processing SCTYPE results.")
  cL_resutls <- do.call("rbind", lapply(unique(seuratObject@meta.data$seurat_clusters), function(cl) {
    es.max.cl <- sort(rowSums(es.max[ ,rownames(seuratObject@meta.data[seuratObject@meta.data$seurat_clusters == cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seuratObject@meta.data$seurat_clusters == cl)), 10)
  }))
  sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  flog.info("SCTYPE results processed.")
  
  # Set low-confident clusters to "Unknown"
  flog.info("Updating low-confidence clusters.")
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- "Unknown"
  
  # Check if metadadata_col_name exists
  if (metedatadata_col_name %in% colnames(seuratObject@meta.data)) {
    if (overwrite) {
      flog.warn(sprintf("Column '%s' already exists. Overwriting it.", metedatadata_col_name))
      seuratObject@meta.data[[metedatadata_col_name]] <- NULL  # Set existing column to NULL
    } else {
      flog.error(sprintf("Column '%s' already exists. Set 'overwrite' to TRUE to overwrite it.", metedatadata_col_name))
      stop(sprintf("Column '%s' already exists. Not overwriting.", metedatadata_col_name))
    }
  }
  
  # Assign new cell type information to the metadata
  flog.info("Assigning new cell type information to Seurat object.")
  seuratObject@meta.data[[metedatadata_col_name]] <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ] 
    seuratObject@meta.data[[metedatadata_col_name]][seuratObject@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
  }
  
  flog.info("Seurat object processing completed.")
  
  return(seuratObject)
}

### ClusterTree Object plots

run_clustree <- function(seuratObject, resolutions = c(0.1, 0.25, 0.5, 1.0, 1.5, 2), log_level = "INFO") {
  library(futile.logger)
  library(clustree)
  
  
  # Set logging level
  flog.threshold(log_level)
  
  clustree_seurat_object <- seuratObject
  
  for (res in resolutions) {
    flog.info(paste("FindClusters Started with resolution", res, "..."))
    
    # Perform clustering
    clustree_seurat_object <- FindClusters(clustree_seurat_object, graph.name = "RNA_snn", resolution = res, algorithm = 1, verbose = FALSE)
    
    flog.info(paste("FindClusters Completed with resolution", res, "..."))
  }
  
  # Create the clustree plot
  plot <- clustree(clustree_seurat_object@meta.data, prefix = "RNA_snn_res.")
  
  return(plot)
}


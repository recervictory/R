# cnv/infercnv.R
# Helper function for setting up inferCNV reference cell annotation.

library(Seurat)

# Create an infer_cnv_reference metadata column by using primaryColumn values
# where they match conditionValues, and secondaryColumn values otherwise.
updateSeuratMetadata <- function(seuratObject, primaryColumn, secondaryColumn, conditionValues) {
  metadata <- seuratObject@meta.data

  if (!primaryColumn %in% colnames(metadata) || !secondaryColumn %in% colnames(metadata)) {
    stop("One or both specified columns do not exist in the Seurat object metadata.")
  }

  seuratObject$infer_cnv_reference <- ifelse(
    metadata[[primaryColumn]] %in% conditionValues,
    metadata[[primaryColumn]],
    metadata[[secondaryColumn]]
  )

  print(table(seuratObject$infer_cnv_reference))
  return(seuratObject)
}

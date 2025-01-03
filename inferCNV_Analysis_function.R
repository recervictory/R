source("https://raw.githubusercontent.com/recervictory/R/refs/heads/master/utils.R")

# Update Seurat Metadata: Create a new metadata column based on conditional logic
updateSeuratMetadata <- function(seuratObject, primaryColumn, secondaryColumn, conditionValues) {
  
  # Retrieve metadata from the Seurat object
  metadata <- seuratObject@meta.data
  
  # Check if the specified columns exist in the metadata
  if (!primaryColumn %in% colnames(metadata) || !secondaryColumn %in% colnames(metadata)) {
    stop("One or both specified columns do not exist in the Seurat object metadata.")
  }
  
  # Create a new column 'infer_cnv_reference' based on the specified condition
  # If the value in the primary column exists in 'conditionValues', use it;
  # otherwise, use the corresponding value from the secondary column
  seuratObject$infer_cnv_reference <- ifelse(metadata[[primaryColumn]] %in% conditionValues, 
                                             metadata[[primaryColumn]], 
                                             metadata[[secondaryColumn]])
  
  # Display a summary table of the newly created column
  print(table(seuratObject$infer_cnv_reference))
  
  # Return the updated Seurat object
  return(seuratObject)
}


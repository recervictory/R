# utils/metadata.R
# Utility functions for manipulating Seurat object metadata columns.

library(Seurat)
library(futile.logger)

# Rename and prefix metadata columns, and remove unwanted columns in one step.
# remove_names:   columns to drop.
# rename_single:  c(old_name, new_name) pair.
# rename_names:   columns that receive rename_prefix as a prefix.
process_seurat_metadata <- function(seurat_obj,
                                    remove_names = c("pANN_0.25_0.09_1834", "RNA_snn_res.0.5"),
                                    rename_single = c("QC_Consensus_Filtered", "QCs_Consensus"),
                                    rename_names  = c("percent_mito", "log10GenesPerUMI", "percent_top50",
                                                      "percent_oxphos", "percent_apop", "percent_dna_repair",
                                                      "percent_ieg", "S.Score", "G2M.Score", "Phase"),
                                    rename_prefix = "QCs") {
  flog.info("😎 Function Name: process_seurat_metadata")

  metadata <- seurat_obj@meta.data
  flog.info("🛠️ Extracted metadata with %d columns and %d rows", ncol(metadata), nrow(metadata))

  # Remove specified columns
  if (!is.null(remove_names)) {
    columns_to_remove <- intersect(colnames(metadata), remove_names)
    if (length(columns_to_remove) > 0) {
      flog.info("❌ Removing columns: %s", paste(columns_to_remove, collapse = ", "))
      metadata <- metadata[, !colnames(metadata) %in% remove_names]
    } else {
      flog.info("⚠️ No matching columns found to remove")
    }
  }

  # Rename a single column
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

  # Add prefix to specified columns
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

  seurat_obj@meta.data <- metadata
  flog.info("✅ Metadata processing completed. Updated metadata has %d columns", ncol(metadata))

  return(seurat_obj)
}


# Create a new metadata column by conditionally copying values from
# source_metadata_colname where they match match_values, otherwise using
# values from destination_metadata_colname.
transferMetadata <- function(seurat_object,
                              source_metadata_colname,
                              destination_metadata_colname,
                              new_metadata_colname,
                              match_values,
                              verbose = TRUE) {
  flog.info("😎 Function Name: transferMetadata")
  flog.info("🔍 Starting metadata transfer from '%s' to '%s'.", source_metadata_colname, new_metadata_colname)

  metadata <- seurat_object@meta.data

  if (!source_metadata_colname %in% colnames(metadata)) {
    flog.error("❌ Source metadata column '%s' does not exist in the Seurat object.", source_metadata_colname)
    stop(sprintf("Source metadata column '%s' not found.", source_metadata_colname))
  }

  if (!destination_metadata_colname %in% colnames(metadata)) {
    flog.error("❌ Destination metadata column '%s' does not exist in the Seurat object.", destination_metadata_colname)
    stop(sprintf("Destination metadata column '%s' not found.", destination_metadata_colname))
  }

  flog.info("✅ Source and destination columns exist.")

  flog.info("📝 Creating new metadata column '%s' based on match values.", new_metadata_colname)
  metadata[[new_metadata_colname]] <- ifelse(
    metadata[[source_metadata_colname]] %in% match_values,
    metadata[[source_metadata_colname]],
    metadata[[destination_metadata_colname]]
  )

  flog.info("🔄 Updating Seurat object with new metadata column: '%s'.", new_metadata_colname)
  seurat_object@meta.data <- metadata

  flog.info("🎉 Metadata transfer completed successfully.")
  return(seurat_object)
}


# Replace NA values in a metadata column with a specified replacement value.
replace_seurat_NA <- function(seurat_obj, col_name, replacement_value = "Unknown") {
  if (!(col_name %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", col_name, "does not exist in the Seurat object metadata."))
  }

  na_count <- sum(is.na(seurat_obj@meta.data[[col_name]]))
  flog.info("NA Found : %s", na_count)

  seurat_obj@meta.data[[col_name]][is.na(seurat_obj@meta.data[[col_name]])] <- replacement_value
  return(seurat_obj)
}


# Map values in old_column to a new_column using a named mapping vector.
# Values not found in mapping are kept as-is.
update_seurat_metadata <- function(seurat_obj, old_column,
                                   new_column = "PBTA_Annotation",
                                   mapping) {
  if (!old_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", old_column, "not found in metadata"))
  }

  seurat_obj[[new_column]] <- seurat_obj@meta.data[[old_column]]

  seurat_obj@meta.data[[new_column]] <- sapply(seurat_obj@meta.data[[new_column]], function(x) {
    if (x %in% names(mapping)) mapping[x] else x
  })

  return(seurat_obj)
}

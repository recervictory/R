# utils/io.R
# I/O helpers for loading Seurat objects and reading CSV configuration files.

library(Seurat)
library(data.table)
library(futile.logger)

# Load a Seurat RDS file from data_dir.
# The file name is constructed as "<project_name>.<batch_name>.<file_name>" or
# "<project_name>.<file_name>" when batch_name is NA.
load_seurat_object <- function(file_name, project_name = NA, batch_name = NA, data_dir = "data") {
  flog.info("Project Name: %s", project_name)
  flog.info("Batch Name: %s", batch_name)

  project_name_cleaned <- gsub(" ", "_", project_name)

  if (is.na(batch_name)) {
    file_name <- paste(project_name_cleaned, file_name, sep = ".")
  } else {
    batch_name_cleaned <- gsub(" ", "_", batch_name)
    file_name <- paste(project_name_cleaned, batch_name_cleaned, file_name, sep = ".")
  }

  flog.info("Seurat Object File Name: %s", file_name)

  file_path <- file.path(data_dir, file_name)
  flog.info("Loading Seurat Object from: %s", file_path)

  seurat_object <- readRDS(file_path)
  flog.info("Loading Seurat Object from *** %s *** Completed.", file_path)

  return(seurat_object)
}


# Read a two-column CSV (columns: filter, value) and assign each row as a
# global variable named after the filter column value.
read_csv_and_create_variables <- function(file_name) {
  flog.info("😎 Function Name: read_csv_and_create_variables")
  data <- fread(file_name)

  for (i in seq_len(nrow(data))) {
    assign(data$filter[i], data$value[i], envir = .GlobalEnv)
  }

  flog.info("Variables have been recreated from the CSV file '%s'.", file_name)
}

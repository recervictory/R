# core/qc.R
# Quality control filtering functions for Seurat objects.

library(Seurat)
library(futile.logger)

# Annotate cells in a Seurat object with QC consensus filter labels.
# Adds a QC_Consensus_Filtered column indicating the reason a cell failed (or "Passed").
filter_qc <- function(seurat_obj,
                      min_nFeature_RNA = 200, max_nFeature_RNA = 9000,
                      min_nCount_RNA = 1000, max_nCount_RNA = 45000,
                      mito_high_cutoff = 25, ...) {
  seurat_obj@meta.data$QC_Consensus_Filtered <- "Passed"

  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$nCount_RNA <= min_nCount_RNA] <- "Low nCOUNT"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$nFeature_RNA <= min_nFeature_RNA] <- "Low Genes per Cell"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$QCs_miQC_keep == FALSE] <- "MiQC+"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$nFeature_RNA >= max_nFeature_RNA] <- "High Genes per Cell"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$nCount_RNA >= max_nCount_RNA] <- "High nCOUNT"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$QCs_doubletFinder == "doublet"] <- "DoubletFinder+"
  seurat_obj@meta.data$QC_Consensus_Filtered[seurat_obj@meta.data$QCs_scDblFinder_class == "doublet"] <- "DblFinder+"
  seurat_obj@meta.data$QC_Consensus_Filtered[
    seurat_obj@meta.data$QCs_doubletFinder == "doublet" &
    seurat_obj@meta.data$QCs_scDblFinder_class == "doublet"
  ] <- "Doublet"

  return(seurat_obj)
}


# Subset a Seurat object based on QC thresholds, log removed cells,
# and write a CSV of removed UMIs with their filter reasons.
filterSeurat <- function(seurat_obj,
                         min_nFeature_RNA = 200,
                         max_nFeature_RNA = 9000,
                         min_nCount_RNA = 1000,
                         max_nCount_RNA = 45000,
                         mito_high_cutoff = 25,
                         samples = c(),
                         filter_dir = "../01_metadata/") {
  flog.info("😎 Function Name: filterSeurat")
  flog.info("🛠️  Parameters - min_nFeature_RNA: %d, max_nFeature_RNA: %d, min_nCount_RNA: %d, max_nCount_RNA: %d, mito_high_cutoff: %d",
            min_nFeature_RNA, max_nFeature_RNA, min_nCount_RNA, max_nCount_RNA, mito_high_cutoff)

  if (!dir.exists(filter_dir)) {
    dir.create(filter_dir)
    flog.info("📁 Directory created: %s", filter_dir)
  } else {
    flog.info("📁 Directory already exists: %s", filter_dir)
  }

  pre_filter <- seurat_obj@meta.data
  flog.info("🔍 Initial number of cells before filtering: %d", nrow(pre_filter))

  flog.info("⚙️  Applying filters on nFeature_RNA, nCount_RNA, percent_mito, QC_Consensus_Filtered, and samples")
  seurat_obj <- seurat_obj %>%
    subset(nFeature_RNA >= min_nFeature_RNA & nFeature_RNA <= max_nFeature_RNA &
             nCount_RNA > min_nCount_RNA & nCount_RNA < max_nCount_RNA &
             percent_mito <= mito_high_cutoff &
             QC_Consensus_Filtered != "Doublet" &
             !seurat_obj@meta.data$sample %in% samples)

  post_filter <- seurat_obj@meta.data
  flog.info("✅ Filtering completed. Number of cells after filtering: %d", nrow(post_filter))

  removed_cells <- pre_filter[!rownames(pre_filter) %in% rownames(post_filter), ]
  flog.info("📉 Number of cells removed during filtering: %d", nrow(removed_cells))

  removed_UMI <- data.frame(
    cells = rownames(removed_cells),
    filterReason = ifelse(
      removed_cells$nFeature_RNA < min_nFeature_RNA | removed_cells$nFeature_RNA > max_nFeature_RNA, "nFeature_RNA",
      ifelse(
        removed_cells$nCount_RNA < min_nCount_RNA | removed_cells$nCount_RNA > max_nCount_RNA, "nCount_RNA",
        ifelse(removed_cells$percent_mito > mito_high_cutoff, "percent_mito", "sample_removal")
      )
    )
  )

  out_file <- file.path(filter_dir, "04_removed_umis.csv")
  flog.info("💾 Writing details of removed cells to file: %s", out_file)
  write.csv(removed_UMI, out_file, row.names = FALSE)
  flog.info("📁 File saved: %s", out_file)

  flog.info("🚀 Filter process completed. Returning filtered Seurat object.")
  return(seurat_obj)
}

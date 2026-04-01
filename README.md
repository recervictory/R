# scRNA Analysis Toolkit

A collection of R functions for single-cell RNA-seq analysis built around
[Seurat](https://seuratproject.org/).

---

## Project Structure

```
R/
├── index.R                  # Entry point — sources all modules
│
├── core/                    # Core Seurat analysis workflows
│   ├── seurat_workflow.R    # runSeuratWorkflow, process_seurat_object, processSeuratForIntegrated
│   ├── integration.R        # integrateAndClusterSeurat (with retry logic)
│   └── qc.R                 # filter_qc, filterSeurat
│
├── annotation/              # Cell-type annotation
│   ├── cluster_annotation.R # name_seurat_cluster, annotate_seurat_cluster, create_annotation
│   ├── sctype.R             # scType_cell_prediction, run_clustree
│   └── mappings.R           # mapping_singleR_cluster, mapping_scTypeBrain,
│                            #   mapping_panglaoBrain, mapping_Tirosh_MP
│
├── plotting/                # Visualization helpers
│   ├── save_plot.R          # save_plot (PNG/PDF/JPEG/TIFF, Cairo-aware)
│   ├── stacked_bar.R        # plot_metadata_stacked_bar
│   ├── cluster_plots.R      # purity_df, normalized_table,
│   │                        #   cluster_purity_split_plots, clusters_meta_visualization
│   └── gene_trends.R        # plot_gene_trends
│
├── data/                    # Static reference data
│   ├── canonical_markers.R  # pbta_marker_genes — curated marker gene lists
│   └── color_codes.R        # Named color vectors + create_color_palette_for_seurat,
│                            #   set_annotation_colors
│
├── utils/                   # General-purpose utilities
│   ├── io.R                 # load_seurat_object, read_csv_and_create_variables
│   ├── metadata.R           # process_seurat_metadata, transferMetadata,
│   │                        #   replace_surat_NA, update_seurat_metadata
│   ├── filter.R             # FilterSeuratForIntegration, splitSeuratObjectToList
│   └── consensus.R          # create_consensus_CT_clusters, create_consensus_Mps
│
└── cnv/                     # Copy-number variation analysis
    ├── infercnv.R           # updateSeuratMetadata (inferCNV reference setup)
    └── scevan.R             # plotCNA_withAnnotCells (SCEVAN heatmap)
```

---

## Setup

Install required R packages before using this toolkit:

```r
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork",
                   "futile.logger", "Cairo", "data.table",
                   "ggtext", "tibble", "tidyr", "RColorBrewer"))

# Bioconductor packages
BiocManager::install(c("scCustomize", "clustree", "infercnv"))

# GitHub packages
remotes::install_github("AntonioDeFalco/SCEVAN")
```

---

## Usage

Source the full toolkit with a single call:

```r
source("index.R")
```

Or source individual modules as needed:

```r
source("core/seurat_workflow.R")
source("plotting/save_plot.R")
```

### Quick Example

```r
source("index.R")

# Load and process a Seurat object
seu <- load_seurat_object("processed.rds", project_name = "MyProject", data_dir = "data/")
seu <- process_seurat_object(seu)
seu <- integrateAndClusterSeurat(seu, integrationMethod = HarmonyIntegration)

# Annotate clusters
seu <- name_seurat_cluster(seu, cluster_list = list(Neurons = c(0, 2), Microglia = c(1)))

# Save a UMAP plot
save_plot(DimPlot(seu), plot_dir = "plots/", file_name = "umap", format = "png")
```

---

## Module Descriptions

| Folder | Purpose |
|--------|---------|
| `core/` | Essential Seurat preprocessing, integration and QC filtering |
| `annotation/` | Cell-type prediction, cluster labelling and mapping dictionaries |
| `plotting/` | Plot generation and file output helpers |
| `data/` | Curated marker gene lists and color palettes |
| `utils/` | I/O, metadata manipulation, filtering and consensus calling |
| `cnv/` | Copy-number variation workflows (SCEVAN, inferCNV) |

# index.R
# Entry point: sources all modules in dependency order.
# Usage: source("index.R")

message("Loading scRNA analysis toolkit...")

# ── Core Seurat workflows ──────────────────────────────────────────────────────
source(file.path(dirname(sys.frame(1)$ofile), "core/seurat_workflow.R"))
source(file.path(dirname(sys.frame(1)$ofile), "core/integration.R"))
source(file.path(dirname(sys.frame(1)$ofile), "core/qc.R"))

# ── Annotation ────────────────────────────────────────────────────────────────
source(file.path(dirname(sys.frame(1)$ofile), "annotation/cluster_annotation.R"))
source(file.path(dirname(sys.frame(1)$ofile), "annotation/sctype.R"))
source(file.path(dirname(sys.frame(1)$ofile), "annotation/mappings.R"))

# ── Plotting ──────────────────────────────────────────────────────────────────
source(file.path(dirname(sys.frame(1)$ofile), "plotting/save_plot.R"))
source(file.path(dirname(sys.frame(1)$ofile), "plotting/stacked_bar.R"))
source(file.path(dirname(sys.frame(1)$ofile), "plotting/cluster_plots.R"))
source(file.path(dirname(sys.frame(1)$ofile), "plotting/gene_trends.R"))

# ── Reference data ────────────────────────────────────────────────────────────
source(file.path(dirname(sys.frame(1)$ofile), "data/canonical_markers.R"))
source(file.path(dirname(sys.frame(1)$ofile), "data/color_codes.R"))

# ── Utilities ─────────────────────────────────────────────────────────────────
source(file.path(dirname(sys.frame(1)$ofile), "utils/io.R"))
source(file.path(dirname(sys.frame(1)$ofile), "utils/metadata.R"))
source(file.path(dirname(sys.frame(1)$ofile), "utils/filter.R"))
source(file.path(dirname(sys.frame(1)$ofile), "utils/consensus.R"))

# ── CNV analysis ──────────────────────────────────────────────────────────────
source(file.path(dirname(sys.frame(1)$ofile), "cnv/infercnv.R"))
source(file.path(dirname(sys.frame(1)$ofile), "cnv/scevan.R"))

message("✅ All modules loaded.")

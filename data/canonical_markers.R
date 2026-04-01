# data/canonical_markers.R
# Reference marker gene lists for brain tumour cell type annotation.



## Versions 1 (23-04-2025)
pbta_marker_genes <- list(
  B_CELL = c("MS4A1", "CD79A", "CD79B", "CD19","IGKC","IGHG1"),
  T_CELL = c("CD3D", "CD3E", "CD3G", "IL7R"),
  Plasma = c("JCHAIN", "MZB1", "TNFRSF17"),
  Monocytes = c("CD14", "LYZ", "FCGR3A", "S100A9"),
  Macrophages = c("CD68", "CD163", "CSF1R", "MRC1"),
  Microglia = c("P2RY12", "TREM2", "CSF1R", "ITGAM", "CX3CR1", "SIGLEC1", "TMEM119"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "KDR"),
  Pericyte = c("RGS5", "ESAM", "MEF2C"),
  Fibroblast = c("COL1A1", "COL1A2", "PDGFRA", "DCN","LUM"),
  Neuronal = c("RBFOX3","MAP2", "SYN1", "ENO2", "CALB2", "DCX", "SMARCB1")
)
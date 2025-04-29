###
library(scCustomize)


create_color_palette_for_seurat <- function(seurat_object, metadata_column, palette = "ditto_seq", named_vector = TRUE) {
  # Extract unique values from the specified metadata column
  unique_values <- unique(seurat_object@meta.data[[metadata_column]])

  # Initialize an empty palette vector
  color_palette <- c()

  # Check if "Normal" is in the unique values and set its color to gray
  if ("Normal" %in% unique_values) {
    color_palette["Normal"] <- "#CCCCCC"
    unique_values <- setdiff(unique_values, "Normal") # Remove "Normal" from the list
  }

  # If there are any clone values, create a gradient for them
  if (length(unique_values) > 0) {
    # Generate a color palette using DiscretePalette_scCustomize
    gradient_colors <- DiscretePalette_scCustomize(num_colors = length(unique_values), palette = palette)
    color_palette <- c(color_palette, setNames(gradient_colors, unique_values))
  }

  # Return named vector if requested
  if (named_vector) {
    return(color_palette)
  } else {
    return(as.list(color_palette))
  }
}


set_annotation_colors <- function(annotation_column, palette = "varibow") {
  unique_levels <- unique(annotation_column)
  color_palette <- DiscretePalette_scCustomize(length(unique_levels), palette = palette)
  return(setNames(color_palette, unique_levels))
}


### Colors for Annotation

colors_cancer_type.V1 <- c(
  "ETMR" = "#9B59B6",
  "Medulloblastoma" = "#228658",
  
  "pLGG" = "#AD7700",
  "HGG" = "#ff7f00",
  "ATRT" = "#00969e",
  
  "CPC" = "#D5C711",
  "CPP" = "#0080ff",
  "aCPP" = "#ff618c",

  
  "PF-EPN" = "#ff88a9",
  "SP-EPN" = "#1C91D4",
  "ST-EPN" = "#64a649",
  "SUB-EPN" = "#0072B2"
)

colors_cancer_type.v2 <- c(
  "ETMR" = "#64a649",
  "Medulloblastoma" = "#228658",
  
  "pLGG" = "#5288ff",
  "HGG" = "#0080ff",
  "ATRT" = "#00969e",
  
  "CPC" = "#2a38ee",
  "CPP" = "#0f1bbd",
  
  "PF-EPN" = "#ff8d4d",
  "SP-EPN" = "#EE3B3B",
  "ST-EPN" = "#ff88a9",
  "SUB-EPN" = "#ff618c"
)



colors_qcs_consensus <- c(
  "Low nCOUNT" = "#33a02c",
  "Low Genes per Cell" = "#b2df8a",
  "High nCOUNT" = "#e31a1c",
  "High Genes per Cell" = "#fb9a99",
  "Doublet" = "#1f78b4",
  "Passed" = "#E5E5E5",
  "MiQC+" = "#ff7f00",
  "DoubletFinder+" = "#6a3d9a",
  "DblFinder+" = "#cab2d6"
)

colors_purity_category <- c(
  "Pure" = "#e31a1c",
  "Mixed" = "#ff7f00",
  "Impure" = "#66C2A5"
)

colors_qcs_phase <- c(
  "G1" = "#E5E5E5",
  "S" = "#44C2A5",
  "G2M" = "#e11a1c"
)

colors_pbta_annotation <- c(
  "Fibroblasts" = "#009E73",
  "CAFs" = "#548B54", 
  "Endothelial" = "#FF6961",
  "Immunes" = "#9B59B6",
  "HemoProgen" = "#F88DCE",
  "Macrophages/Monocytes" = "#56B4E9",
  "Microglial cells" = "#62B8D6",
  "T_cells" = "#DD79A7",
  "B_cells" = "#e11a1c",
  
  "Malignant" = "#ff7f00",
  "Neurons" = "#F0E442",
  "NeuroGlial" = "#E69F00",
  "Unknown" = "#EAEAE5",
  "Stem Cells Like" = "#D55E00"
)


colors_scATOMIC_pan_cancer <- c(
  "Cancer" = "#FF6961",
  "Normal" = "#E5E5E5"
)


colors_scATOMIC_copyKat <- c(
  aneuploid = "#FF6961",
  Unknown = "#E5E5E5",
  diploid = "#66C2A5"
)


colors_scTypeBrain <- c(
  # Fibroblasts
  "Fibroblast" = "#00A55C", # A green shade similar to Fibroblasts
  "Vascular-Related CAFs" = "#D55E00", # Orange similar to CAFs
  "Inflammatory/Immune-Related CAFs" = "#E74411", # A slightly darker orange than CAFs
  "Tumor-Associated/Matrix CAFs" = "#F85522", # Light orange-red similar to CAFs

  # T_cells
  "Naive CD8+ T cells" = "#C77D9F", # A soft pink-purple for T_cells
  "Naive CD4+ T cells" = "#D687B3", # A pastel pink for T_cells
  "Memory CD8+ T cells" = "#CE6A96", # A muted pinkish tone for T_cells
  "Memory CD4+ T cells" = "#F5A2C2", # A pale pink for T_cells
  "Effector CD8+ T cells" = "#D56D99", # A deeper pink for T_cells
  "Effector CD4+ T cells" = "#F88DCE", # A soft magenta for T_cells
  "γδ-T cells" = "#D67F99", # A medium pinkish purple for T_cells
  "CD8+ NKT-like cells" = "#D27A92", # A muted pinkish tone for T_cells
  "CD4+ NKT-like cells" = "#F29AC7", # A lighter magenta for T_cells

  # Macrophages/Monocytes
  "Classical Monocytes" = "#6497B1", # A muted blue for Macrophages
  "Non-classical monocytes" = "#8EACC6", # A softer blue for Macrophages
  "Intermediate monocytes" = "#A4C4D9", # A pastel blue for Macrophages
  "Macrophages" = "#6DA8D6", # A light blue for Macrophages

  # Endothelial
  "Endothelial cells" = "#FF7D6C", # A muted red-orange for Endothelial
  "Endothelial" = "#FF6A53", # A deeper red for Endothelial
  "Pericyte" = "#D95F4B", # A darker red for Endothelial

  # B_cells
  "Pro-B cells" = "#D7261E", # A deep red for B_cells
  "Pre-B cells" = "#BF1A1A", # A bright red for B_cells
  "Naive B cells" = "#E1422F", # A soft red for B_cells
  "Memory B cells" = "#B02D2A", # A darker red for B_cells
  "Plasma B cells" = "#C94C3A", # A medium red for B_cells

  # Neurons
  "Astrocytes" = "#F0A500", # A warm yellow for Neurons
  "Cholinergic neurons" = "#E8A500", # A more golden yellow for Neurons
  "Dopaminergic neurons" = "#F1B700", # A richer yellow for Neurons
  "GABAergic neurons" = "#F2C800", # A soft yellow-orange for Neurons
  "Glutamatergic neurons" = "#F8D200", # A brighter yellow for Neurons
  "Immature neurons" = "#F7D900", # A fresh yellow for Neurons
  "Mature neurons" = "#F8E400", # A pale yellow for Neurons
  "Neurons" = "#F9E600", # A light yellow for Neurons
  "Neuron" = "#F9E600", # A light yellow for Neurons
  "Serotonergic neurons" = "#F3D400", # A deeper yellow for Neurons

  # Malignant
  "Cancer cells" = "#E0C900", # A golden yellow for Malignant
  "Cancer stem cells" = "#D5B600", # A mustard yellow for Malignant

  # Immunes
  "Immune system cells" = "#9B59B6", # A purple for Immunes
  "ISG expressing immune cells" = "#9D4C97", # A deep purple for Immunes
  "Natural killer cells" = "#B35A9C", # A medium purple for Immunes
  "Eosinophils" = "#A25591", # A muted purple for Immunes
  "Neutrophils" = "#AD4A88", # A soft purple for Immunes
  "Basophils" = "#B66F9C", # A soft lilac for Immunes
  "Mast cells" = "#9A4A83", # A deep lilac for Immunes
  "Platelets" = "#B383B0", # A lavender for Immunes
  "Myeloid Dendritic cells" = "#8B4B95", # A dark purple for Immunes
  "Plasmacytoid Dendritic cells" = "#9F57A4", # A rich purple for Immunes
  "Granulocytes" = "#9D56A2", # A muted purple for Immunes
  "Microglial cells" = "#62B8D6",

  # Neural & Glial Cells (NeuroGlial)
   # A warm yellow-orange for NeuroGlial
  "Myelinating Schwann cells" = "#E09E00", # A golden yellow for NeuroGlial
  "Neural Progenitor cells" = "#D59F00", # A deep golden yellow for NeuroGlial
  "Neural Stem Cells" = "#866500", # A mustard yellow for NeuroGlial
  "Neuroblasts" = "#C7A300", # A darker yellow-orange for NeuroGlial
  "Neuroepithelial" = "#D0A400", # A soft orange-yellow for NeuroGlial
  "Non myelinating Schwann" = "#E0A400", # A warm golden yellow for NeuroGlial
  "Oligodendrocyte precursor" = "#E99F00", # A vibrant yellow-orange for NeuroGlial
  "Oligodendrocytes" = "#F1B500", # A bright golden-yellow for NeuroGlial
  "Radial glial cells" = "#F2B800", # A rich yellow for NeuroGlial
  "Schwann precursor cells" = "#F4B200", # A deep golden-yellow for NeuroGlial
  "Tanycytes" = "#D99E00", # A soft yellow for NeuroGlial
  "Oligodendrocyte" = "#E9A600", # A warm yellow for NeuroGlial
  "OPC" = "#E6A800", # A deeper yellow for NeuroGlial

  # Hematopoietic & Progenitor Cells (HemoProgen)
  "Progenitor cells" = "#F7A400", # A vibrant yellow for HemoProgen
  "Erythroid precursor cells" = "#F8A600", # A golden yellow for HemoProgen
  "HSC/MPP cells" = "#E29F00", # A dark yellow for HemoProgen

  "Unknown" = "#EAEAE5" # A neutral gray for Unknown
)


colors_panglaoBrain <- c(
  # Neural & Glial Cells (NeuroGlial)
  "Chondrocytes" = "#00A55C", # A greenish tone for NeuroGlial
  "Stromal cells" = "#1D7D56", # A dark green for NeuroGlial
  "Adrenergic neurons" = "#FF9F5B", # An orange for NeuroGlial
  "Astrocytes" = "#f16000", # A golden yellow for NeuroGlial
  "Bergmann glia" = "#D49F00", # A warm golden yellow for NeuroGlial
  "Cajal-Retzius cells" = "#FFCC00", # A yellow for NeuroGlial
  "Cholinergic neurons" = "#F8D200", # A bright yellow for NeuroGlial
  "Choroid plexus cells" = "#be851c", # A soft yellow-orange for NeuroGlial
  "Dopaminergic neurons" = "#F1B700", # A rich yellow for Neurons
  "Ependymal cells" = "#D0A400", # A yellow-orange for NeuroGlial
  "GABAergic neurons" = "#F7D100", # A deep yellow for Neurons
  "Glutaminergic neurons" = "#F8E100", # A vibrant yellow for Neurons
  "Immature neurons" = "#F7D900", # A fresh yellow for NeuroGlial
  "Interneurons" = "#F9E600", # A light yellow for Neurons
  "Meningeal cells" = "#D99E00", # A soft yellow for NeuroGlial
  
  "Motor neurons" = "#F3D400", # A deeper yellow for Neurons
  "Neural stem/precursor cells" = "#E0A500", # A vibrant yellow-orange for NeuroGlial
  "Neuroblasts" = "#c75a00", # A darker yellow-orange for NeuroGlial
  "Neuroendocrine cells" = "#F5C400", # A golden yellow for NeuroGlial
  "Neuroepithelial cells" = "#d1990a", # A soft golden-yellow for NeuroGlial
  "Neurons" = "#dff138", # A light yellow for Neurons
  "Noradrenergic neurons" = "#F6E100", # A pale yellow for Neurons
  "Oligodendrocyte progenitor cells" = "#8adf2a", # A warm yellow for NeuroGlial
  "Oligodendrocytes" = "#83672c", # A bright yellow for NeuroGlial
  "Pyramidal cells" = "#F2B800", # A rich yellow for Neurons
  "Radial glia cells" = "#c02424", # A deep golden-yellow for NeuroGlial
  "Retinal ganglion cells" = "#F8D200", # A golden yellow for Neurons
  "Satellite glial cells" = "#D99E00", # A soft yellow for NeuroGlial
  "Schwann cells" = "#F3B500", # A soft yellow-orange for NeuroGlial
  "Serotonergic neurons" = "#F2D700", # A yellow-orange for Neurons
  "Purkinje neurons" = "#F7D500", # A rich yellow for Neurons
  "Trigeminal neurons" = "#F9D800", # A golden-yellow for Neurons
  "Chromaffin cells" = "#E7A200", # A warm yellow-orange for NeuroGlial

  # Immune Cells (Immunes)
  "B cells" = "#9B59B6", # A purple for B_cells
  "B cells memory" = "#7d3e80", # A deeper purple for B_cells
  "B cells naive" = "#e7199f", # A soft purple for B_cells
  "Basophils" = "#8D6E92", # A muted purple for Immunes
  "Dendritic cells" = "#9B58B6", # A rich purple for Immunes
  "Eosinophils" = "#9A4F88", # A deep purple for Immunes
  "Gamma delta T cells" = "#8E44AD", # A purple for T_cells
  "Macrophages" = "#56B4E9", # A light blue for Macrophages/Monocytes
  "Mast cells" = "#9B59B6", # A soft purple for Immunes
  "Megakaryocytes" = "#7B4F9B", # A dark purple for Immunes
  "Monocytes" = "#3D8DAE", # A muted blue for Macrophages/Monocytes
  "Myeloid-derived suppressor cells" = "#A59ACF", # A pale blue for Immunes
  "NK T cells" = "#9B59B6", # A deep purple for T_cells
  "Neutrophils" = "#A6A6FF", # A light blue for Immunes
  "NK cells" = "#8E44AD", # A purple for Immunes
  "Plasma cells" = "#C94C3A", # A medium red for B_cells
  "Plasmacytoid dendritic cells" = "#9F57A4", # A rich purple for Immunes
  "Platelets" = "#A1A1A1", # A soft grey for Immunes
  "T cells" = "#DD79A7", # A soft pink for T_cells
  "T cytotoxic cells" = "#D56D99", # A deeper pink for T_cells
  "T follicular helper cells" = "#F88DCE", # A soft magenta for T_cells
  "T helper cells" = "#F7C1D7", # A pale pink for T_cells
  "T memory cells" = "#E59BCF", # A medium pink for T_cells
  "T regulatory cells" = "#D37F8D", # A muted pink for T_cells
  "Microglia" = "#62B8D6", # A warm yellow-orange for NeuroGlial

  # Epithelial and Muscle Cells (Miscellaneous)
  "Basal cells" = "#2D6A4F", # A darker green for Fibroblasts
  "Epithelial cells" = "#4A9D6F", # A softer green for Fibroblasts
  "Mesothelial cells" = "#5A7E5C", # A muted green for Fibroblasts
  "Airway smooth muscle cells" = "#6D9E79", # A medium green for Fibroblasts
  "Myoepithelial cells" = "#3B7A53", # A deep green for Fibroblasts
  "Smooth Muscle cells" = "#6B8C6A", # A soft green for Fibroblasts
  "Fibroblasts" = "#009E73", # A green similar to Fibroblasts

  # Stem Cells (Pluripotent)
  "Embryonic stem cells" = "#F8D100", # A bright yellow for HemoProgen
  "Epiblast cells" = "#F9D400", # A golden yellow for HemoProgen
  "Germ cells" = "#F6C300", # A rich yellow for HemoProgen
  "Pluripotent stem cells" = "#F7D200", # A warm yellow for HemoProgen

  # Endothelial & Pericytes (Endothelial)
  "Endothelial cells" = "#FF6A53", # A deeper red for Endothelial
  "Endothelial cells (blood brain barrier)" = "#FF7D6C", # A muted red-orange for Endothelial
  "Pericytes" = "#D95F4B", # A darker red for Endothelial

  # Hematopoietic Cells (HemoProgen)
  "Erythroblasts" = "#E1A800", # A warm yellow for HemoProgen
  "Erythroid-like and erythroid precursor cells" = "#D89F00", # A deep yellow for HemoProgen
  "Hematopoietic stem cells" = "#F0A400", # A bright yellow for HemoProgen

  "Unknown" = "#EAEAE5" # A neutral gray for Unknown
)


colors_mapping_Tirosh_MP <- c(
  # Malignant Group Mapping
  "Malignant:Cell Cycle - G2/M" = "#F4A300", # A vibrant yellow for Malignant
  "Malignant:Cell Cycle - G1/S" = "#F6B800", # A bright yellow-orange for Malignant
  "Malignant:Cell Cycle HMG-rich" = "#F1B400", # A soft yellow for Malignant
  "Malignant:Chromatin" = "#F2C700", # A deeper yellow for Malignant
  "Malignant:Stress" = "#FF6A53", # A soft red-orange for Malignant
  "Malignant:Hypoxia" = "#FF8F65", # A warm red for Malignant
  "Malignant:Stress (in vitro)" = "#E64A19", # A vibrant orange-red for Malignant
  "Malignant:Proteasomal degradation" = "#D44F00", # A deeper orange for Malignant
  "Malignant:Unfolded protein response" = "#FF9F5B", # A golden-orange for Malignant
  "Malignant:Protein maturation" = "#FFB74D", # A soft orange for Malignant
  "Malignant:Translation initiation" = "#F2A51D", # A yellow-orange for Malignant
  "Malignant:EMT-I" = "#F4C400", # A deep yellow for Malignant
  "Malignant:EMT-II" = "#F7D300", # A brighter yellow for Malignant
  "Malignant:EMT-III" = "#F9E100", # A pale yellow for Malignant
  "Malignant:EMT-IV" = "#FFCC00", # A soft yellow for Malignant
  "Malignant:MES (glioma)" = "#F5A600", # A warm golden yellow for Malignant
  "Malignant:Interferon/MHC-II (I)" = "#F8C500", # A rich yellow-orange for Malignant
  "Malignant:Interferon/MHC-II (II)" = "#F3C300", # A darker yellow for Malignant
  "Malignant:Epithelial Senescence" = "#F7D400", # A soft yellow for Malignant
  "Malignant:MYC" = "#E0A400", # A golden yellow for Malignant
  "Malignant:Respiration" = "#FFB400", # A golden-orange for Malignant
  "Malignant:Secreted I" = "#F1D500", # A rich yellow for Malignant
  "Malignant:Secreted II" = "#FFDB5C", # A light golden yellow for Malignant
  "Malignant:Cilia" = "#D98F00", # A darker yellow for Malignant
  "Malignant:Astrocytes" = "#E79A00", # A warm yellow for Malignant
  "Malignant:NPC Glioma" = "#FF8800", # A soft orange for Malignant
  "Malignant:Oligo Progenitor" = "#FF9600", # A yellow-orange for Malignant
  "Malignant:Oligo normal" = "#F58500", # A rich orange for Malignant
  "Malignant:NPC/OPC" = "#D97300", # A deeper orange for Malignant
  "Malignant:PDAC-classical" = "#F8B200", # A golden yellow for Malignant
  "Malignant:Alveolar" = "#F9C200", # A soft yellow for Malignant
  "Malignant:Skin-pigmentation" = "#F9D400", # A light yellow for Malignant
  "Malignant:RBCs" = "#D85F1A", # A reddish orange for Malignant
  "Malignant:Platelet-activation" = "#F2B500", # A rich yellow for Malignant
  "Malignant:Hemato-related-I" = "#FF6A00", # A deep orange for Malignant
  "Malignant:IG" = "#D95F00", # A dark orange for Malignant
  "Malignant:Hemato-related-II" = "#F38B00", # A golden orange for Malignant
  "Malignant:Glutathione" = "#F1A100", # A soft golden yellow for Malignant
  "Malignant:Metal-response" = "#FF8B00", # A deeper orange for Malignant
  "Malignant:PDAC-related" = "#F0A800", # A rich yellow-orange for Malignant
  "Malignant:Unassigned" = "#D8B100", # A pale yellow for Malignant
  
  # B cells
  "B cells:Plasma" = "#9B59B6", # A purple for B_cells
  "B cells:MHC-II" = "#8E44AD", # A deep purple for B_cells
  "B cells:Cell Cycle" = "#A56D91", # A muted purple for B_cells
  "B cells:Stress" = "#9B58B6", # A rich purple for B_cells
  "B cells:Memory" = "#A36B95", # A soft purple for B_cells
  "B cells:Metabolism/MYC" = "#A59BCF", # A medium purple for B_cells
  "B cells:Germinal Center" = "#9D4C97", # A deep purple for B_cells
  "B cells:Interferon" = "#D37F8D", # A muted pink for B_cells
  "B cells:Progenitor" = "#9A5B94", # A deep magenta for B_cells
  "B cells:B-cells1" = "#C77D9F", # A light purple for B_cells
  "B cells:Respiration" = "#9E6F92", # A soft pink for B_cells
  "B cells:HSP/Stress" = "#9D4C97", # A rich purple for B_cells
  
  # CD4 T cells
  "CD4 T cells:T_reg" = "#DD79A7", # A soft pink for T_cells
  "CD4 T cells:Naive1" = "#D56D99", # A muted pink for T_cells
  "CD4 T cells:Cell Cycle" = "#F88DCE", # A bright pink for T_cells
  "CD4 T cells:Cytotoxic" = "#F7C1D7", # A soft pink for T_cells
  "CD4 T cells:Dysfunction" = "#F7D300", # A deeper pink for T_cells
  "CD4 T cells:Interferon" = "#F29AC7", # A soft magenta for T_cells
  "CD4 T cells:Glycolysis/MYC" = "#F7C1D7", # A light pink for T_cells
  "CD4 T cells:Naive2" = "#F29AC7", # A soft magenta for T_cells
  "CD4 T cells:Unassigned" = "#F7D300", # A soft pink for T_cells
  "CD4 T cells:Stress/HSP" = "#E5A3A1", # A light pink for T_cells
  
  # CD8 T cells
  "CD8 T cells:Cytotoxic" = "#F56B8A", # A red for T_cells
  "CD8 T cells:Cell Cycle" = "#F8A5B8", # A soft pink for T_cells
  "CD8 T cells:Memory/Naive1" = "#F0A9C5", # A light pink for T_cells
  "CD8 T cells:Interferon" = "#F29AC7", # A soft magenta for T_cells
  "CD8 T cells:Unassigned1" = "#F1B1C8", # A pastel pink for T_cells
  "CD8 T cells:Naive2" = "#D56D99", # A muted pink for T_cells
  "CD8 T cells:Glycolysis/MYC" = "#F7D300", # A soft yellow for T_cells
  "CD8 T cells:Chromatin" = "#D87D97", # A deep pink for T_cells
  "CD8 T cells:Unassigned2" = "#F6C1D9", # A soft pink for T_cells
  "CD8 T cells:Heat_shock" = "#F7D200", # A rich yellow for T_cells
  "CD8 T cells:Naive3" = "#F7C1D7", # A light pink for T_cells
  "CD8 T cells:Dysfunction" = "#F56A89", # A dark red for T_cells
  
  # Endothelial
  "Endothelial:Notch-signaling" = "#E74C3C", # A soft red for Endothelial
  "Endothelial:HEV1" = "#D9534F", # A muted red for Endothelial
  "Endothelial:HEV2" = "#C0392B", # A deeper red for Endothelial
  "Endothelial:Endo1" = "#FF6A53", # A warm red-orange for Endothelial
  "Endothelial:Endo2" = "#D95F4B", # A muted red for Endothelial
  "Endothelial:Endo3" = "#FF9F7F", # A soft red-orange for Endothelial
  "Endothelial:Endo4" = "#C0392B", # A darker red for Endothelial
  "Endothelial:Endo5" = "#F25D3D", # A warm red for Endothelial
  "Endothelial:Endo6" = "#FF7F50", # A light red-orange for Endothelial
  "Endothelial:Endo7" = "#FF6347", # A soft tomato red for Endothelial
  "Endothelial:Stress" = "#FF4500", # A deep orange-red for Endothelial
  "Endothelial:Cell Cycle" = "#F44336", # A vivid red for Endothelial
  "Endothelial:Interferon" = "#E74C3C", # A red for Endothelial
  "Endothelial:NF-kb" = "#FF5733", # A warm red-orange for Endothelial
  "Endothelial:Coagulation" = "#F44336", # A strong red for Endothelial
  
  # Macrophages
  "Macrophages:Lipid-associated" = "#56B4E9", # A light blue for Macrophages
  "Macrophages:Monocyte/Secreted" = "#3D8DAE", # A muted blue for Macrophages
  "Macrophages:Cell-cycle" = "#A6D1F7", # A soft blue for Macrophages
  "Macrophages:Interferon" = "#A59BCF", # A light purple for Macrophages
  "Macrophages:MES/Glycolysis" = "#A2C1E3", # A soft blue for Macrophages
  "Macrophages:MAC1" = "#81A7D7", # A cool blue for Macrophages
  "Macrophages:MAC2" = "#90B9D9", # A light blue for Macrophages
  "Macrophages:MAC3" = "#70A2D5", # A muted blue for Macrophages
  "Macrophages:Stress/HSP" = "#D0D9F0", # A pastel blue for Macrophages
  "Macrophages:Proteasomal-degradation" = "#A1B1D4", # A soft blue for Macrophages
  "Macrophages:MYC/Mitochondria" = "#88B0C8", # A soft blue for Macrophages
  "Macrophages:Unfolded-protein-response" = "#A0B1C9", # A cool blue for Macrophages
  "Macrophages:Respiration" = "#80A0C6", # A deeper blue for Macrophages
  
  # Fibroblasts
  "Fibroblasts :CAF1" = "#009E73", # A green for Fibroblasts
  "Fibroblasts :CAF2" = "#1F7A4A", # A darker green for Fibroblasts
  "Fibroblasts :CAF3" = "#56B4E9", # A soft blue for Fibroblasts
  "Fibroblasts :CAF4" = "#8D6E5B", # A warm brown for Fibroblasts
  "Fibroblasts :CAF5" = "#B48B3A", # A golden yellow for Fibroblasts
  "Fibroblasts :CAF6" = "#D09E00", # A mustard yellow for Fibroblasts
  "Fibroblasts :CAF7" = "#B7C54C", # A yellow-green for Fibroblasts
  "Fibroblasts :CAF8" = "#93A62C", # A muted yellow for Fibroblasts
  "Fibroblasts :CAF9" = "#6C8A3C", # A greenish yellow for Fibroblasts
  "Fibroblasts :CAF10" = "#A0C25C", # A soft yellow-green for Fibroblasts
  "Fibroblasts :Cell-cycle" = "#88B14D", # A fresh green for Fibroblasts
  "Fibroblasts :Stress" = "#A7D65C", # A fresh yellow-green for Fibroblasts
  "Fibroblasts :Hypoxia" = "#99C82F", # A yellow-green for Fibroblasts
  "Fibroblasts :Complement" = "#8D9A31", # A muted yellow for Fibroblasts
  "Fibroblasts :Interferon" = "#6B7A28", # A deeper yellow for Fibroblasts
  "Fibroblasts :Myofibroblasts" = "#70852E", # A darker yellow-green for Fibroblasts
  "Fibroblasts :PI16+" = "#4F7741", # A muted green for Fibroblasts
  "Fibroblasts :Metal-response" = "#7C8A30", # A rich yellow-green for Fibroblasts
  "Fibroblasts :MHC-II/Cytokine" = "#A2BB2B", # A bright yellow-green for Fibroblasts
  "Fibroblasts :Lipid-metabolism" = "#91A324", # A fresh yellow-green for Fibroblasts
  "Fibroblasts :MHC-II" = "#6E7F1E", # A dark yellow-green for Fibroblasts
  "Fibroblasts :Pericyte-like" = "#5E6F1A", # A muted yellow-green for Fibroblasts
  
  "Unknown" = "#EAEAE5" # A neutral gray for Unknown
)

colors_singleR_cluster <- c(
  # NeuroGlial
  "Neuroepithelial_cell" = "#F4A500", 
  "Astrocyte" = "#E79A00", 
  
  # Neurons
  "Neurons" = "#F9E600", # A light yellow for Neurons
  
  # Endothelial
  "Endothelial_cells" = "#FF7D6C", # A warm red-orange for Endothelial
  
  # Macrophages/Monocytes
  "Macrophage" = "#56B4E9", # A soft blue for Macrophages/Monocytes
  "Monocyte" = "#6497B1",
  
  # Tissue and Stem Cells
  "Tissue_stem_cells" = "#009E73",
  "Smooth_muscle_cells" = "#56B4E9", 
  
  # MSC (Mesenchymal Stem Cells)
  "MSC" = "#009E73", # A green for Fibroblasts (MSC-like cells)
  
  # T_cells
  "T_cells" = "#DD79A7", # A soft pink for T_cells
  
  # B_cells
  "B_cell" = "#9B59B6", # A purple for B_cells
  
  # Immunes
  "iPS_cells" = "#9B59B6",
  "Neutrophils" = "#A6A6FF",
  "Unknown" = "#EAEAE5"
)



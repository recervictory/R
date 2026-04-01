# annotation/mappings.R
# Cell type mapping dictionaries for harmonizing annotation labels
# across different annotation tools (SingleR, scType, PanglaoDB, Tirosh MPs).
mapping_singleR_cluster <- c(
  # NeuroGlial
  "Neuroepithelial_cell" = "NeuroGlial",
  "Astrocyte" = "NeuroGlial",

  # Neurons
  "Neurons" = "Neurons",

  # Endothelial
  "Endothelial_cells" = "Endothelial",

  # Macrophages/Monocytes
  "Macrophage" = "Macrophages/Monocytes",
  "Monocyte" = "Macrophages/Monocytes", # Myeloid


  # Tissue and Stem Cells
  "Tissue_stem_cells" = "Fibroblasts", # Assuming these cells function as fibroblast-like
  "Smooth_muscle_cells" = "Fibroblasts", # These could be grouped under fibroblasts due to their structural role

  # MSC (Mesenchymal Stem Cells)
  "MSC" = "Fibroblasts", # MSCs can be grouped under fibroblasts or related stromal cells

  # T_cells
  "T_cells" = "T_cells",

  # B_cells
  "B_cell" = "B_cells",

  # Immunes
  "iPS_cells" = "Immunes",
  "Neutrophils" = "Immunes",
  "Unknown" = "Unknown"
)



mapping_scTypeBrain <- c(
  # Fibroblasts
  "Fibroblast" = "Fibroblasts",
  "Vascular-Related CAFs" = "CAFs",
  "Inflammatory/Immune-Related CAFs" = "CAFs",
  "Tumor-Associated/Matrix CAFs" = "CAFs",
  "Pericyte" = "Fibroblasts",

  # T_cells
  "Naive CD8+ T cells" = "T_cells",
  "Naive CD4+ T cells" = "T_cells",
  "Memory CD8+ T cells" = "T_cells",
  "Memory CD4+ T cells" = "T_cells",
  "Effector CD8+ T cells" = "T_cells",
  "Effector CD4+ T cells" = "T_cells",
  "γδ-T cells" = "T_cells",
  "CD8+ NKT-like cells" = "T_cells",
  "CD4+ NKT-like cells" = "T_cells",

  # Macrophages/Monocytes
  "Classical Monocytes" = "Macrophages/Monocytes",
  "Non-classical monocytes" = "Macrophages/Monocytes",
  "Intermediate monocytes" = "Macrophages/Monocytes",
  "Macrophages" = "Macrophages/Monocytes",

  # Endothelial
  "Endothelial cells" = "Endothelial",
  "Endothelial" = "Endothelial",


  # B_cells
  "Pro-B cells" = "B_cells",
  "Pre-B cells" = "B_cells",
  "Naive B cells" = "B_cells",
  "Memory B cells" = "B_cells",
  "Plasma B cells" = "B_cells",

  # Neurons
  "Astrocytes" = "Neurons",
  "Cholinergic neurons" = "Neurons",
  "Dopaminergic neurons" = "Neurons",
  "GABAergic neurons" = "Neurons",
  "Glutamatergic neurons" = "Neurons",
  "Immature neurons" = "Neurons",
  "Mature neurons" = "Neurons",
  "Neurons" = "Neurons",
  "Neuron" = "Neurons",
  "Serotonergic neurons" = "Neurons",

  # Malignant
  "Cancer cells" = "Malignant",
  "Cancer stem cells" = "Malignant",


  # Immunes
  "Immune system cells" = "Immunes",
  "ISG expressing immune cells" = "Immunes",
  "Natural killer cells" = "Immunes",
  "Eosinophils" = "Immunes",
  "Neutrophils" = "Immunes",
  "Basophils" = "Immunes",
  "Mast cells" = "Immunes",
  "Platelets" = "Immunes",
  "Myeloid Dendritic cells" = "Immunes",
  "Plasmacytoid Dendritic cells" = "Immunes",
  "Granulocytes" = "Immunes",
  "Microglial cells" = "Microglial cells",

  # Neural & Glial Cells (NeuroGlial)
  "Myelinating Schwann cells" = "NeuroGlial",
  "Neural Progenitor cells" = "NeuroGlial",
  "Neural Stem Cells" = "NeuroGlial",
  "Neuroblasts" = "NeuroGlial",
  "Neuroepithelial" = "NeuroGlial",
  "Non myelinating Schwann" = "NeuroGlial",
  "Oligodendrocyte precursor" = "NeuroGlial",
  "Oligodendrocytes" = "NeuroGlial",
  "Radial glial cells" = "NeuroGlial",
  "Schwann precursor cells" = "NeuroGlial",
  "Tanycytes" = "NeuroGlial",
  "Oligodendrocyte" = "NeuroGlial",
  "OPC" = "NeuroGlial",

  # Hematopoietic & Progenitor Cells (HemoProgen)
  "Progenitor cells" = "HemoProgen",
  "Erythroid precursor cells" = "HemoProgen",
  "HSC/MPP cells" = "HemoProgen",
  "Unknown" = "Unknown"
)



mapping_panglaoBrain <- c(
  # Neural & Glial Cells (NeuroGlial)
  "Chondrocytes" = "NeuroGlial",
  "Stromal cells" = "NeuroGlial",
  "Adrenergic neurons" = "NeuroGlial",
  "Astrocytes" = "NeuroGlial",
  "Bergmann glia" = "NeuroGlial",
  "Cajal-Retzius cells" = "NeuroGlial",
  "Cholinergic neurons" = "NeuroGlial",
  "Choroid plexus cells" = "NeuroGlial",
  "Immature neurons" = "NeuroGlial",
  "Ependymal cells" = "NeuroGlial",
  "Neural stem/precursor cells" = "NeuroGlial",
  "Oligodendrocyte progenitor cells" = "NeuroGlial",
  "Meningeal cells" = "NeuroGlial",

  "Fibroblasts" = "Fibroblasts",
  "Pericytes" = "Fibroblasts",

  
  "Dopaminergic neurons" = "Neurons",

  "GABAergic neurons" = "Neurons",
  "Glutaminergic neurons" = "Neurons",
  
  "Interneurons" = "Neurons",
  
  "Motor neurons" = "Neurons",
  
  
  "Neurons" = "Neurons",
  "Neuron" = "Neurons",
  "Noradrenergic neurons" = "Neurons",
  
  "Oligodendrocytes" = "NeuroGlial",
  "Pyramidal cells" = "Neurons",
  "Radial glia cells" = "NeuroGlial",
  "Retinal ganglion cells" = "Neurons",
  "Satellite glial cells" = "NeuroGlial",
  "Schwann cells" = "NeuroGlial",
  "Serotonergic neurons" = "Neurons",
  "Purkinje neurons" = "Neurons",
  "Trigeminal neurons" = "Neurons",
  "Chromaffin cells" = "NeuroGlial",

  # Immune Cells (Immunes)
  "B cells" = "B_cells",
  "B cells memory" = "B_cells",
  "B cells naive" = "B_cells",
  "Basophils" = "Immunes",
  "Dendritic cells" = "Immunes",
  "Eosinophils" = "Immunes",
  "Gamma delta T cells" = "T_cells",
  "Macrophages" = "Macrophages/Monocytes",
  "Mast cells" = "Immunes",
  "Megakaryocytes" = "Immunes",
  "Monocytes" = "Macrophages/Monocytes",
  "Myeloid-derived suppressor cells" = "Immunes",
  "NK T cells" = "T_cells",
  "Neutrophils" = "Immunes",
  "NK cells" = "Immunes",
  "Plasma cells" = "B_cells",
  "Plasmacytoid dendritic cells" = "Immunes",
  "Platelets" = "Immunes",
  "T cells" = "T_cells",
  "T cytotoxic cells" = "T_cells",
  "T follicular helper cells" = "T_cells",
  "T helper cells" = "T_cells",
  "T memory cells" = "T_cells",
  "T regulatory cells" = "T_cells",
  "Microglia" = "Microglial cells",

  # Epithelial and Muscle Cells (Miscellaneous)
  "Basal cells" = "Fibroblasts",
  "Epithelial cells" = "Fibroblasts",
  "Mesothelial cells" = "Fibroblasts",
  "Airway smooth muscle cells" = "Fibroblasts",
  "Myoepithelial cells" = "Fibroblasts",
  "Smooth Muscle cells" = "Fibroblasts",

  # Stem Cells (Pluripotent)
  "Embryonic stem cells" = "Stem Cells Like",
  "Epiblast cells" = "Stem Cells Like",
  "Germ cells" = "Stem Cells Like",
  "Pluripotent stem cells" = "Stem Cells Like",

  # Endothelial & Pericytes (Endothelial)
  "Endothelial cells" = "Endothelial",
  "Endothelial cells (blood brain barrier)" = "Endothelial",


  # Hematopoietic Cells (HemoProgen)
  "Erythroblasts" = "HemoProgen",
  "Erythroid-like and erythroid precursor cells" = "HemoProgen",
  "Hematopoietic stem cells" = "HemoProgen",
  "Unknown" = "Unknown"
)


mapping_Tirosh_MP <- c(
  # Malignant Group Mapping
  "Malignant:Cell Cycle - G2/M" = "Malignant",
  "Malignant:Cell Cycle - G1/S" = "Malignant",
  "Malignant:Cell Cycle HMG-rich" = "Malignant",
  "Malignant:Chromatin" = "Malignant",
  "Malignant:Stress" = "Malignant",
  "Malignant:Hypoxia" = "Malignant",
  "Malignant:Stress (in vitro)" = "Malignant",
  "Malignant:Proteasomal degradation" = "Malignant",
  "Malignant:Unfolded protein response" = "Malignant",
  "Malignant:Protein maturation" = "Malignant",
  "Malignant:Translation initiation" = "Malignant",
  "Malignant:EMT-I" = "Malignant",
  "Malignant:EMT-II" = "Malignant",
  "Malignant:EMT-III" = "Malignant",
  "Malignant:EMT-IV" = "Malignant",
  "Malignant:MES (glioma)" = "Malignant",
  "Malignant:Interferon/MHC-II (I)" = "Malignant",
  "Malignant:Interferon/MHC-II (II)" = "Malignant",
  "Malignant:Epithelial Senescence" = "Malignant",
  "Malignant:MYC" = "Malignant",
  "Malignant:Respiration" = "Malignant",
  "Malignant:Secreted I" = "Malignant",
  "Malignant:Secreted II" = "Malignant",
  "Malignant:Cilia" = "Malignant",
  "Malignant:Astrocytes" = "Malignant",
  "Malignant:NPC Glioma" = "Malignant",
  "Malignant:Oligo Progenitor" = "Malignant",
  "Malignant:Oligo normal" = "Malignant",
  "Malignant:NPC/OPC" = "Malignant",
  "Malignant:PDAC-classical" = "Malignant",
  "Malignant:Alveolar" = "Malignant",
  "Malignant:Skin-pigmentation" = "Malignant",
  "Malignant:RBCs" = "Malignant",
  "Malignant:Platelet-activation" = "Malignant",
  "Malignant:Hemato-related-I" = "Malignant",
  "Malignant:IG" = "Malignant",
  "Malignant:Hemato-related-II" = "Malignant",
  "Malignant:Glutathione" = "Malignant",
  "Malignant:Metal-response" = "Malignant",
  "Malignant:PDAC-related" = "Malignant",
  "Malignant:Unassigned" = "Malignant",

  # B cells
  "B cells:Plasma" = "B_cells",
  "B cells:MHC-II" = "B_cells",
  "B cells:Cell Cycle" = "B_cells",
  "B cells:Stress" = "B_cells",
  "B cells:Memory" = "B_cells",
  "B cells:Metabolism/MYC" = "B_cells",
  "B cells:Germinal Center" = "B_cells",
  "B cells:Interferon" = "B_cells",
  "B cells:Progenitor" = "B_cells",
  "B cells:B-cells1" = "B_cells",
  "B cells:Respiration" = "B_cells",
  "B cells:HSP/Stress" = "B_cells",

  # CD4 T cells
  "CD4 T cells:T_reg" = "T_cells",
  "CD4 T cells:Naive1" = "T_cells",
  "CD4 T cells:Cell Cycle" = "T_cells",
  "CD4 T cells:Cytotoxic" = "T_cells",
  "CD4 T cells:Dysfunction" = "T_cells",
  "CD4 T cells:Interferon" = "T_cells",
  "CD4 T cells:Glycolysis/MYC" = "T_cells",
  "CD4 T cells:Naive2" = "T_cells",
  "CD4 T cells:Unassigned" = "T_cells",
  "CD4 T cells:Stress/HSP" = "T_cells",

  # CD8 T cells
  "CD8 T cells:Cytotoxic" = "T_cells",
  "CD8 T cells:Cell Cycle" = "T_cells",
  "CD8 T cells:Memory/Naive1" = "T_cells",
  "CD8 T cells:Interferon" = "T_cells",
  "CD8 T cells:Unassigned1" = "T_cells",
  "CD8 T cells:Naive2" = "T_cells",
  "CD8 T cells:Glycolysis/MYC" = "T_cells",
  "CD8 T cells:Chromatin" = "T_cells",
  "CD8 T cells:Unassigned2" = "T_cells",
  "CD8 T cells:Heat_shock" = "T_cells",
  "CD8 T cells:Naive3" = "T_cells",
  "CD8 T cells:Dysfunction" = "T_cells",

  # Endothelial
  "Endothelial:Notch-signaling" = "Endothelial",
  "Endothelial:HEV1" = "Endothelial",
  "Endothelial:HEV2" = "Endothelial",
  "Endothelial:Endo1" = "Endothelial",
  "Endothelial:Endo2" = "Endothelial",
  "Endothelial:Endo3" = "Endothelial",
  "Endothelial:Endo4" = "Endothelial",
  "Endothelial:Endo5" = "Endothelial",
  "Endothelial:Endo6" = "Endothelial",
  "Endothelial:Endo7" = "Endothelial",
  "Endothelial:Stress" = "Endothelial",
  "Endothelial:Cell Cycle" = "Endothelial",
  "Endothelial:Interferon" = "Endothelial",
  "Endothelial:NF-kb" = "Endothelial",
  "Endothelial:Coagulation" = "Endothelial",

  # Macrophages
  "Macrophages:Lipid-associated" = "Macrophages/Monocytes",
  "Macrophages:Monocyte/Secreted" = "Macrophages/Monocytes",
  "Macrophages:Cell-cycle" = "Macrophages/Monocytes",
  "Macrophages:Interferon" = "Macrophages/Monocytes",
  "Macrophages:MES/Glycolysis" = "Macrophages/Monocytes",
  "Macrophages:MAC1" = "Macrophages/Monocytes",
  "Macrophages:MAC2" = "Macrophages/Monocytes",
  "Macrophages:MAC3" = "Macrophages/Monocytes",
  "Macrophages:Stress/HSP" = "Macrophages/Monocytes",
  "Macrophages:Proteasomal-degradation" = "Macrophages/Monocytes",
  "Macrophages:MYC/Mitochondria" = "Macrophages/Monocytes",
  "Macrophages:Unfolded-protein-response" = "Macrophages/Monocytes",
  "Macrophages:Respiration" = "Macrophages/Monocytes",

  # Fibroblasts
  "Fibroblasts :CAF1" = "Fibroblasts",
  "Fibroblasts :CAF2" = "Fibroblasts",
  "Fibroblasts :CAF3" = "Fibroblasts",
  "Fibroblasts :CAF4" = "Fibroblasts",
  "Fibroblasts :CAF5" = "Fibroblasts",
  "Fibroblasts :CAF6" = "Fibroblasts",
  "Fibroblasts :CAF7" = "Fibroblasts",
  "Fibroblasts :CAF8" = "Fibroblasts",
  "Fibroblasts :CAF9" = "Fibroblasts",
  "Fibroblasts :CAF10" = "Fibroblasts",
  "Fibroblasts :Cell-cycle" = "Fibroblasts",
  "Fibroblasts :Stress" = "Fibroblasts",
  "Fibroblasts :Hypoxia" = "Fibroblasts",
  "Fibroblasts :Complement" = "Fibroblasts",
  "Fibroblasts :Interferon" = "Fibroblasts",
  "Fibroblasts :Myofibroblasts" = "Fibroblasts",
  "Fibroblasts :PI16+" = "Fibroblasts",
  "Fibroblasts :Metal-response" = "Fibroblasts",
  "Fibroblasts :MHC-II/Cytokine" = "Fibroblasts",
  "Fibroblasts :Lipid-metabolism" = "Fibroblasts",
  "Fibroblasts :MHC-II" = "Fibroblasts",
  "Fibroblasts :Pericyte-like" = "Fibroblasts",
  "Unknown" = "Unknown"
)

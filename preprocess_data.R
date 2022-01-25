
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)

tmp_dir <- 'D:/MiamiProject/sci_scRNAseq/'
app_data_fraction <- 0.35
  
# Extract sci expression and dimreduc data -----------------------------------

sci <- readRDS(file = paste0(tmp_dir, 'data/sci.rds'))
DefaultAssay(sci) <- 'RNA'
sci[['SCT']] <- NULL
sci[['integrated']] <- NULL
sci <- NormalizeData(sci)
vars_sci <- sci[['RNA']]@meta.features
log_x_sci <- sci[['RNA']]@data
obs_sci <- sci@meta.data
umap_sci <- sci[['umap']]@cell.embeddings
colnames(umap_sci) <- c('sci_UMAP_1','sci_UMAP_2')

# # Features that have very low expression are not saved or loaded at all. Counts
# # sum cutoff equal to 
# format(object.size(log_x_sci[rownames(log_x_sci)[sparseMatrixStats::rowSums2(log_x_sci) > 0],]), 'Mb')

format(object.size(log_x_sci[rownames(log_x_sci)[sparseMatrixStats::rowSums2(log_x_sci) > 15],]), 'Mb')

# 
# # If total counts across all cells less than 5, dont even bother displaying
# almostzero_feats <- rownames(log_x_sci)[sparseMatrixStats::rowSums2(log_x_sci) <= 5]
# almostzero_feats_log_x_sci <- log_x_sci[almostzero_feats,]
# 
# # Variable Features should be loaded for every session
# sci <- FindVariableFeatures(sci, nfeatures = 6000)
# var_feats <- VariableFeatures(sci) # 230 genes are variable and almost zero
# table(var_feats %in% almostzero_feats)
# var_feats <- var_feats[!var_feats %in% almostzero_feats]
# var_feats_log_x_sci <- log_x_sci[var_feats,]
# 
# # All other genes. Need to split because many genes and total size is > 1Gb.
# rest_feats <- rownames(log_x_sci)[!rownames(log_x_sci) %in% union(var_feats, almostzero_feats)]
# rest_feats_log_x_sci <- log_x_sci[rest_feats,]
# split_genes <- split(1:nrow(rest_feats_log_x_sci), 1:3)
# rest_feats_1 <- rest_feats[split_genes[[1]]]
# rest_feats_2 <- rest_feats[split_genes[[2]]]
# rest_feats_3 <- rest_feats[split_genes[[3]]]
# rest_feats_log_x_sci_1 <- rest_feats_log_x_sci[rest_feats_1,]
# rest_feats_log_x_sci_2 <- rest_feats_log_x_sci[rest_feats_2,]
# rest_feats_log_x_sci_3 <- rest_feats_log_x_sci[rest_feats_3,]


# Extract myeloid metadata and dimreduc -----------------------------------

myeloid <- readRDS(file = paste0(tmp_dir, 'data/myeloid.rds'))
DefaultAssay(myeloid) <- 'RNA'
myeloid[['RNAcorrected']] <- NULL
myeloid[['integrated']] <- NULL
myeloid$Layer1_compartment <- 'Myeloid'
myeloid$Layer2_celltype <- myeloid$celltype
myeloid$Layer3_subtype <- myeloid$myeloid_subcluster
myeloid$preprint_subtype <- myeloid$old_subcluster
myeloid$myeloid_seurat_res.0.35 <- myeloid$integrated_snn_res.0.35
obs_myeloid <- c('Layer1_compartment','Layer2_celltype','Layer3_subtype','celltype','myeloid_subcluster','preprint_subtype','myeloid_seurat_res.0.35')
umap_myeloid <- myeloid[['umap']]@cell.embeddings
colnames(umap_myeloid) <- c('myeloid_UMAP_1','myeloid_UMAP_2')
obs_myeloid <- cbind(myeloid@meta.data[, obs_myeloid], umap_myeloid)
for (i in 1:ncol(obs_myeloid)) {
  if (class(obs_myeloid[[i]]) == 'factor') {
    obs_myeloid[[i]] <- as.character(obs_myeloid[[i]])
  }
}

# Extract vascular metadata and dimreduc -----------------------------------

vascular <- readRDS(file = paste0(tmp_dir, 'data/vascular.rds'))
DefaultAssay(vascular) <- 'RNA'
vascular[['RNAcorrected']] <- NULL
vascular[['integrated']] <- NULL
vascular$Layer1_compartment <- 'Vascular'
vascular$Layer2_celltype <- vascular$celltype
vascular$Layer3_subtype <- vascular$vascular_subcluster
vascular$preprint_subtype <- plyr::mapvalues(
  x = vascular$integrated_snn_res.0.4,
  from = 0:8,
  to = c('C1-Endothelial', 'C2-Endothelial','Tip Cell','A-Endothelial','U-Vascular','V-Endothelial','Fibroblast','Pericyte','VSMC')
)
vascular$vascular_seurat_res.0.4 <- vascular$integrated_snn_res.0.4
obs_vascular <- c('Layer1_compartment','Layer2_celltype','Layer3_subtype','celltype','vascular_subcluster','preprint_subtype','vascular_seurat_res.0.4')
umap_vascular <- vascular[['umap']]@cell.embeddings
colnames(umap_vascular) <- c('vascular_UMAP_1','vascular_UMAP_2')
obs_vascular <- cbind(vascular@meta.data[, obs_vascular], umap_vascular)
for (i in 1:ncol(obs_vascular)) {
  if (class(obs_vascular[[i]]) == 'factor') {
    obs_vascular[[i]] <- as.character(obs_vascular[[i]])
  }
}

# Extract macroglia metadata and dimreduc -----------------------------------

macroglia <- readRDS(file = paste0(tmp_dir, 'data/macroglia.rds'))
DefaultAssay(macroglia) <- 'RNA'
macroglia[['RNAcorrected']] <- NULL
macroglia[['integrated']] <- NULL
macroglia$Layer1_compartment <- 'Macroglia'
macroglia$Layer2_celltype <- macroglia$celltype
macroglia$Layer3_subtype <- macroglia$macroglia_subcluster
macroglia$preprint_subtype <- macroglia$macroglia_subcluster
macroglia$macroglia_seurat_res.0.4 <- macroglia$integrated_snn_res.0.4
obs_macroglia <- c('Layer1_compartment','Layer2_celltype','Layer3_subtype','celltype','macroglia_subcluster','preprint_subtype','macroglia_seurat_res.0.4')
umap_macroglia <- macroglia[['umap']]@cell.embeddings
colnames(umap_macroglia) <- c('macroglia_UMAP_1','macroglia_UMAP_2')
obs_macroglia <- cbind(macroglia@meta.data[, obs_macroglia], umap_macroglia)
for (i in 1:ncol(obs_macroglia)) {
  if (class(obs_macroglia[[i]]) == 'factor') {
    obs_macroglia[[i]] <- as.character(obs_macroglia[[i]])
  }
}

# Combine metadata ---------------------------------------------------------

# convenience function to pull myeloid/vascular/macroglia metadata in order of 
# barcodes in sci@meta.data rows
match_metadata <- function(order.char, obs.ls, col.char) {
  meta_dat <- rep(x = NA, times = length(order.char))
  for (i in 1:length(obs.ls)) {
    if (col.char %in% colnames(obs.ls[[i]])) {
      meta_match <- match(x = rownames(obs.ls[[i]]),
                          table = order.char)
      meta_dat[meta_match] <- obs.ls[[i]][, col.char]
    }
  }
  return(meta_dat)
}
cols_transfer <- c('Layer1_compartment', 'Layer2_celltype', 'Layer3_subtype',
                   'preprint_subtype', 'myeloid_UMAP_1', 'myeloid_UMAP_2',
                   'vascular_UMAP_1', 'vascular_UMAP_2', 'macroglia_UMAP_1',
                   'macroglia_UMAP_2', 'myeloid_subcluster', 
                   'vascular_subcluster', 'macroglia_subcluster')
obs_transfer <- data.frame(
  x = lapply(
    X = cols_transfer,
    FUN = match_metadata,
    order.char = rownames(obs_sci),
    obs.ls = list(obs_myeloid, obs_vascular, obs_macroglia)
  )
)
colnames(obs_transfer) <- cols_transfer

sci_meta_retain <- c('sample_id','time','orig.ident','nCount_RNA','nFeature_RNA','S.Score','G2M.Score','Phase','CC.Difference','percent_mt','percent_rp','doublet_scores','dissociationMethod','chemistry', 'celltype')
obs_sci <- do.call(
  what = cbind, 
  args = list(obs_transfer, obs_sci[, sci_meta_retain], umap_sci)
)


# Tidying metadata --------------------------------------------------------

# Rename neurons and lymphocyte compartments
obs_sci$Layer1_compartment[obs_sci$celltype == 'Neuron'] <- 'Neural'
obs_sci$Layer1_compartment[obs_sci$celltype == 'Lymphocyte'] <- 'Lymphocyte'
obs_sci$Layer2_celltype[obs_sci$celltype == 'Neuron'] <- 'Neuron'
obs_sci$Layer2_celltype[obs_sci$celltype == 'Lymphocyte'] <- 'Lymphocyte'
obs_sci$Layer3_subtype[obs_sci$celltype == 'Neuron'] <- 'Neuron'
obs_sci$Layer3_subtype[obs_sci$celltype == 'Lymphocyte'] <- 'Lymphocyte'
obs_sci$preprint_subtype[obs_sci$celltype == 'Neuron'] <- 'Neuron'
obs_sci$preprint_subtype[obs_sci$celltype == 'Lymphocyte'] <- 'Lymphocyte'

# Rename some metadata variables and set types
rename_meta <- c(
  Layer1_compartment = 'Compartment',
  Layer2_celltype = 'Celltype',
  Layer3_subtype = 'Subtype',
  time = 'InjuryTimePoint',
  myeloid_subcluster = 'Myeloid_Subcluster',
  vascular_subcluster = 'Vascular_Subcluster',
  macroglia_subcluster = 'Macroglia_Subcluster',
  sample_id = 'Sample_ID',
  percent_mt = 'Mito_Percent',
  percent_rp = 'Ribo_Percent',
  doublet_scores = 'DoubletScore',
  dissociationMethod = 'DissociationMethod',
  chemistry = 'Chemistry_10X',
  celltype = 'Celltype',
  Phase = 'CellCyclePhase',
  preprint_subtype = 'Preprint_Subtype'
)
colnames(obs_sci) <- plyr::mapvalues(
  x = colnames(obs_sci),
  from = names(rename_meta),
  to = rename_meta
)
meta_order <- c(
  # These are essential
  # 'Layer1_Compartment', 'Layer2_Celltype', 'Layer3_Subtype',
  'Compartment', 'Celltype', 'Subtype',
  'InjuryTimePoint', 'Sample_ID', 'CellCyclePhase', 
  # 'Celltype',
  # 'Myeloid_Subcluster',
  # 'Vascular_Subcluster', 
  # 'Macroglia_Subcluster',
  'Preprint_Subtype', 
  'DissociationMethod', 
  'Chemistry_10X', 
  # 'nCount_RNA', 
  # 'nFeature_RNA', 
  # 'Mito_Percent',
  # 'Ribo_Percent', 
  # 'S.Score', 
  # 'G2M.Score', 
  # 'DoubletScore', 
  'sci_UMAP_1', 
  'sci_UMAP_2', 
  'myeloid_UMAP_1', 
  'myeloid_UMAP_2', 
  'vascular_UMAP_1', 
  'vascular_UMAP_2',
  'macroglia_UMAP_1', 
  'macroglia_UMAP_2'
)
obs_sci <- obs_sci[meta_order]

obs_sci$Compartment <- factor(
  x = obs_sci$Compartment,
  levels = c('Myeloid', 'Vascular', 'Macroglia', 'Neural', 'Lymphocyte')
)
obs_sci$Celltype <- factor(
  x = obs_sci$Celltype,
  levels = c('Neutrophil','Monocyte','Macrophage','Dendritic','Microglia',
             'Div-Myeloid','Fibroblast','Endothelial','Pericyte','OPC',
             'Oligodendrocyte', 'Astrocyte', 'Ependymal', 'Neuron', 
             'Lymphocyte')
)
obs_sci$Subtype <- factor(
  x = obs_sci$Subtype,
  levels = c("Neutrophil", "Monocyte", "Chemotaxis-Inducing Mac", "Inflammatory Mac", "Border-Associated Mac", "Dendritic", "Dividing Myeloid", "Homeostatic Microglia", "Inflammatory Microglia", "Dividing Microglia", "Migrating Microglia", "Interferon Myeloid", "A-Endothelial", "C-Endothelial", "V-Endothelial", "Tip Cell", "Pericyte", "VSMC", "Fibroblast", "U-Vascular", "Ependymal-A", "Ependymal-B", "Astroependymal", "Astrocyte", "OPC-A", "OPC-B", "Div-OPC", "Pre-Oligo", "Oligodendrocyte")
)
obs_sci$Preprint_Subtype <- factor(
  x = obs_sci$Preprint_Subtype,
  levels = c("Neutrophil", "Monocyte", "Macrophage-A", "Macrophage-B", "BA-Macrophage", "Dendritic", "Div-Myeloid", "H-Microglia", "DAM-A", "DAM-B", "DAM-C", "IFN-Myeloid", "A-Endothelial", "C1-Endothelial", "C2-Endothelial", "V-Endothelial", "Tip Cell", "Pericyte", "VSMC", "Fibroblast", "U-Vascular", "Ependymal-A", "Ependymal-B", "Astroependymal", "Astrocyte", "OPC-A", "OPC-B", "Div-OPC", "Pre-Oligo", "Oligodendrocyte")
)

rm(myeloid, vascular, macroglia, obs_myeloid, obs_vascular, obs_macroglia,
   umap_myeloid, umap_vascular, umap_macroglia, obs_transfer, sci_meta_retain,
   meta_order, rename_meta, cols_transfer, match_metadata, i)

# Gene-level metadata -----------------------------------------------------

# Including which variables genes were used for clustering in each compartment
# might be helpful information but it's pretty low priority. Will push back for
# a later date.



# Other related data ---------------------------------------------------

# dataset key-value pairs for subsetting data
dataset_dict <- c(
  All_SCI = 'sci', 
  Myeloid = 'myeloid',
  Vascular = 'vascular',
  Macroglia = 'macroglia'
)

# Categorical variables by which cells can be grouped.
categorical_vars <- sapply(obs_sci, function(x) !is.numeric(x))
cell_groupings <- colnames(obs_sci)[categorical_vars]
for (i in 1:ncol(obs_sci)) {
  if (class(obs_sci[[i]]) == 'character') {
    obs_sci[[i]] <- factor(obs_sci[[i]])
  }
}

# All possible genes to query
all_features <- rownames(vars_sci)

# Cell Counts
cell_counts <- obs_sci %>% count(InjuryTimePoint)

# dimplot label coordinates
umap_vars <- grep(pattern = 'UMAP', colnames(obs_sci), value = TRUE)
label_coords <- obs_sci[c(cell_groupings, umap_vars)] %>%
  pivot_longer(cols = !contains(match = 'UMAP'), 
               names_to = 'variables', 
               values_to = 'values') %>%
  group_by(variables, values) %>% 
  summarise(sci_UMAP_1 = median(sci_UMAP_1, na.rm = TRUE),
            sci_UMAP_2 = median(sci_UMAP_2, na.rm = TRUE),
            myeloid_UMAP_1 = median(myeloid_UMAP_1, na.rm = TRUE),
            myeloid_UMAP_2 = median(myeloid_UMAP_2, na.rm = TRUE),
            vascular_UMAP_1 = median(vascular_UMAP_1, na.rm = TRUE),
            vascular_UMAP_2 = median(vascular_UMAP_2, na.rm = TRUE),
            macroglia_UMAP_1 = median(macroglia_UMAP_1, na.rm = TRUE),
            macroglia_UMAP_2 = median(macroglia_UMAP_2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'variables', values_from = 'values') %>% 
  as.data.frame()
# # Testing label coordinates
# ggplot(data = label_coords) +
#   geom_text(mapping = aes(x = myeloid_UMAP_1,
#                           y = myeloid_UMAP_2,
#                           label = Sample_ID))

rm(categorical_vars, umap_vars)


# Website details ---------------------------------------------------------

window_title = "Single-Cell RNAseq of Mouse Spinal Cord Injury"
title_link_text = "by Jae Lee lab"
title_link_url = "https://www.jaeleelab.com"

# Save data ---------------------------------------------------------------

# Since free version of ShinyApps server only allows 1Gb of memory, need to 
# reduce the number of cells for which data is displayed. This is set by the 
# `app_data_fraction` variable set at the top of this script.

keep_cells <- sample(
  x = colnames(log_x_sci), 
  size = trunc(ncol(log_x_sci) * app_data_fraction),
  replace = FALSE
)
log_x_sci <- log_x_sci[, keep_cells]
obs_sci <- obs_sci[keep_cells,]

save(log_x_sci,
     obs_sci,
     vars_sci,
     label_coords,
     dataset_dict,
     cell_groupings,
     all_features,
     cell_counts,
     window_title,
     title_link_text,
     title_link_url,
     file = './data/appdata.RData')

# # Save a key for when server requests genes
# feature_split_dict <- list(
#   var_feats = var_feats,
#   almostzero_feats = almostzero_feats,
#   rest_feats_1 = rest_feats_1,
#   rest_feats_2 = rest_feats_2,
#   rest_feats_3 = rest_feats_3
# )
# 
# saveRDS(rest_feats_log_x_sci_1, file = './data/rest_feats_log_x_sci_1.rds')
# saveRDS(rest_feats_log_x_sci_2, file = './data/rest_feats_log_x_sci_2.rds')
# saveRDS(rest_feats_log_x_sci_3, file = './data/rest_feats_log_x_sci_3.rds')
# saveRDS(var_feats_log_x_sci, file = './data/var_feats_log_x_sci.rds')
# saveRDS(almostzero_feats_log_x_sci, file = './data/almostzero_feats_log_x_sci.rds')
# saveRDS(feature_split_dict, file = './data/feature_split_dict.rds')
# saveRDS(obs_sci, file = './data/obs_sci.rds')
# saveRDS(vars_sci, file = './data/vars_sci.rds')
# saveRDS(label_coords, file = './data/label_coords.rds')
# saveRDS(log_x_sci, file = './data/log_x_sci.rds') # will not really need this

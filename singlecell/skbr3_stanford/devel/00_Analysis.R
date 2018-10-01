install.packages('Seurat')
packageVersion("Seurat")
library(Seurat)
library(dplyr)
library(tibble)

setwd("/Users/charlie/setup/datahub/public/skbr3_stanford/devel")

#for single cell data we map each assay to a patient
#                     we map each cell to a sample

#required attributes for assay
assay_data = list(
  list(PATIENT_ID="skbr3_control", DATA=Read10X(data.dir = "skbr3_control/"), ASSAY_CONDITION="CONTROL"),
  list(PATIENT_ID="skbr3_treated", DATA=Read10X(data.dir = "skbr3_treated/"), ASSAY_CONDITION="TEMPEST")
)
for(i in seq_along(assay_data)) { assay_data[[i]][["ASSAY_ID"]] = assay_data[[i]]$PATIENT_ID }
for(i in seq_along(assay_data)) { assay_data[[i]][["N_CELLS"]] = ncol(assay_data[[i]]$DATA) }
for(i in seq_along(assay_data)) { assay_data[[i]][["ID_CELLS"]] = paste(assay_data[[i]]$PATIENT_ID, colnames(assay_data[[i]]$DATA), sep="-") }
for(i in seq_along(assay_data)) { colnames(assay_data[[i]][["DATA"]]) = assay_data[[i]][["ID_CELLS"]] }
id_assays = sapply(assay_data, "[[", "PATIENT_ID")
n_cells =  sapply(assay_data, "[[", "N_CELLS") #cells per assay
names(assay_data) = id_assays

#optional attributes for assay
#data_clinical_patient.txt
Patient_Header_List = list(
  c("PATIENT_ID", "Assay identifier alias as Patient ID", "sBioPortal", "STRING", 1),
  c("ASSAY_ID", "Assay identifier", "sBioPortal", "STRING", 1),
  c("ASSAY_CONDITION", "Assay Condition", "sBioPortal", "STRING", 1)
)
data_clinical_patient = new_tibble(list(), nrow=length(assay_data)) %>%
  add_column(PATIENT_ID = sapply(assay_data, "[[", "PATIENT_ID")) %>%
  add_column(ASSAY_ID = sapply(assay_data, "[[", "ASSAY_ID")) %>%
  add_column(ASSAY_CONDITION = sapply(assay_data, "[[", "ASSAY_CONDITION"))
data_clinical_patient_file = "~/setup/datahub/public/skbr3_stanford/devel/data_clinical_patient.txt"
con <- file(data_clinical_patient_file, open="wt")
comment = sapply(2:5, function(x, l, c='#') {
  cells = sapply(l,"[[", x)
  return(paste(c, paste(cells, collapse='\t'), sep=''))
}, Patient_Header_List)
writeLines(comment, con) 
write.table(data_clinical_patient, row.names=F, sep='\t', quote=F, con) 
close(con)

#data_clinial_sample.txt
Sample_Header_List = list(
  c("PATIENT_ID", "Assay identifier", "sBioPortal", "STRING", 1),                
  c("SAMPLE_ID", "Cell identifier", "sBioPortal", "STRING", 1)
)
data_clinical_sample = new_tibble(list(), nrow=sum(n_samples)) %>%
  add_column(PATIENT_ID = rep(id_assays, n_cells)) %>%
  add_column(SAMPLE_ID = unlist(sapply(assay_data, "[[", "ID_CELLS")) )
data_clinical_sample_file = "~/setup/datahub/public/skbr3_stanford/devel/data_clinical_sample.txt"
con <- file(data_clinical_sample_file, open="wt")
comment = sapply(2:5, function(x, l, c='#') {
  cells = sapply(l,"[[", x)
  return(paste(c, paste(cells, collapse='\t'), sep=''))
}, Sample_Header_List)
writeLines(comment, con) 
write.table(data_clinical_sample, row.names=F, sep='\t', quote=F, con) 
close(con)

#data_expression_file.txt
data_expression_file = "~/setup/datahub/public/skbr3_stanford/devel/data_expression_file.txt"
con = file(data_expression_file, open="wt")
comment = paste(c("Hugo_Symbol", do.call("c", sapply(assay_data, "[[", "ID_CELLS"))), collapse="\t")
writeLines(comment, con) 
write.table(as(do.call("cbind", sapply(assay_data, "[[", "DATA")), "matrix"), #convert to ordinary matrix
                  row.names=T, col.names = F, sep='\t', quote=F, con)
close(con)
#76MB -> 500MB, 7 times waste

# Set up control object
# colnames(x = pbmc_ctrl_8.data) <- paste('ctrl_8', colnames(x = pbmc_ctrl_8.data), sep = '_')
# ctrl <- CreateSeuratObject(raw.data = pbmc_ctrl_8.data, project = "media", min.cells = 3)
# ctrl@meta.data$stim <- "CTRL"
# ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 200, high.thresholds = Inf)
# ctrl <- NormalizeData(ctrl)
# ctrl <- ScaleData(ctrl, display.progress = F)
# # Set up stimulated object
# colnames(x = pbmc_tpst_8.data) <- paste('tpst_8', colnames(x = pbmc_tpst_8.data), sep = '_')
# stim <- CreateSeuratObject(raw.data = pbmc_tpst_8.data, project = "tpst", min.cells = 3)
# stim@meta.data$stim <- "STIM"
# stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 200, high.thresholds = Inf)
# stim <- NormalizeData(stim)
# stim <- ScaleData(stim, display.progress = F)
# 
# ctrl <- FindVariableGenes(ctrl, do.plot = F)
# stim <- FindVariableGenes(stim, do.plot = F)
# 
# length(x = ctrl@var.genes)
# length(x = stim@var.genes)
# 
# g.1 <- head(rownames(ctrl@hvg.info), 1000)
# g.2 <- head(rownames(stim@hvg.info), 1000)
# genes.use <- unique(c(g.1, g.2))
# genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
# genes.use <- intersect(genes.use, rownames(stim@scale.data))
# 
# immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)
# # visualize results of CCA plot CC1 versus CC2 and look at a violin plot
# p1 <- DimPlot(object = immune.combined, reduction.use = "cca", group.by = "stim", 
#               pt.size = 0.5, do.return = TRUE)
# p2 <- VlnPlot(object = immune.combined, features.plot = "CC1", group.by = "stim", 
#               do.return = TRUE)
# plot_grid(p1, p2)
# 
# PrintDim(object = immune.combined, reduction.type = "cca", dims.print = 1:2, 
#          genes.print = 10)
# 
# immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", 
#                                  dims.align = 1:20)
# 
# p1 <- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "stim", 
#               do.return = TRUE)
# p2 <- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "stim", 
#               do.return = TRUE)
# plot_grid(p1, p2)
# 
# 
# immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, 
#                            do.fast = T)
# immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", 
#                                 resolution = 0.6, dims.use = 1:20)
# 
# p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
# p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
# plot_grid(p1, p2)
# 
# FeaturePlot(object = immune.combined, features.plot = c("CD3D", "SELL", "CREM", 
#                                                         "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), min.cutoff = "q9", cols.use = c("lightgrey", 
#                                                                                                                                             "blue"), pt.size = 0.5)
# 
# immune.combined <- BuildClusterTree(object = immune.combined, do.reorder = TRUE, reorder.numeric = TRUE)
# 
# TSNEPlot(object = immune.combined, do.label = T)
# 
# immune.combined.marker<- FindAllMarkers(object = immune.combined, only.pos = TRUE, min.pct = 0.10, thresh.use = 0.4)
# 
# immune.types <- c("CD3D","CD79A", "CD14", "CD8A", "GNLY", "IL1B", "CREM","IL17F", "CD69")
# 
# FeaturePlot(object = immune.combined, features.plot = immune.types, no.legend = FALSE, min.cutoff = "q10", max.cutoff = "q90", cols.use = c("lightgrey","red"), pt.size = 0.5)
# 
# table(immune.combined@ident, immune.combined@meta.data$orig.ident)
# 
# top5.cca <-immune.combined.marker %>% group_by(cluster) %>% top_n(5, avg_diff)
# DoHeatmap(object = SubsetData(object = immune.combined, max.cells.per.ident = 100), genes.use = top5.cca$gene, slim.col.label = TRUE, group.label.rot = TRUE)
# 
# VlnPlot(object = immune.combined, features.plot = c("nGene", "nUMI"), nCol = 2, point.size.use = 0.1)
# 
# ####test
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = immune.combined@data), value = TRUE)
# 
# percent.mito <- Matrix::colSums(immune.combined@raw.data[mito.genes, ])/Matrix::colSums(immune.combined@raw.data)
# 
# immune.combined <- AddMetaData(object = immune.combined, metadata = percent.mito, col.name = "percent.mito")
# 
# VlnPlot(object = immune.combined, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, point.size.use = 0.1)
# 
# FeatureHeatmap(immune.combined, features.plot = c("FLT3LG", "CSF2", "CXCL1", "IFNG", "CCL7", "IL13", "IL1RN", "LTA"), group.by = "stim", pt.size = 0.25, key.position = "top", max.exp = 3)
#                max.exp = 3)
# immune.combined@meta.data$celltype.stim <- paste0(immune.combined@ident, "_", immune.combined@meta.data$stim)
# head(immune.combined@meta.data)
# 
# immune.combined <- StashIdent(immune.combined, save.name = "celltype")
# 
# immune.combined <- SetAllIdent(immune.combined, id = "celltype.stim")
# immune.combined <- SetAllIdent(immune.combined, id = "tree.ident")
# tpst_8.response <- FindMarkers(immune.combined, ident.1 = "1_STIM", ident.2 = "1_CTRL", print.bar = FALSE)
# head(tpst_8.response, 20)
# View(tpst_8.response)
# write.table(tpst_8.response, file ="tpst_8_mono_response.txt", col.names = NA, quote = FALSE, sep="\t")
# saveRDS(immune.combined, file ="/mnt/ix1/Projects/M048_180606_Tempest/pbmc_IA_aggr_8/6314_pbmc_IA_aggr_8/outs/filtered_gene_bc_matrices_mex/GRCh38/immune.combined_8.rds")
# 
# TSNEPlot(object = immune.combined, do.label = T)
# immune.combined <- StashIdent(immune.combined, save.name = "celltype")
# immune.combined <- SetAllIdent(immune.combined, id = "orig.ident")
# VlnPlot(object = immune.combined, features.plot = c("FLT3LG", "CSF2", "CXCL1", "IFNG", "CCL7", "IL13", "IL1RN", "LTA"), point.size.use = 0.1)
# 

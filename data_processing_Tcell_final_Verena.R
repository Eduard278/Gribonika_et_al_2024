library(Seurat)
library(stringr)
library(Matrix)
library(gplots)
library(ggplot2)
library(cellhashR)

Tcells_aggr_readin <- Read10X("gene_expression_data/filtered_feature_bc_matrix/")


#Sumarize HTO hashtags by tissue
m <- Matrix(nrow = 3, ncol = 27818, data = 0, sparse = TRUE)
m <- as(m, "dgTMatrix") # by default, Matrix() returns dgCMatrix
m@Dimnames[[1]] <- c("skin", "LN", "spleen")
m@Dimnames[[2]] <- Tcells_aggr_readin$Custom@Dimnames[[2]]

m[1,] <- Tcells_aggr_readin$Custom['Mouse1_F_skin',] + Tcells_aggr_readin$Custom['Mouse2_F_skin',] +
  Tcells_aggr_readin$Custom['Mouse3_F_skin',] + Tcells_aggr_readin$Custom['Mouse4_F_skin',] + Tcells_aggr_readin$Custom['Mouse5_M_skin',] +
  Tcells_aggr_readin$Custom['Mouse6_M_skin',] + Tcells_aggr_readin$Custom['Mouse7_M_skin',] + Tcells_aggr_readin$Custom['Mouse8_M_skin',] +
  Tcells_aggr_readin$Custom['Mouse9_M_skin',] + Tcells_aggr_readin$Custom['Mouse10_M_skin',]
m[2,] <- Tcells_aggr_readin$Custom['Mouse1_F_LN',] + Tcells_aggr_readin$Custom['Mouse2_F_LN',] + Tcells_aggr_readin$Custom['Mouse3_F_LN',] +
  Tcells_aggr_readin$Custom['Mouse4_F_LN',] + Tcells_aggr_readin$Custom['Mouse5_M_LN',] + Tcells_aggr_readin$Custom['Mouse6_M_LN',]
m[3,] <- Tcells_aggr_readin$Custom['Mouse1_F_spleen',] + Tcells_aggr_readin$Custom['Mouse2_F_spleen',] + Tcells_aggr_readin$Custom['Mouse3_F_spleen',] +
  Tcells_aggr_readin$Custom['Mouse4_F_spleen',] + Tcells_aggr_readin$Custom['Mouse5_M_spleen',] + Tcells_aggr_readin$Custom['Mouse6_M_spleen',]


Tcells_aggr <- CreateSeuratObject(counts = Tcells_aggr_readin[["Gene Expression"]], min.cells = 3, min.features = 200)
Tcells_aggr[["HTO"]] <- CreateAssayObject(Tcells_aggr_readin[["Custom"]][, colnames(x = Tcells_aggr)])
Tcells_aggr[["HTO_combined"]] <- CreateAssayObject(m[, colnames(x = Tcells_aggr)])

Idents(Tcells_aggr) <- "orig.ident"
Tcells_aggr[["percent.mt"]] <- PercentageFeatureSet(Tcells_aggr, assay="RNA",pattern = "mt-")
VlnPlot(Tcells_aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
Tcells_aggr_filtered <- subset(Tcells_aggr, subset = nFeature_RNA  > 200 & nFeature_RNA < 3000 & percent.mt < 5 & nCount_RNA < 10000)
VlnPlot(Tcells_aggr_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)


Tcells_aggr_filtered <- NormalizeData(Tcells_aggr_filtered, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000)
Tcells_aggr_filtered <- NormalizeData(Tcells_aggr_filtered, assay = "HTO",normalization.method = "CLR")
Tcells_aggr_filtered <- NormalizeData(Tcells_aggr_filtered, assay = "HTO_combined",normalization.method = "CLR")


Tcells_aggr_filtered <- HTODemux(Tcells_aggr_filtered, assay = "HTO_combined", positive.quantile = 0.99)


RidgePlot(Tcells_aggr_filtered, assay = "HTO_combined", features = rownames(Tcells_aggr_filtered[["HTO_combined"]]), ncol = 2)

barcodeData_Tcells <- ProcessCountMatrix(rawCountData = 'HTO_data/filtered_feature_bc_matrix/', 
                                         minCountPerCell = 5, datatypeName = 'Custom')


pdf("plots_from_cellhash_calls_Tcells_combined_HTO.pdf")
calls_Tcells_combined <- GenerateCellHashingCalls(barcodeMatrix = barcodeData_Tcells, methods = c('multiseq', 'htodemux', 'bff_cluster', 'bff_raw', 'dropletutils'))
dev.off()


#For discordant calls use the multiseq results
discordant_Tcells <- subset(calls_Tcells_combined, calls_Tcells_combined[['consensuscall']] == "Discordant")
non_discordant_Tcells <- subset(calls_Tcells_combined, calls_Tcells_combined[['consensuscall']] != "Discordant")

discordant_Tcells$consensuscall <- discordant_Tcells$multiseq
discordant_Tcells$consensuscall.global <- rep('Singlet', length(discordant_Tcells$consensuscall))

calls_combined_Tcells <- rbind(discordant_Tcells, non_discordant_Tcells)


#Add classification to Seurat object
rownames(calls_combined_Tcells) <- calls_combined_Tcells$cellbarcode
Tcells_aggr_filtered <- AddMetaData(object = Tcells_aggr_filtered,metadata = calls_combined_Tcells['consensuscall'],col.name = 'consensuscall_combined')
Tcells_aggr_filtered <- AddMetaData(object = Tcells_aggr_filtered,metadata = calls_combined_Tcells['consensuscall.global'],col.name = 'consensuscall_combined.global')

#filter data
Tcells_aggr_filtered_demux_tissue <- subset(Tcells_aggr_filtered, HTO_combined_classification.global == "Singlet")

DefaultAssay(Tcells_aggr_filtered_demux_tissue) <- "RNA"

Tcells_aggr_filtered_demux_tissue <- NormalizeData(Tcells_aggr_filtered_demux_tissue)
Tcells_aggr_filtered_demux_tissue <- FindVariableFeatures(Tcells_aggr_filtered_demux_tissue)
Tcells_aggr_filtered_demux_tissue <- ScaleData(Tcells_aggr_filtered_demux_tissue)
Tcells_aggr_filtered_demux_tissue <- RunPCA(Tcells_aggr_filtered_demux_tissue, verbose = FALSE)
Tcells_aggr_filtered_demux_tissue <- FindNeighbors(Tcells_aggr_filtered_demux_tissue, dims = 1:17)
Tcells_aggr_filtered_demux_tissue <- FindClusters(Tcells_aggr_filtered_demux_tissue, resolution = 0.8, verbose = FALSE)
Tcells_aggr_filtered_demux_tissue <- RunUMAP(Tcells_aggr_filtered_demux_tissue, dims = 1:14)

#Rename samples
Tcells_aggr_filtered_demux_tissue@meta.data$sample <- Tcells_aggr_filtered_demux_tissue@meta.data$HTO_combined_maxID
Tcells_aggr_filtered_demux_tissue@meta.data$origin <- Tcells_aggr_filtered_demux_tissue@meta.data$sample
Tcells_aggr_filtered_demux_tissue@meta.data$origin <- str_remove(Tcells_aggr_filtered_demux_tissue@meta.data$origin, "Mouse.+-.-")

#Remove outliers
Tcells_aggr_filtered_demux_tissue_filtered <- subset(Tcells_aggr_filtered_demux_tissue, seurat_clusters != "7" & seurat_clusters !="13")

Tcells_aggr_filtered_demux_tissue_filtered <- NormalizeData(Tcells_aggr_filtered_demux_tissue_filtered)
Tcells_aggr_filtered_demux_tissue_filtered <- FindVariableFeatures(Tcells_aggr_filtered_demux_tissue_filtered)
Tcells_aggr_filtered_demux_tissue_filtered <- ScaleData(Tcells_aggr_filtered_demux_tissue_filtered)
Tcells_aggr_filtered_demux_tissue_filtered <- RunPCA(Tcells_aggr_filtered_demux_tissue_filtered, verbose = FALSE)
Tcells_aggr_filtered_demux_tissue_filtered <- FindNeighbors(Tcells_aggr_filtered_demux_tissue_filtered, dims = 1:17)
Tcells_aggr_filtered_demux_tissue_filtered <- FindClusters(Tcells_aggr_filtered_demux_tissue_filtered, resolution = 0.8, verbose = FALSE)

pct <- Tcells_aggr_filtered_demux_tissue_filtered[["pca"]]@stdev / sum(Tcells_aggr_filtered_demux_tissue_filtered[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1

#Determine where the percent change in variation between consectuve PCs is less than 0.1%
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2

Tcells_aggr_filtered_demux_tissue_filtered <- RunUMAP(Tcells_aggr_filtered_demux_tissue_filtered, dims = 1:15)
DimPlot(Tcells_aggr_filtered_demux_tissue_filtered, label = TRUE)

#Remove cluster 18 - keratinocytes
Tcells_aggr_filtered_demux_tissue_filtered2 <- subset(Tcells_aggr_filtered_demux_tissue_filtered, seurat_clusters != "18")

Tcells_aggr_filtered_demux_tissue_filtered2 <- NormalizeData(Tcells_aggr_filtered_demux_tissue_filtered2)
Tcells_aggr_filtered_demux_tissue_filtered2 <- FindVariableFeatures(Tcells_aggr_filtered_demux_tissue_filtered2)
Tcells_aggr_filtered_demux_tissue_filtered2 <- ScaleData(Tcells_aggr_filtered_demux_tissue_filtered2)
Tcells_aggr_filtered_demux_tissue_filtered2 <- RunPCA(Tcells_aggr_filtered_demux_tissue_filtered2, verbose = FALSE)
Tcells_aggr_filtered_demux_tissue_filtered2 <- FindNeighbors(Tcells_aggr_filtered_demux_tissue_filtered2, dims = 1:17)
Tcells_aggr_filtered_demux_tissue_filtered2 <- FindClusters(Tcells_aggr_filtered_demux_tissue_filtered2, resolution = 0.8, verbose = FALSE)

pct <- Tcells_aggr_filtered_demux_tissue_filtered2[["pca"]]@stdev / sum(Tcells_aggr_filtered_demux_tissue_filtered2[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1

#Determine where the percent change in variation between consectuve PCs is less than 0.1%
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2

Tcells_aggr_filtered_demux_tissue_filtered2 <- RunUMAP(Tcells_aggr_filtered_demux_tissue_filtered2, dims = 1:14)

#Cluster assignments
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$seurat_clusters
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^0$", "Naive I")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^1$", "activated memory Th17-like")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^2$", "Follicular-like Tregs")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^3$", "Tregs")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^4$", "Naive II")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^5$", "tTregs")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^6$", "Th1-like TFH")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^7$", "Naive III")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^8$", "Th17-like TFH")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^9$", "memory-like TFH")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^10$", "TFR")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^11$", "resident memory T cells")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^12$", "Central memory")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^13$", "Th1 cells I")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^14$", "effector memory T cells")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^15$", "activated T cells")
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes <- str_replace(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes, "^16$", "Th1 cells II")


Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes_summary <- Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes
Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes_summary <- str_remove(Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$celltypes_summary, " I.*")


#Add VDJ data to this object
#Add TCR data
TCR_data <- read.csv("filtered_contig_annotations_TCR.csv")
#Just keep cells that were not filtered out during scRNA_seq QC
TCR_data_filtered <- subset(TCR_data, TCR_data$barcode %in% names(Tcells_aggr_filtered_demux_tissue_filtered2@active.ident))

#Add samples and treatment to VDJ data
sample_info <- data.frame(barcode = colnames(Tcells_aggr_filtered_demux_tissue_filtered2), 
                          samples = Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$HTO_combined_classification, 
                          ID = Tcells_aggr_filtered_demux_tissue_filtered2@meta.data$consensuscall)

TCR_data_with_info <- merge(TCR_data_filtered, sample_info, by = "barcode")


combined <- combineTCR(TCR_data_filtered, 
                       samples  = c("TCR"))
combined$TCR$barcode <- str_remove(combined$TCR$barcode, "TCR_")

Tcells_aggr_filtered_demux_tissue_filtered2 <- combineExpression(combined, Tcells_aggr_filtered_demux_tissue_filtered2, 
                                                                 cloneCall="gene", 
                                                                 group.by = "sample", 
                                                                 proportion = FALSE, 
                                                                 cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

#Save data and read in
# saveRDS(Tcells_aggr_filtered_demux_tissue_filtered2, file = "Seurat_Tcells_final_object_with_VDJ.rds")
# Tcells_aggr_filtered_demux_tissue_filtered2 <- readRDS(file = "Seurat_Tcells_final_object_with_VDJ.rds")

#Figure 2F
Idents(Tcells_aggr_filtered_demux_tissue_filtered2) <- "celltypes_summary"
DimPlot(Tcells_aggr_filtered_demux_tissue_filtered2, pt.size = 4,
        cols = c("#FFFF00", "#bfef45", "#000000", "#3cb44b", "#469990", "#8A3324", "#FF00FF", "#42d4f4", 
                 "#000075", "#FF0e0e", "#C154C1", "#89cff0", "#AB3F00", "#F88379"))

#Figure S2E
gene_list <- c("Cd40lg", "Icos", "Pdcd1", "Bcl6", "Foxp3", "Il2ra", "Tigit", "Gzmb", "Il10", "Nrp1", "Ikzf2", "Tbx21",
               "Ifng", "Id3", "Tnf", "Cxcr4", "Il12rb2", "Gata3", "Il4", "Rorc", "Rora", "Il23r", "Ccr6", "Cd44", 
               "Cd69", "Ccr7", "Sell")
DotPlot(Tcells_aggr_filtered_demux_tissue_filtered2, features = gene_list) + theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "firebrick")

#Figure S2F
Idents(Tcells_aggr_filtered_demux_tissue_filtered2) <- "HTO_combined_maxID"
DimPlot(Tcells_aggr_filtered_demux_tissue_filtered2, pt.size = 4, 
        cols = c("#c68642", "#00008B", "#FF0000"))

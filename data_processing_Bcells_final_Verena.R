library(Seurat)

Tcells_aggr_readin <- Read10X("gene_expression_data/filtered_feature_bc_matrix/")

#Extract LTa-/-  B cells
scRNA_LTa_Bcells <- NormalizeData(scRNA_LTa_Bcells)
scRNA_LTa_Bcells.genes <- rownames(scRNA_LTa_Bcells)
scRNA_LTa_Bcells <- FindVariableFeatures(scRNA_LTa_Bcells, selection.method ="vst", nfeatures = 2000)
scRNA_LTa_Bcells <- ScaleData(scRNA_LTa_Bcells, features = scRNA_LTa_Bcells.genes)
scRNA_LTa_Bcells <- RunPCA(scRNA_LTa_Bcells, features = VariableFeatures(object = scRNA_LTa_Bcells))
scRNA_LTa_Bcells <- FindNeighbors(scRNA_LTa_Bcells, dims = 1:10)
scRNA_LTa_Bcells <- FindClusters(scRNA_LTa_Bcells, resolution = 1)
scRNA_LTa_Bcells <- RunUMAP(scRNA_LTa_Bcells, dims = 1:10)
Idents(scRNA_LTa_Bcells) <- "seurat_clusters"
DimPlot(scRNA_LTa_Bcells, reduction = "umap", pt.size = 2, label = TRUE)

##################cell annotations annotations ################################

celltypes_manual <- scRNA_LTa_Bcells@meta.data$seurat_clusters
celltypes_manual <- str_replace(celltypes_manual, '^0', "P0 - Pre-GCB-like")
celltypes_manual <- str_replace(celltypes_manual, '^1', "P1 - Activated B-like")
celltypes_manual <- str_replace(celltypes_manual, '^2', "P2 - Mature B-like")
celltypes_manual <- str_replace(celltypes_manual, '^3', "P3 - Follicular B-like")
celltypes_manual <- str_replace(celltypes_manual, '^4', "P4 - Germinal Center B-like")
celltypes_manual <- str_replace(celltypes_manual, '^5', "P5 - Cycling B-like")
celltypes_manual <- str_replace(celltypes_manual, '^6', "P6 - Marginal zone B-like")
celltypes_manual <- str_replace(celltypes_manual, '^7', "P7 - Immunoregulatory B-like")
celltypes_manual <- str_replace(celltypes_manual, '^8', "P8 - Memory B-like")
###############################################################################


############################ Figure 5F #######################################
scRNA_LTa_Bcells <- AddMetaData(scRNA_LTa_Bcells, celltypes_manual, col.name = "celltypes_manual")
Idents(scRNA_LTa_Bcells) <- "celltypes_manual"
DimPlot(scRNA_LTa_Bcells, reduction = "umap", pt.size = 4, 
        cols = c("#469990", "#000075", "#ffd8b1", "#e6194B", "#f58231", "#3cb44b", "#42d4f4", "#f032e6", "#000000"), split.by = "condition") 
###############################################################################



############################ Figure S5J #######################################
###### Number of cells per cluster for B cells
table_celltypes_cond <- as.data.frame(table(scRNA_LTa_Bcells@meta.data$celltypes_manual, scRNA_LTa_Bcells@meta.data$condition))
table_celltypes_cond_TA <- subset(table_celltypes_cond, table_celltypes_cond$Var2 == "TA")
table_celltypes_cond_TSB <- subset(table_celltypes_cond, table_celltypes_cond$Var2 == "TSB")

table_celltypes_cond_TA$Perc <- (table_celltypes_cond_TA$Freq/sum(table_celltypes_cond_TA$Freq)) * 100
table_celltypes_cond_TSB$Perc <- (table_celltypes_cond_TSB$Freq/sum(table_celltypes_cond_TSB$Freq)) * 100

table_celltypes_cond <- rbind(table_celltypes_cond_TA, table_celltypes_cond_TSB)

custom.color.fun_raw <- colorRampPalette(colors = c("white", "lightblue", "turquoise", "blue4"), space = "rgb")
ggplot(table_celltypes_cond, aes(x = Var1, y = Var2, size = Perc, col = Freq)) + geom_point() + 
  scale_colour_gradientn(colours = custom.color.fun_raw(255)) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(title = "Bcells")
###############################################################################


############################ Figure 5G #######################################
gene_list_for_publication <- c("Cd74", "Ezr", "Trim25", "Vps37b", "Pou2af1", "Limd2", "Nfkbia", "Tnfrsf13c", "Dusp2",
                               "Junb", "Egr1", "Myc", "Gpr183", "Cd83", "Bcl2a1b", "Cxcr5", "Vim", "Ucp2", "Plac8", 
                               "Crip1", "Ighm", "Mzb1", "Id3", "S1pr4", "Cxcr4", "Klf2", "Ccr7", "Cd69", "Fos", 
                               "mt-Cytb", "Gm42418", "AY036118", "Lars2", "Lgals3", "Lgals1", "Igkc", "Iglc3", "Iglc2", 
                               "Ltb", "Cd79b", "Ms4a1", "Fcer2a", "Foxo1", "Klf3", "Nfkbiz", "Ddx21", "Cd22", "Malat1", 
                               "Bach2", "Bank1", "Ighd", "Cd38")


rev_genes_to_confirm_manual_ass <- rev(genes_to_confirm_manual_ass)
rev_gene_list_for_publication <- rev(gene_list_for_publication)
Idents(scRNA_LTa_Bcells) <- "celltypes_manual"
DotPlot(scRNA_LTa_Bcells, features = rev_gene_list_for_publication) + theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "firebrick")
###############################################################################


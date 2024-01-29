options(stringsAsFactors = FALSE)
Sys.getenv('R_MAX_VSIZE')


library(sctransform)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(clusterProfiler)
library(patchwork)
library(KEGG.db)
library(DOSE)
library(org.Mm.eg.db)
library(colorRamps)
library(RColorBrewer)
library(tidyr)
library(SeuratWrappers)
library(doBy)
library(fgsea)
library(data.table)
library(pheatmap)
library(viridis)
library(gprofiler2)
library(scDblFinder)
library(SoupX)
library(lsa)
library(msigdbr)

###functions to use###############
pc_select <- function(seu){
  .seu <- seu
  .pct <- .seu[["pca"]]@stdev / sum(.seu[["pca"]]@stdev) * 100
  .cumu <- cumsum(.pct)
  .co1 <- which(.cumu > 90 & .pct < 5)[1]
  .co2 <- sort(which((.pct[1:length(.pct) - 1] - .pct[2:length(.pct)]) > 0.1), decreasing = T)[1] + 1
  .pc_num <- min(.co1, .co2)
  return(.pc_num)
}



#####Integration with soupX decontamination########################
soupx_lib1 <- load10X("cellranger_output_for_soupX/lib1_outs/")
soupx_lib1 <- autoEstCont(soupx_lib1)
soupx_lib1 <- adjustCounts(soupx_lib1)

set.seed(1234)
lib1 <- CreateSeuratObject(counts = soupx_lib1, project = "lib1")
dblet <- scDblFinder(GetAssayData(lib1, slot="counts"))
lib1$scDblFinder.class <- dblet$scDblFinder.class
lib1##13102 cells
lib1$condition <- "GFP_CFC"
lib1[["percent.mt"]] <- PercentageFeatureSet(lib1, pattern = "^mt-")
lib1 <- RenameCells(lib1, add.cell.id = "lib1")

soupx_lib2 <- load10X("cellranger_output_for_soupX/lib2_outs/")
soupx_lib2 <- autoEstCont(soupx_lib2)
soupx_lib2 <- adjustCounts(soupx_lib2)

set.seed(1234)
lib2 <- CreateSeuratObject(counts = soupx_lib2, project = "lib2")
dblet <- scDblFinder(GetAssayData(lib2, slot="counts"))
lib2$scDblFinder.class <- dblet$scDblFinder.class
lib2##11792 cells
lib2$condition <- "CRE_CFC"
lib2[["percent.mt"]] <- PercentageFeatureSet(lib2, pattern = "^mt-")
lib2 <- RenameCells(lib2, add.cell.id = "lib2")

soupx_lib3 <- load10X("cellranger_output_for_soupX/lib3_outs/")
soupx_lib3 <- autoEstCont(soupx_lib3)
soupx_lib3 <- adjustCounts(soupx_lib3)

set.seed(1234)
lib3 <- CreateSeuratObject(counts = soupx_lib3, project = "lib3")
dblet <- scDblFinder(GetAssayData(lib3, slot="counts"))
lib3$scDblFinder.class <- dblet$scDblFinder.class
lib3##11197 cells
lib3$condition <- "GFP_none"
lib3[["percent.mt"]] <- PercentageFeatureSet(lib3, pattern = "^mt-")
lib3 <- RenameCells(lib3, add.cell.id = "lib3")


soupx_lib4 <- load10X("cellranger_output_for_soupX/lib4_outs/")
soupx_lib4 <- autoEstCont(soupx_lib4)
soupx_lib4 <- adjustCounts(soupx_lib4)


set.seed(1234)
lib4 <- CreateSeuratObject(counts = soupx_lib4, project = "lib4")
dblet <- scDblFinder(GetAssayData(lib4, slot="counts"))
lib4$scDblFinder.class <- dblet$scDblFinder.class
lib4##16301 cells
lib4$condition <- "CRE_none"
lib4[["percent.mt"]] <- PercentageFeatureSet(lib4, pattern = "^mt-")
lib4 <- RenameCells(lib4, add.cell.id = "lib4")


soupx_lib5 <- load10X("cellranger_output_for_soupX/lib5_outs/")
#soupx_lib5 <- autoEstCont(soupx_lib5)
soupx_lib5 = setContaminationFraction(soupx_lib5, 0.2)
soupx_lib5 <- adjustCounts(soupx_lib5)

set.seed(1234)
lib5 <- CreateSeuratObject(counts = soupx_lib5, project = "lib5")
dblet <- scDblFinder(GetAssayData(lib5, slot="counts"))
lib5$scDblFinder.class <- dblet$scDblFinder.class
lib5##9950 cells
lib5$condition <- "CRE_none"
lib5[["percent.mt"]] <- PercentageFeatureSet(lib5, pattern = "^mt-")
lib5 <- RenameCells(lib5, add.cell.id = "lib5")




####filtration############################
lib1 <- subset(lib1, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 &percent.mt < 5 & scDblFinder.class == "singlet")##8254
lib2 <- subset(lib2, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & scDblFinder.class == "singlet")##7568
lib2
lib3 <- subset(lib3, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & scDblFinder.class == "singlet")##8147
lib3
lib4 <- subset(lib4, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & scDblFinder.class == "singlet")##12210
lib4
lib5 <- subset(lib5, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 5 & scDblFinder.class == "singlet")##7962
lib5




list3 <- c(lib1,lib2,lib3,lib4,lib5)

list3 <- lapply(X = list3, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, verbose = FALSE)
  x <- SCTransform(x, vst.flavor = "v2", verbose = TRUE) %>% RunPCA(npcs = 50, verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = list3, nfeatures = 3000)
list3 <- PrepSCTIntegration(object.list = list3, anchor.features = features)


TLR9.fourSamples.withDbltRemoved.withSoupX.anchors <- FindIntegrationAnchors(object.list = list3, normalization.method = "SCT",
                                                                             anchor.features = features)

TLR9.fourSamples.withDbltRemoved.withSoupX.sct <- IntegrateData(anchorset = TLR9.fourSamples.withDbltRemoved.withSoupX.anchors, normalization.method = "SCT")

saveRDS(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, file = "TLR9.fourSamples.withDbltRemoved.withSoupX.sct.afterIntegrateData.rds")
TLR9.fourSamples.withDbltRemoved.withSoupX.sct <- readRDS("TLR9.fourSamples.withDbltRemoved.withSoupX.sct.afterIntegrateData.rds")


TLR9.fourSamples.withDbltRemoved.withSoupX.sct <- RunPCA(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, npcs = 50, verbose = FALSE)

pc_num <- pc_select(TLR9.fourSamples.withDbltRemoved.withSoupX.sct)

TLR9.fourSamples.withDbltRemoved.withSoupX.sct <- TLR9.fourSamples.withDbltRemoved.withSoupX.sct %>%
  RunUMAP(reduction = "pca", dims = 1:pc_num, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:pc_num) %>%
  FindClusters(resolution = 0.5)

saveRDS(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, file = "TLR9.fourSamples.withDbltRemoved.withSoupX.sct.Final.rds")
TLR9.fourSamples.withDbltRemoved.withSoupX.sct <- readRDS("TLR9.fourSamples.withDbltRemoved.withSoupX.sct.Final.rds")


####Cluster markers######################################
cluster.ALLmarkers <- FindAllMarkers(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, assay = "RNA", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

ordered.ALLmarkers <- cluster.ALLmarkers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)

write.table(ordered.ALLmarkers, "soupX_corrected_analysis/cluster_ALLmarkers.csv", sep = ",",row.names = F)



filtered.markers <- cluster.ALLmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 25, order_by = avg_log2FC)

write.table(filtered.markers, "soupX_corrected_analysis/cluster_ALLmarkers_top20.csv", sep = ",",row.names = F)




####Differential analysis#######################################
DefaultAssay(TLR9.fourSamples.withDbltRemoved.withSoupX.sct) <- "RNA"
TLR9.fourSamples.withDbltRemoved.withSoupX.sct.GFPpos <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, subset= `EGFP-gn` > 0)


table(TLR9.fourSamples.withDbltRemoved.withSoupX.sct$orig.ident,TLR9.fourSamples.withDbltRemoved.withSoupX.sct$seurat_clusters)
table(TLR9.fourSamples.withDbltRemoved.withSoupX.sct.GFPpos$orig.ident,TLR9.fourSamples.withDbltRemoved.withSoupX.sct.GFPpos$seurat_clusters)



for (clst in c(0:4,6:8,10:29)) {
  
  std.combined.cluster <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct,ident=clst)
  
  Idents(std.combined.cluster) <- "orig.ident"
  
  if (sum(std.combined.cluster$orig.ident == "lib1") > 3 & sum(std.combined.cluster$orig.ident == "lib3") > 3){
    cluster.marker.by.condition1 <- FindMarkers(std.combined.cluster,  assay = "RNA",ident.1 = "lib1", ident.2 = "lib3",verbose = TRUE,
                                                logfc.threshold = 0,min.pct = 0.01)
  }
  
  if (sum(std.combined.cluster$orig.ident == "lib2") > 3 & sum(std.combined.cluster$orig.ident == "lib1") > 3){
    cluster.marker.by.condition2 <- FindMarkers(std.combined.cluster,  assay = "RNA",ident.1 = "lib2", ident.2 = "lib1",verbose = TRUE,
                                                logfc.threshold = 0,min.pct = 0.01)
  }
  
  if (sum(std.combined.cluster$orig.ident == "lib4") > 3 & sum(std.combined.cluster$orig.ident == "lib3") > 3){
    cluster.marker.by.condition3 <- FindMarkers(std.combined.cluster,  assay = "RNA",ident.1 = "lib4", ident.2 = "lib3",verbose = TRUE,
                                                logfc.threshold = 0,min.pct = 0.01)
  }
  
  
  write.csv(cluster.marker.by.condition1, paste0("soupX_corrected_analysis/differential_analysis/cluster_",clst,"_Lib1vsLib3_full_diff.csv"))
  write.csv(cluster.marker.by.condition2, paste0("soupX_corrected_analysis/differential_analysis/cluster_",clst,"_Lib2vsLib1_full_diff.csv"))
  write.csv(cluster.marker.by.condition3, paste0("soupX_corrected_analysis/differential_analysis/cluster_",clst,"_Lib4vsLib3_full_diff.csv"))
  
  diff_sig1 <- cluster.marker.by.condition1 %>% filter(abs(avg_log2FC) > 0.585 & p_val_adj < 0.01)
  diff_sig2 <- cluster.marker.by.condition2 %>% filter(abs(avg_log2FC) > 0.585 & p_val_adj < 0.01)
  diff_sig3 <- cluster.marker.by.condition3 %>% filter(abs(avg_log2FC) > 0.585 & p_val_adj < 0.01)
  
  if (length(diff_sig1$p_val) != 0){
    write.csv(diff_sig1,paste0("soupX_corrected_analysis/differential_analysis/cluster_",clst,"_Lib1vsLib3_SIG_diff.csv"),row.names = T)
  }
  
  if (length(diff_sig2$p_val) != 0){
    write.csv(diff_sig2,paste0("soupX_corrected_analysis/differential_analysis/cluster_",clst,"_Lib2vsLib1_SIG_diff.csv"),row.names = T)
  }
  
  if (length(diff_sig3$p_val) != 0){
    write.csv(diff_sig3,paste0("soupX_corrected_analysis/differential_analysis/cluster_",clst,"_Lib4vsLib3_SIG_diff.csv"),row.names = T)
  }
  
  print(paste0("Finished for cluster ",clst))
  
}





for (clst in c(0:4,6:8,10:29)) {
  
  std.combined.cluster <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct.GFPpos,ident=clst)
  
  Idents(std.combined.cluster) <- "orig.ident"
  
  if (sum(std.combined.cluster$orig.ident == "lib1") > 3 & sum(std.combined.cluster$orig.ident == "lib3") > 3){
    cluster.marker.by.condition1 <- FindMarkers(std.combined.cluster,  assay = "RNA",ident.1 = "lib1", ident.2 = "lib3",verbose = TRUE,
                                                logfc.threshold = 0,min.pct = 0.01)
  }
  
  if (sum(std.combined.cluster$orig.ident == "lib2") > 3 & sum(std.combined.cluster$orig.ident == "lib1") > 3){
    cluster.marker.by.condition2 <- FindMarkers(std.combined.cluster,  assay = "RNA",ident.1 = "lib2", ident.2 = "lib1",verbose = TRUE,
                                                logfc.threshold = 0,min.pct = 0.01)
  }
  
  if (sum(std.combined.cluster$orig.ident == "lib4") > 3 & sum(std.combined.cluster$orig.ident == "lib3") > 3){
    cluster.marker.by.condition3 <- FindMarkers(std.combined.cluster,  assay = "RNA",ident.1 = "lib4", ident.2 = "lib3",verbose = TRUE,
                                                logfc.threshold = 0,min.pct = 0.01)
  }
  
  
  write.csv(cluster.marker.by.condition1, paste0("soupX_corrected_analysis/differential_analysis_GFPpos/cluster_",clst,"_Lib1vsLib3_full_diff.csv"))
  write.csv(cluster.marker.by.condition2, paste0("soupX_corrected_analysis/differential_analysis_GFPpos/cluster_",clst,"_Lib2vsLib1_full_diff.csv"))
  write.csv(cluster.marker.by.condition3, paste0("soupX_corrected_analysis/differential_analysis_GFPpos/cluster_",clst,"_Lib4vsLib3_full_diff.csv"))
  
  diff_sig1 <- cluster.marker.by.condition1 %>% filter(abs(avg_log2FC) > 0.585 & p_val_adj < 0.1)
  diff_sig2 <- cluster.marker.by.condition2 %>% filter(abs(avg_log2FC) > 0.585 & p_val_adj < 0.1)
  diff_sig3 <- cluster.marker.by.condition3 %>% filter(abs(avg_log2FC) > 0.585 & p_val_adj < 0.1)
  
  if (length(diff_sig1$p_val) != 0){
    write.csv(diff_sig1,paste0("soupX_corrected_analysis/differential_analysis_GFPpos/cluster_",clst,"_Lib1vsLib3_SIG_diff.csv"),row.names = T)
  }
  
  if (length(diff_sig2$p_val) != 0){
    write.csv(diff_sig2,paste0("soupX_corrected_analysis/differential_analysis_GFPpos/cluster_",clst,"_Lib2vsLib1_SIG_diff.csv"),row.names = T)
  }
  
  if (length(diff_sig3$p_val) != 0){
    write.csv(diff_sig3,paste0("soupX_corrected_analysis/differential_analysis_GFPpos/cluster_",clst,"_Lib4vsLib3_SIG_diff.csv"),row.names = T)
  }
  
  print(paste0("Finished for cluster ",clst))
  
}



####Pathway analysis#####################################
######GO#####
msigdbr_go <- msigdbr("Mus musculus", "C5")
msigdbr_go <- msigdbr_go %>% filter(gs_subcat %in% c("GO:BP","GO:CC","GO:MF"))

msigdbr_list2 = split(x = msigdbr_go$gene_symbol, f = msigdbr_go$gs_name)

for (path in c(1:length(msigdbr_list2))){
  msigdbr_list2[[path]] <- unique(msigdbr_list2[[path]])
  
}


for (clst in c(0:4,6:8,10:29)) {
  diff_table_combined <- read.csv(paste0("soupX_corrected_analysis/differential_analysis_allCells/cluster_",clst,"_Lib2vsLib1_full_diff.csv"))
  
  ##GESA pathway analysis
  ranks <- diff_table_combined$avg_log2FC
  names(ranks) <- diff_table_combined$X
  
  fgseaRes1 <- fgsea(msigdbr_list2, ranks, minSize=15, maxSize = 500)
  
  for (n in c(1:nrow(fgseaRes1))){
    fgseaRes1$leadingEdgeSize[n] <- length(fgseaRes1[,8][[1]][[n]])
  }
  fgseaRes1<- fgseaRes1[,c(1:7,9,8)] %>% arrange(pval)
  fwrite(fgseaRes1,file = paste0("soupX_corrected_analysis/pathway_analysis_allCells/cluster_",clst,"_Lib2vsLib1_GSEA_FULL_GO_pathway.csv"),sep = ",",sep2 = c(""," ",""))
  
  gsea_output1 <- fgseaRes1[which(fgseaRes1$padj < 0.05)][order(pval, -abs(NES)), ]
  
  ##print out significant regulated pathways
  if (length(gsea_output1$pathway) != 0){
    fwrite(gsea_output1, paste0("soupX_corrected_analysis/pathway_analysis_allCells/cluster_",clst,"_Lib2vsLib1_GSEA_SIG_GO_pathway.csv"),sep = ",",sep2 = c(""," ",""))
  }
  
  print(paste0("Finished for cluster ",clst))
}




for (clst in c(0:4,6:8,10:29)) {
  
  diff_table_combined <- read.csv(paste0("soupX_corrected_analysis/differential_analysis_allCells/cluster_",clst,"_Lib1vsLib3_full_diff.csv"))
  
  ##GESA pathway analysis
  ranks <- diff_table_combined$avg_log2FC
  names(ranks) <- diff_table_combined$X
  
  fgseaRes1 <- fgsea(msigdbr_list2, ranks, minSize=15, maxSize = 500)
  
  for (n in c(1:nrow(fgseaRes1))){
    fgseaRes1$leadingEdgeSize[n] <- length(fgseaRes1[,8][[1]][[n]])
  }
  fgseaRes1<- fgseaRes1[,c(1:7,9,8)] %>% arrange(pval)
  fwrite(fgseaRes1,file = paste0("soupX_corrected_analysis/pathway_analysis_allCells/cluster_",clst,"_Lib1vsLib3_GSEA_FULL_GO_pathway.csv"),sep = ",",sep2 = c(""," ",""))
  
  gsea_output1 <- fgseaRes1[which(fgseaRes1$padj < 0.05)][order(pval, -abs(NES)), ]
  
  ##print out significant regulated pathways
  if (length(gsea_output1$pathway) != 0){
    fwrite(gsea_output1, paste0("soupX_corrected_analysis/pathway_analysis_allCells/cluster_",clst,"_Lib1vsLib3_GSEA_SIG_GO_pathway.csv"),sep = ",",sep2 = c(""," ",""))
  }
  
  print(paste0("Finished for cluster ",clst))
}



####plots and tables#########
######Stats for up/down regulated genes per cluster###################################
table <- data.frame()
for (clst in c(0:4,6:8,10:29)) {
  
  if(file.exists(paste0("soupX_corrected_analysis/differential_analysis_allCells/cluster_",clst,"_Lib1vsLib3_SIG_diff.csv"))){
    diff1 <- read.csv(paste0("soupX_corrected_analysis/differential_analysis_allCells/cluster_",clst,"_Lib1vsLib3_SIG_diff.csv"))
    
    up1 <- length(diff1 %>% filter(avg_log2FC > 0) %>% pull(X))
    down1 <- length(diff1 %>% filter(avg_log2FC < 0) %>% pull(X))
    
    table <- rbind(table,c(clst,up1,down1))
  }
  
  else{
    table <- rbind(table,c(clst,0,0))
  }
  
  print(paste0("Finished for cluster ",clst))
}


table2 <- data.frame()
for (clst in c(0:4,6:8,10:29)) {
  
  if(file.exists(paste0("soupX_corrected_analysis/differential_analysis_allCells/cluster_",clst,"_Lib2vsLib1_SIG_diff.csv"))){
    diff1 <- read.csv(paste0("soupX_corrected_analysis/differential_analysis_allCells/cluster_",clst,"_Lib2vsLib1_SIG_diff.csv"))
    
    up1 <- length(diff1 %>% filter(avg_log2FC > 0) %>% pull(X))
    down1 <- length(diff1 %>% filter(avg_log2FC < 0) %>% pull(X))
    
    table2 <- rbind(table2,c(clst,up1,down1))
  }
  
  else{
    table2 <- rbind(table2,c(clst,0,0))
  }
  
  print(paste0("Finished for cluster ",clst))
}

colnames(table) <- c("cluster","up", "down")
colnames(table2) <- c("cluster","up", "down")


write.csv(table,"soupX_corrected_analysis/lib1vs3_significant_diff_stats.csv")
write.csv(table2,"soupX_corrected_analysis/lib2vs1_significant_diff_stats.csv")

######Violin plots of Hsp90b1################################
lib1 <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, subset=orig.ident == "lib1")
lib2 <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, subset=orig.ident == "lib2")
lib3 <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, subset=orig.ident == "lib3")

lib1 <- subset(lib1, subset=seurat_clusters != "9")
lib2 <- subset(lib2, subset=seurat_clusters != "9")
lib3 <- subset(lib3, subset=seurat_clusters != "9")

p1 <- VlnPlot(lib3, features = c("Hsp90b1"), pt.size = 0) + BoldTitle() + ylab("lib3") + ggtitle("Hsp90b1") + SeuratAxes() +
  theme(legend.position = "none",axis.title.y=element_text(angle=0,vjust=0.5), axis.title.x = element_blank(),axis.text.x=element_blank())
p2 <- VlnPlot(lib1, features = c("Hsp90b1"), pt.size = 0) + BoldTitle() + ylab("lib1") + SeuratAxes() +
  theme(legend.position = "none",axis.title.x = element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),plot.title = element_blank())
p3 <- VlnPlot(lib2, features = c("Hsp90b1"), pt.size = 0) + SeuratAxes() + ylab("lib2") + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),plot.title = element_blank())


pdf("soupX_corrected_analysis/Vlnplot_Hsp90b1.pdf",width = 12, height = 8)
wrap_plots(p1,p2,p3,ncol=1)
dev.off()


z1 <- FeaturePlot(lib1, features = "Hsp90b1") + ggtitle("Lib1")
z2 <- FeaturePlot(lib2, features = "Hsp90b1")+ ggtitle("Lib2")
z3 <- FeaturePlot(lib3, features = "Hsp90b1")+ ggtitle("Lib3")


pdf("soupX_corrected_analysis/FeaturePlot_Hsp90b1.pdf",width = 12, height = 8)
wrap_plots(z3,z1,z2,ncol=3)
dev.off()


###split violinplots
cluster19 <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, subset=seurat_clusters == "19")
cluster20 <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, subset=seurat_clusters == "20")

cluster19_sample123 <- subset(cluster19, subset=orig.ident %in% c("lib1","lib2","lib3"))
cluster20_sample123 <- subset(cluster20, subset=orig.ident %in% c("lib1","lib2","lib3"))

table(cluster19$seurat_clusters)
table(cluster19$orig.ident)
table(cluster19_sample123$orig.ident)
table(cluster20_sample123$orig.ident)


cluster19_sample123$orig.ident <- factor(cluster19_sample123$orig.ident,levels = c("lib3","lib1","lib2"))
cluster20_sample123$orig.ident <- factor(cluster20_sample123$orig.ident,levels = c("lib3","lib1","lib2"))

pdf("soupX_corrected_analysis/Hsp90b1_plots/cluster19_vlnPlot_Hsp90b1.pdf",width = 8, height = 6)
VlnPlot(cluster19_sample123, group.by = "orig.ident",features = c("Hsp90b1"), pt.size = 0.5)
dev.off()

pdf("soupX_corrected_analysis/Hsp90b1_plots/cluster20_vlnPlot_Hsp90b1.pdf",width = 8, height = 6)
VlnPlot(cluster20_sample123, group.by = "orig.ident",features = c("Hsp90b1"), pt.size = 0.5)
dev.off()

dittoPlot(cluster19_sample123, "Hsp90b1", group.by = "orig.ident",
          plots = c("vlnplot", "jitter"))

dittoBoxPlot(cluster19_sample123, "Hsp90b1", group.by = "orig.ident")


###heatmap on selected genes
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install dittoSeq
BiocManager::install("dittoSeq")
library(dittoSeq)

gelist <- c("Hsp90b1","Hspa5","Atp6v0c","Bsg","Cck","Itm2b","Grina","Pcsk1n","Ly6h","Apoe","Ttr","Cdh13","Galntl6","Unc5d","Nrg1","Pcdh15","Sorcs1")

pdf("soupX_corrected_analysis/Hsp90b1_plots/cluster19_heatmap.pdf",width = 8, height = 6)
dittoHeatmap(cluster19_sample123, assay = "RNA", slot = "data", gelist,group.by = "orig.ident",annot.by = "orig.ident")
dev.off()

pdf("soupX_corrected_analysis/Hsp90b1_plots/cluster20_heatmap.pdf",width = 8, height = 6)
dittoHeatmap(cluster20_sample123, assay = "RNA", slot = "data", gelist,group.by = "orig.ident",annot.by = "orig.ident")
dev.off()



######correlation between Hsp90b1 and Atp6v0c########################
cluster19_sample123$Hsp90b1 <- cluster19_sample123@assays$RNA@data["Hsp90b1",]
cluster19_sample123$Atp6v0c <- cluster19_sample123@assays$RNA@data["Atp6v0c",]

cluster20_sample123$Hsp90b1 <- cluster20_sample123@assays$RNA@data["Hsp90b1",]
cluster20_sample123$Atp6v0c <- cluster20_sample123@assays$RNA@data["Atp6v0c",]

test <- data.frame(cluster19_sample123$Hsp90b1,cluster19_sample123$Atp6v0c)

cosine(cluster19_sample123$Hsp90b1,cluster19_sample123$Atp6v0c)##0.608
cosine(cluster20_sample123$Hsp90b1,cluster20_sample123$Atp6v0c)##0.591

###cosine similarity and permutation test 
cs_list <- c()
gene_pool <- read.csv("soupX_corrected_analysis/differential_analysis_allCells/cluster_19_Lib1vsLib3_full_diff.csv")
gene_pool <- data.frame(gene_pool$X)
for (i in c(1:10000)){
  rand <- sample(1:15745,2)
  rand1 <- rand[1]
  rand2 <- rand[2]
  
  cs <- cosine(cluster19_sample123@assays$RNA@data[gene_pool[rand1,],],cluster19_sample123@assays$RNA@data[gene_pool[rand2,],])
  
  cs_list <- c(cs_list, cs)
}

table(cs_list > 0.608)


cs_list <- c()
gene_pool <- read.csv("soupX_corrected_analysis/differential_analysis_allCells/cluster_20_Lib1vsLib3_full_diff.csv")
gene_pool <- data.frame(gene_pool$X)
for (i in c(1:10000)){
  rand <- sample(1:15745,2)
  rand1 <- rand[1]
  rand2 <- rand[2]
  
  cs <- cosine(cluster20_sample123@assays$RNA@data[gene_pool[rand1,],],cluster20_sample123@assays$RNA@data[gene_pool[rand2,],])
  
  cs_list <- c(cs_list, cs)
}

table(cs_list > 0.591)#0.008

table(cluster19_sample123$Hsp90b1 > 0)##0.45
table(cluster19_sample123$Atp6v0c > 0)##0.41
table(cluster19_sample123$Hsp90b1 > 0 & cluster19_sample123$Atp6v0c > 0)##0.26


table(cluster20_sample123$Hsp90b1 > 0)##0.44
table(cluster20_sample123$Atp6v0c > 0)##0.39
table(cluster20_sample123$Hsp90b1 > 0 & cluster20_sample123$Atp6v0c > 0)##0.23

cosine(cluster20_sample123$Hsp90b1,cluster20_sample123$Atp6v0c)##0.59


######Rad51,Trp53bp1,Cetn2#####################################

#########lib4 vs lib3##################
lib3_4 <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, subset=orig.ident %in% c("lib3","lib4"))


##Rad51, Trp53BP1, and centrin2
lib3_4$Rad51 <- lib3_4@assays$RNA@data["Rad51",]
lib3_4$Trp53bp1 <- lib3_4@assays$RNA@data["Trp53bp1",]
lib3_4$Cetn2 <- lib3_4@assays$RNA@data["Cetn2",]

data.summary <- lib3_4@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    Rad51 = mean(Rad51),
    Trp53bp1 = mean(Trp53bp1),
    Cetn2 = mean(Cetn2)
  )


DefaultAssay(lib3_4) <- "RNA"

z1 <- VlnPlot(lib3_4, features = ("Rad51"), split.by = "orig.ident")
z2 <- VlnPlot(lib3_4, features = ("Trp53bp1"), split.by = "orig.ident")
z3 <- VlnPlot(lib3_4, features = ("Cetn2"), split.by = "orig.ident")

pdf("soupX_corrected_analysis/libr3_vs_lib4/VlnPlot_Rad51_Trp53bp1_Cetn2.pdf",width = 12, height = 6)
wrap_plots(z1,z2,z3,ncol=1)
dev.off()



#######lib2 vs lib1##################
lib1_2 <- subset(TLR9.fourSamples.withDbltRemoved.withSoupX.sct, subset=orig.ident %in% c("lib1","lib2"))


##Rad51, Trp53BP1, and centrin2
lib1_2$Rad51 <- lib1_2@assays$RNA@data["Rad51",]
lib1_2$Trp53bp1 <- lib1_2@assays$RNA@data["Trp53bp1",]
lib1_2$Cetn2 <- lib1_2@assays$RNA@data["Cetn2",]

data.summary <- lib1_2@meta.data %>%
  group_by(seurat_clusters,orig.ident) %>%
  summarise(
    Rad51 = mean(Rad51),
    Trp53bp1 = mean(Trp53bp1),
    Cetn2 = mean(Cetn2)
  )






z1 <- VlnPlot(lib1_2, features = ("Rad51"), split.by = "orig.ident") + ylim(-0.1,2.5)
z2 <- VlnPlot(lib1_2, features = ("Trp53bp1"), split.by = "orig.ident")+ ylim(-0.1,2.5)
z3 <- VlnPlot(lib1_2, features = ("Cetn2"), split.by = "orig.ident")+ ylim(-0.1,2.5)

p1 <- VlnPlot(lib1_2, features = c("Rad51"),split.by = "orig.ident") + BoldTitle() + ylab("Rad51")  + SeuratAxes() + ylim(-0.1,2.5) +
  theme(legend.position = "top",axis.title.x = element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),plot.title = element_blank())+ scale_fill_manual(values=c("blue", "green3"))
p2 <- VlnPlot(lib1_2, features = c("Trp53bp1"),split.by = "orig.ident") + BoldTitle() + ylab("Trp53bp1") + SeuratAxes() + ylim(-0.1,2.5) +
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),plot.title = element_blank())+ scale_fill_manual(values=c("blue", "green3"))
p3 <- VlnPlot(lib1_2, features = c("Cetn2"),split.by = "orig.ident") + BoldTitle() + ylab("Cetn2") + SeuratAxes() + ylim(-0.1,2.5) +
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),plot.title = element_blank())+ scale_fill_manual(values=c("blue", "green3"))





pdf("soupX_corrected_analysis/libr3_vs_lib4/lib1vs2_VlnPlot_Rad51_Trp53bp1_Cetn2_v2.pdf",width = 12, height = 6)
wrap_plots(p1,p2,p3,ncol=1)
dev.off()


my_table <- data.frame()
for (clst in c(0:4,6:8,10:29)) {
  
  diff_table <- read.csv(paste0("soupX_corrected_analysis/differential_analysis_allCells/cluster_",clst,"_Lib2vsLib1_full_diff.csv"))
  
  diff_table <- diff_table %>% filter(X %in% c("Rad51","Trp53bp1","Cetn2"))
  diff_table$cluster <- clst
  
  my_table <- rbind(my_table,diff_table)
}

write.csv(my_table, "soupX_corrected_analysis/libr3_vs_lib4/Lib2vsLib1_Rad51_Trp53bp1_Cetn2_diff.csv")

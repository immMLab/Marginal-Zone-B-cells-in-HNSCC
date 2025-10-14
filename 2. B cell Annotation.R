# subset B cells -------------
setwd('~/Project/5Dataset')
options(warn=-1) # turn off warning message globally
rm(list=ls());gc()
options(stringsAsFactors = F) 
options(as.is = T)
width.ppi=6.5
height.ppi=5


scRNA <- readRDS('3-CellType/sce.all.singleR.rds')

dir.create("4-Bcells")
setwd("4-Bcells")
table(scRNA$celltype_DICE)
Bcell.cut <- subset(scRNA,subset=celltype_DICE=="B cells, naive")
dim(Bcell.cut)   # 37488 31310
table(Bcell.cut$orig.ident)
Bcell.cut <- SplitObject(Bcell.cut, split.by="dataset")
Bcell.list.sct <- lapply(X = Bcell.cut, FUN = function(x) {
  SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})
Bcell.list.pca <- lapply(X = Bcell.list.sct, FUN = function(y) {
  RunPCA(y, npcs = 30, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = Bcell.list.pca, nfeatures = 3000)
Bcell.list.pca <- PrepSCTIntegration(object.list = Bcell.list.pca, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = Bcell.list.pca, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.3)


p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "dataset")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
              repel = TRUE)
p3 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "orig.ident", label = TRUE,
              repel = TRUE) +NoLegend()
p1 | p2 | p3

Bcell = immune.combined.sct


p1 <- DimPlot(Bcell, group.by = "seurat_clusters", raster=FALSE);p1
ggsave('Bcell_dimplot.png',p1,height = height.ppi,width = width.ppi)
p2 <- DimPlot(Bcell, group.by = "orig.ident", raster=FALSE,repel = FALSE) +NoLegend();p2
ggsave('Bcell_no_batch_effect.png',p2,height = height.ppi,width = width.ppi)



###### B cell features ######
Bcell <- AddModuleScore(Bcell, list(c("IGHA1", "IGHA2")), name='IgA.RNA', assay = "RNA")
Bcell <- AddModuleScore(Bcell, list(c("IGHG1", "IGHG2", "IGHG3", "IGHG4")), name='IgG.RNA', assay = "RNA")
Bcell <- AddModuleScore(Bcell, list(c('MME', 'BCL6', 'BCL7A', "PCNA")), name='GC.RNA', assay = "RNA")
Bcell <- AddModuleScore(Bcell, list(c("CD1C", "PLD4",'CDH17', 'DOCK11', 'PTK2B', 'DLL1', 'LFNG', 'MFNG', 'NOTCH2', 'DOCK10', 'BCL3')), name='MZB.RNA', assay = "RNA") #"CD27", 'IGHM',"CR2"


rna_features_output_split_dot <- c("IGHD", "IGHM", "CD27", "CD1C", "PLD4", "CD5","CR2","FCER2", "CCR7", "MME", "MZB1", "IGLL5", "CD24", "CD38", "COCH", "CXCR5", "ZEB2", "MX1", "BCL6", "BCL7A", "PCNA","MZB.RNA1", "GC.RNA1", "IgG.RNA1", "IgA.RNA1")

plot_alldot <- DotPlot(Bcell, features = rna_features_output_split_dot, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"));plot_alldot
ggsave("dotplot_B_celltype2.pdf", plot_alldot, width=width.ppi*1.5 ,height=height.ppi)  # -> find marker genes to compose Figure 2b

saveRDS(Bcell,'Bcell.rds')

sce = Bcell
sce <- FindClusters(object = sce,resolution = c(seq(from=.1,to=1.6,by=.2))) 
colnames(sce@meta.data)
p<-clustree(sce@meta.data, prefix = "integrated_snn_res.");p
ggsave("Bcell_clustertree.png",p,width = width.ppi*2,height = height.ppi*2)
#clustree(sce, prefix = "SCT_snn_res.",node_colour = "CD1C", node_colour_aggr = "median")

tab.1=table(sce@meta.data$SCT_snn_res.0.5,sce@meta.data$SCT_snn_res.0.8)



sce <- SetIdent(sce, value = 'integrated_snn_res.0.3')
p1 <- DotPlot(sce, features = rna_features_output_split_dot, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"));p1
p2=DimPlot(sce, reduction = "umap",label=T,raster=FALSE);p2
ggsave("Bcell_0.3_markers.png",p1,width = width.ppi,height = height.ppi)
ggsave("Bcell_resolution_0.3.png",p2,width = width.ppi,height = height.ppi)

sce <- SetIdent(sce, value = 'integrated_snn_res.0.5')
p3 <- DotPlot(sce, features = rna_features_output_split_dot, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"));p3
p4=DimPlot(sce, reduction = "umap",label=T,raster=FALSE);p4
ggsave("Bcell_0.5_markers.png",p3,width = width.ppi,height = height.ppi)
ggsave("Bcell_resolution_0.5.png",p4,width = width.ppi,height = height.ppi)

sce <- SetIdent(sce, value = 'integrated_snn_res.0.1')
p5 <- DotPlot(sce, features = rna_features_output_split_dot, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"));p5
p6=DimPlot(sce, reduction = "umap",label=T,raster=FALSE);p6
ggsave("Bcell_0.1_markers.png",p5,width = width.ppi,height = height.ppi)
ggsave("Bcell_resolution_0.1.png",p6,width = width.ppi,height = height.ppi)

saveRDS(sce,'Bcell.resolution0.1~1.3.rds')

#### SingleR ####
setwd("../")
dir.create("5-annotationB")
setwd("5-annotationB")
library(SingleR)
library(celldex)

##set the width and height of plots
width.ppi=6.5
height.ppi=5


## subsets color
Idents(Bcell) <- 'celltype'
cell_type_cols <- c(brewer.pal(12, "Set3"), 
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")  
COLORS_NAMES <- levels(Bcell$celltype)
names(cell_type_cols) <- COLORS_NAMES

DimPlot(Bcell, label = TRUE, repel = TRUE,cols = cell_type_cols)

##### HPCA #####

refdata <- HumanPrimaryCellAtlasData()
Bcell <- readRDS('../4-Bcells/Bcell.resolution0.1~1.3.rds')

#refdata$ {"label.fine"}
testdata <- GetAssayData(Bcell, slot="data")
clusters <- Bcell@meta.data$integrated_snn_res.0.1

cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellpred$pruned.labels
a <- plotScoreHeatmap(cellpred, clusters=cellpred@rownames, fontsize.row = 9,show_colnames = T)
ggsave("HPCA_heatmap.png", a, width=7 ,height=3)


celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"celltype_HPCA.csv",row.names = F)
for(i in 1:nrow(celltype)){
  Bcell@meta.data[which(Bcell@meta.data$integrated_snn_res.0.1 == celltype$ClusterID[i]),'celltype_HPCA'] <- celltype$celltype[i]}

a1 <- UMAPPlot(object = Bcell, group.by = "celltype_HPCA", pt.size = 0.3, repel=T, label=T, label.size=5);a1
ggsave("UMAP_celltype_HPCA.png", a, width=10 ,height=6)


#####  BlueprintEncodeData  #####
refdata <- BlueprintEncodeData()
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellpred$pruned.labels
a <- plotScoreHeatmap(cellpred, clusters=cellpred@rownames, fontsize.row = 9,show_colnames = T)
ggsave("BlueprintEncodeData_heatmap.png", a, width=7 ,height=3)

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"BlueprintEncodeData.csv",row.names = F)
for(i in 1:nrow(celltype)){
  Bcell@meta.data[which(Bcell@meta.data$integrated_snn_res.0.1 == celltype$ClusterID[i]),'celltype_BlueprintEncodeData'] <- celltype$celltype[i]}

a2 <- UMAPPlot(object = Bcell, group.by = "celltype_BlueprintEncodeData", pt.size = 0.3, repel=T, label=T, label.size=5);a2
ggsave("UMAP_celltype_HPCA.png", a, width=10 ,height=6)

a1+NoLegend()|a2+NoLegend()
ggsave("UMAP_celltype_HPCA_BPE.png", width=10 ,height=6)

#### GSE193868 Reference ####

#### Tissue
Tissue <- readRDS("~/Project/GSE193868-MZB two subsets/plots/Tissue/10xtissue.annotated.rds")

refdata1 <- Seurat::as.SingleCellExperiment(Tissue, assay = c("integrated","RNA"))
colData(refdata1)$cell_type1 <- Tissue$subset_names
rowData(refdata1)$feature_symbol <- rownames(refdata1)
#refdata <- scater::logNormCounts(refdata)

testdata <- GetAssayData(Bcell[["SCT"]], slot="data")
clusters <- Bcell@meta.data$integrated_snn_res.0.3
cellpred <- SingleR(test = testdata, ref = refdata1, labels = refdata1$cell_type1,clusters = clusters 
                    ,assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellpred$pruned.labels #predictions$pruned.labels
sum(is.na(cellpred$pruned.labels))
cellpred$pruned.labels[which(is.na(cellpred$pruned.labels))] <- "ambiguous"
a <- plotScoreHeatmap(cellpred, clusters=cellpred@rownames, fontsize.row = 9,show_colnames = T)
#ggsave("SingleR_heatmap_10xBlood.png", a, width=7 ,height=3)

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"SingleR_GSE193868_Tissue.csv",row.names = F)
for(i in 1:nrow(celltype)){
  Bcell@meta.data[which(Bcell@meta.data$integrated_snn_res.0.3 == celltype$ClusterID[i]),'celltype_GSE193_Tissue'] <- celltype$celltype[i]}

a1 <- UMAPPlot(object = Bcell, group.by = "celltype_GSE193_Tissue", pt.size = 0.3, repel=T, label=T, label.size=5);a1
ggsave("UMAP_SingleR_GSE193_Tissue.png", a, width=10 ,height=6)



#### Blood 
Blood <- readRDS("~/Project/GSE193868-MZB two subsets/10xblood_Bcell.rds")
refdata2 <- Seurat::as.SingleCellExperiment(Blood, assay = c("integrated","RNA"))
colData(refdata2)$cell_type1 <- Blood$subset_names
rowData(refdata2)$feature_symbol <- rownames(refdata2)

cellpred <- SingleR(test = testdata, ref = refdata2, labels = refdata2$cell_type1,clusters = clusters 
                    ,assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellpred$pruned.labels #predictions$pruned.labels
sum(is.na(cellpred$pruned.labels))
cellpred$pruned.labels[which(is.na(cellpred$pruned.labels))] <- "ambiguous"

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"SingleR_GSE193868_Blood.csv",row.names = F)
for(i in 1:nrow(celltype)){
  Bcell@meta.data[which(Bcell@meta.data$integrated_snn_res.0.3 == celltype$ClusterID[i]),'celltype_GSE193_Blood'] <- celltype$celltype[i]}

a2 <- UMAPPlot(object = Bcell, group.by = "celltype_GSE193_Blood", pt.size = 0.3, repel=T, label=T, label.size=5);a2


a1|a2
ggsave("UMAP_SingleR_GSE193_Blood_Tissue.png", width=10 ,height=6)


#### resolution = 0.5 
clusters <- Bcell@meta.data$integrated_snn_res.0.5
cellpred1 <- SingleR(test = testdata, ref = refdata1, labels = refdata1$cell_type1,clusters = clusters 
                    ,assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred1), celltype=cellpred1$labels, stringsAsFactors = F)
# write.csv(celltype,"SingleR_GSE193868_Tissue.csv",row.names = F)
for(i in 1:nrow(celltype)){
  Bcell@meta.data[which(Bcell@meta.data$integrated_snn_res.0.5 == celltype$ClusterID[i]),'celltype_GSE193_Tissue_0.5'] <- celltype$celltype[i]}
a1 <- UMAPPlot(object = Bcell, group.by = "celltype_GSE193_Tissue_0.5", pt.size = 0.3, repel=T, label=T, label.size=5);a1


cellpred2 <- SingleR(test = testdata, ref = refdata2, labels = refdata2$cell_type1,clusters = clusters 
                    ,assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred2), celltype=cellpred2$labels, stringsAsFactors = F)
# write.csv(celltype,"SingleR_GSE193868_Blood.csv",row.names = F)
for(i in 1:nrow(celltype)){
  Bcell@meta.data[which(Bcell@meta.data$integrated_snn_res.0.5 == celltype$ClusterID[i]),'celltype_GSE193_Blood_0.5'] <- celltype$celltype[i]}

a2 <- UMAPPlot(object = Bcell, group.by = "celltype_GSE193_Blood_0.5", pt.size = 0.3, repel=T, label=T, label.size=5);a2

a1|a2
ggsave("UMAP_SingleR_GSE193_Tissue_Blood_res0.5.png", width=10 ,height=6)


#### resolution = 1.3
clusters <- Bcell@meta.data$integrated_snn_res.1.3
cellpred1 <- SingleR(test = testdata, ref = refdata1, labels = refdata1$cell_type1,clusters = clusters 
                     ,assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred1), celltype=cellpred1$labels, stringsAsFactors = F)
# write.csv(celltype,"SingleR_GSE193868_Tissue.csv",row.names = F)
for(i in 1:nrow(celltype)){
  Bcell@meta.data[which(Bcell@meta.data$integrated_snn_res.1.3 == celltype$ClusterID[i]),'celltype_GSE193_Tissue_1.3'] <- celltype$celltype[i]}
a1 <- UMAPPlot(object = Bcell, group.by = "celltype_GSE193_Tissue_1.3", pt.size = 0.3, repel=T, label=T, label.size=5);a1


cellpred2 <- SingleR(test = testdata, ref = refdata2, labels = refdata2$cell_type1,clusters = clusters 
                     ,assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred2), celltype=cellpred2$labels, stringsAsFactors = F)
# write.csv(celltype,"SingleR_GSE193868_Blood.csv",row.names = F)
for(i in 1:nrow(celltype)){
  Bcell@meta.data[which(Bcell@meta.data$integrated_snn_res.1.3 == celltype$ClusterID[i]),'celltype_GSE193_Blood_1.3'] <- celltype$celltype[i]}

a2 <- UMAPPlot(object = Bcell, group.by = "celltype_GSE193_Blood_1.3", pt.size = 0.3, repel=T, label=T, label.size=5);a2

p <- DimPlot(Bcell, reduction = "umap", group.by = "integrated_snn_res.1.3", label = TRUE,
             repel = TRUE);p

a1+NoLegend()|a2+NoLegend()|p+NoLegend()
ggsave("UMAP_SingleR_GSE193_Tissue_Blood_res1.3.png", width=10 ,height=6)

library(scmap)
par(mar=c(13, 4, 2, 0))
barplot(table(cellpred$pruned.labels), las=2)


plot(
  getSankey(
    Bcell$integrated_snn_res.1.3, 
    Bcell$celltype_GSE193_Blood_1.3,
    plot_height = 400
  )
)

plot(
  getSankey(
    Bcell$integrated_snn_res.1.3, 
    Bcell$celltype_GSE193_Tissue_1.3,
    plot_height = 400
  )
)

# when res = 1.3 , MZB-1 = cluster 3, MZB-2 = cluster 8 


##### Kernal Density Estimation ##### 

## MZB
library("Nebulosa")
library("BiocFileCache")

p <- DimPlot(Bcell, label = TRUE, repel = TRUE,cols = cell_type_cols)

p1 <- plot_density(Bcell, c("CD1C", "PLD4","CD27","CR2","IGHM"), joint = TRUE)
#p1 + plot_layout(ncol = 1)
p1[[6]]/p


p2 <- plot_density(Bcell, c("CD1C", "PLD4","CD27"), joint = TRUE)
p2 + plot_layout(ncol = 1)
p2[[4]]/p
ggsave('5-annotationB/MZB_KernalDensity.pdf',width = width.ppi*1.2, height = height.ppi*2)


p3 <- plot_density(Bcell, c("CD1C", "PLD4"), joint = TRUE)
#p1 + plot_layout(ncol = 1)
p3[[3]]/DimPlot(Bcell, label = TRUE, repel = TRUE)

p1[[5]]/p2[[4]]/p3[[3]]


## MZP
p1 <- plot_density(Bcell, c("CD1C", "IGHD","IGHM","PTPRC","ABCB1"), joint = TRUE)
#p1 + plot_layout(ncol = 1)
p1[[6]]+DimPlot(Bcell, label = TRUE, repel = TRUE)
ggsave('5-annotationB/MZP_KernalDensity.png',width = width.ppi*2.2, height = height.ppi*1)
p1[[6]]/DimPlot(Bcell, group.by = "integrated_snn_res.1.3", label = TRUE, repel = TRUE)

## refer to bulk RNA paper
sig.features = c("NR4A1",
                 "NR4A2",
                 "NR4A3",
                 "CD83",
                 "ENTPD1", #CD39
                 "NT5E", #CD73
                 "HLA-G",
                 "TLR7",
                 "TLR10",
                 "CD1A",
                 "CD1C",
                 "CD1D",
                 "IL21R",
                 "TGFB1",
                 "GZMB",
                 "PDCD1", #PD-L1
                 "IL10",
                 "IL10RA","IL10RB" #IL-10R
)
p_list <- plot_density(Bcell, sig.features, joint = TRUE, combine = FALSE)

Bcell <- SetIdent (Bcell, value = "celltype") 
p_list[[length(p_list)]]/DimPlot(Bcell, label = TRUE, repel = TRUE)

p_list[[length(p_list)]]+DimPlot(Bcell, label = TRUE, repel = TRUE)
ggsave('5-annotationB/bulk_reference_NR4A_KernalDensity.pdf',width = width.ppi*2.2, height = height.ppi*1)



#### Manually annotate MZB clusters ####
sel.clust = "integrated_snn_res.1.3"
Bcell <- SetIdent(Bcell, value = sel.clust)
table(Bcell@active.ident) 

##Rearrange order based on cluster tree (hierarchical)
Bcell <- BuildClusterTree(object = Bcell, dims = 1:30)
data.tree <- Tool(Bcell, slot = "BuildClusterTree")
is_tip <- data.tree$edge[,2] <= length(data.tree$tip.label)
ordered_tips <- data.tree$edge[is_tip, 2]
my_levels <- data.tree$tip.label[ordered_tips]
Bcell@active.ident <- factor(x = Bcell@active.ident, levels = my_levels)
table(Bcell@active.ident)

##RNA scaling
Bcell <- NormalizeData(Bcell, assay = "RNA")
all.genes <- rownames(Bcell)
Bcell <- ScaleData(Bcell, features = all.genes, split.by = "sample_type", assay = "RNA")

##visualization
library(viridis)
rna_markers <- c("CD27", "CD38", "IGHM", "IGHD")
#RNA
DefaultAssay(Bcell) <- "RNA"
rna_plot <- FeaturePlot(Bcell, features = rna_markers, min.cutoff = "q15", max.cutoff = "q95", ncol = 1, cols = viridis(200), order = T) & NoAxes() & NoLegend();rna_plot
ggsave(filename = "umap_Bcell_rnaseqmarkers.png", plot = rna_plot, width = 10, height = 40, units = c("cm"))

##RNA features
rna_features_output_split_dot <- c("IGHD", "IGHM", "CD27", "CD1C", "PLD4", "CD5","CR2","FCER2", "CCR7", "MME", "MZB1", "IGLL5", "CD24", "CD38", "COCH", "CXCR5", "ZEB2", "MX1", "BCL6", "BCL7A", "PCNA","MZB.RNA1", "GC.RNA1", "IgG.RNA1", "IgA.RNA1")

plot_alldot <- DotPlot(Bcell, features = rna_features_output_split_dot, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"));plot_alldot
ggsave(filename = "Bcell_dotplot_res0.5.png", plot = plot_alldot, width = 25, height = 15, units = c("cm"))


p <- DimPlot(Bcell, reduction = "umap", group.by = "integrated_snn_res.1.3", label = TRUE,
              repel = TRUE);p


backup = Bcell 
Bcell = backup
Bcell <- RenameIdents(object = Bcell, 
                      `3` = "MZB-1", 
                      `8` = "MZB-2", 
                      `10` = "GCB",
                      `19` = "GCB",
                      `24` = "PB",
                      `22` = "PB",
                      `28` = "PB",
                      `14` = "PB",
                      `4` = "PB",
                      `29` = "Plasmablast",
                      `18` = "PB",
                      `30` = "Plasmablast",
                      `21` = "PB",
                      `27` = "PB",
                      `15` = "DN",
                      `23` = "DN",
                      `5` = "Naive",
                      `1` = "Naive",
                      `11` = "TS",
                      `2` = "TS",
                      `25` = "TS",
                      `16` = "AcB1",
                      `12` = "AcB2",
                      `26` = "AcB3",
                      `9` = "IgM-only",
                      `20` = "Memory",
                      `7` = "Memory",
                      `6` = "Memory",
                      `0` = "Memory",
                      `9` = "Memory",
                      `17` = "Memory",
                      `13` = "Memory")
p <- DimPlot(Bcell, reduction = "umap", label = TRUE,
             repel = TRUE);p
ggsave("UMAP_Bcell.png", width=7 ,height=6)

plot_alldot <- DotPlot(Bcell, features = rna_features_output_split_dot, scale.by = "radius", dot.min = 0, dot.scale = 15, assay = "RNA") +
  RotatedAxis() +
  scale_color_viridis_c() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  guides(size = guide_legend(title = "Percent\nExpressed"), color = guide_colourbar(title = "Average\nExpression"));plot_alldot
ggsave(filename = "Bcell_dotplot.png", plot = plot_alldot, width = 25, height = 15, units = c("cm")) # -> this is Figure 2b



Bcell$celltype <- Bcell@active.ident
saveRDS(Bcell,'Bcell.annotated.rds')


# Markers ----------------------------------------
# selected_cells <- sample(colnames(Bcell), 1000)
# seurat_subset <- subset(Bcell, cells = selected_cells)


Bcell@active.ident <- as.factor(Bcell$celltype)

de_pos_list <- {}
de_neg_list <- {}

for (i in levels(Bcell)){
  condition.diffgenes <- FindMarkers(Bcell, ident.1 = i,
                                     logfc.threshold = 0.25,min.pct = 0.1)
  pos_dat <- subset(condition.diffgenes, avg_log2FC > 0, p_val_adj < 0.05)
  neg_dat <- subset(condition.diffgenes, avg_log2FC < 0, p_val_adj < 0.05)
  pos_marker <- rownames(pos_dat)
  neg_marker <- rownames(neg_dat)
  de_pos_list[[i]] <- pos_marker
  de_neg_list[[i]] <- neg_marker
}

de_pos_summary <- t(plyr::ldply(de_pos_list, rbind, .id = "cluster"))
de_neg_summary <- t(plyr::ldply(de_neg_list, rbind, .id = "cluster"))
out_dir ='5-annotationB/'
file_name <- paste0("pos_degenes_summary.csv")
write.csv(de_pos_summary, file=file.path(out_dir, file_name), col.names = F)
file_name <- paste0("neg_degenes_summary.csv")
write.csv(de_neg_summary, file=file.path(out_dir, file_name), col.names = F)



idents_list <- levels(Bcell$celltype)
for (i in 1:length(idents_list)){
  ident1 <- idents_list[i]
  
  condition.diffgenes <- FindMarkers(Bcell, ident.1 = ident1, only.pos = TRUE)
  
  file_name <- paste0("DEG_", ident1, ".csv")
  write.csv(condition.diffgenes, file=file.path(out_dir, file_name), col.names = F)
}

# HD's MZB DEG list -----
table(Bcell$group)
HD=subset(Bcell,subset=group==c("HD"))
table(HD@active.ident)
HD <- RenameIdents(object = HD, 
                      `AcB3` = "AcB1") 

de_pos_list <- {}
de_neg_list <- {}

for (i in levels(HD)){
  condition.diffgenes <- FindMarkers(HD, ident.1 = i,
                                     logfc.threshold = 0.25,min.pct = 0.1)
  pos_dat <- subset(condition.diffgenes, avg_log2FC > 0, p_val_adj < 0.05)
  neg_dat <- subset(condition.diffgenes, avg_log2FC < 0, p_val_adj < 0.05)
  pos_marker <- rownames(pos_dat)
  neg_marker <- rownames(neg_dat)
  de_pos_list[[i]] <- pos_marker
  de_neg_list[[i]] <- neg_marker
}

de_pos_summary <- t(plyr::ldply(de_pos_list, rbind, .id = "cluster"))
de_neg_summary <- t(plyr::ldply(de_neg_list, rbind, .id = "cluster"))
out_dir ='5-annotationB/'
file_name <- paste0("HD_pos_degenes_summary.csv")
write.csv(de_pos_summary, file=file.path(out_dir, file_name), col.names = F)


# compare with published 2 MZB clusters paper -------

paper.markers <- read.csv('published_tissue_marker_genes.csv',header = T)
View(paper.markers)

de_pos_summary <- read.csv("5-annotationB/pos_degenes_summary.csv",header = T)
colnames(de_pos_summary) <- de_pos_summary[1,]
de_pos_summary <- de_pos_summary[-1,]
rownames(de_pos_summary) <- de_pos_summary[,1]
de_pos_summary <- de_pos_summary[,-1]
View(de_pos_summary)

common_elements1 <- intersect(paper.markers[,'MZB.1'], de_pos_summary[,'MZB-1']);common_elements1
# [1] "LINC01857" "CD1C"      "TNFRSF13B" "AIM2"      "CLECL1"    "CR2"       "SMIM14"   
# [8] "MYC"       "SMARCB1"   "SLA" 
common_elements2 <- intersect(paper.markers[,'MZB.2'], de_pos_summary[,'MZB-2']);common_elements2
# [1] "FABP5" "PSME2" "PRDX1" "ENO1"  "LDHA"  "TPI1"  "RAN"  

common_elements3 <- intersect(paper.markers[,'MZB.1'], de_pos_summary[,'MZB-2']);common_elements3
# [1] "CD1C"      "TNFRSF13B" "PTPN18"    "COTL1"     "CLECL1"    "SMIM14"   
common_elements4 <- intersect(paper.markers[,'MZB.2'], de_pos_summary[,'MZB-1']);common_elements4
# [1] "HSP90AB1"




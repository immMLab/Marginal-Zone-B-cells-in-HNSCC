
GSE234933 <- readRDS("~/Project/5Dataset/GSE234933/GSE234933.rds")
GSE139 <- readRDS("~/Project/GSE139324/GSE139324.rds")
GSE164 <- readRDS("~/Project/GSE164690/GSE164690.rds")

table(GSE181919$patient.id) #patient.id, Gender, Age, tissue.type, subsite, hpv, cell.type
# > table(GSE181919$tissue.type) 
# CA(HNSCC)    LN(LN metastasis)    LP(leukoplakia)    NL(normal) 
# 23088  8204  6527 16420 
table(GSE234933$subclustering) #sample, subclustering, major_state, minor_state, projecTIL_subtype
table(GSE164$sample_type) #group, sample_type, dataset

Idents(GSE164) <- "sample_type"
GSE164 <- RenameIdents(object = GSE164, `CD45p` = "TIL")
GSE164$sample_type = Idents(GSE164)
table(GSE164$sample_type)


# AddMetaData GSE181919 ----------
GSE181919 <- readRDS("~/Project/5Dataset/GSE181919/GSE181919.rds")
GSE181919@meta.data$dataset =  'GSE181919'

# orig.ident
table(GSE181919$orig.ident)
table(GSE181919$patient.id)
GSE181919$orig.ident <- GSE181919$patient.id
GSE181919$orig.ident <- paste0("GSE181_", GSE181919$orig.ident)

patient.id <- c("P4","P6","P7","P8","P9",
                "P15","P16","P17","P21",
                "P22","P26","P30","P31","P33","P38",
                "P43","P46","P51","P57","P59","P60",
                "P84","P86") 
T_stage <- c("T2", "T2", "T2","T1","T2","T1",
             "NA","NA","T4","T2","T2","T2","T2","NA",
             "T4","T3","T2","T2","T1","T3","T4","T2","T2")  
HPV <- c("p16-", "p16-","p16-","p16-","p16-","p16-","p16-","p16-","p16-",
         "P16+","p16-","p16-","p16-","p16-","p16-","P16+","P16+",
         "p16-","P16+","P16+", "p16-","P16+","P16+")  
N_stage <- c("N2a", "N1", "N0","N0","N0","N0",
             "NA","NA","N1","N2","N0","N0","N2","NA",
             "N3b","N0","N1","N1","N0","N2b","N0","N1","N0")  
metadata_df <- data.frame(patient.id, T_stage, HPV, N_stage)
rownames(metadata_df) <- metadata_df$patient.id

GSE181919$T_stage <- metadata_df[GSE181919$patient.id, "T_stage"]
GSE181919$HPV <- metadata_df[GSE181919$patient.id, "HPV"]
GSE181919$N_stage <- metadata_df[GSE181919$patient.id, "N_stage"]

head(GSE181919@meta.data)
table(GSE181919$N_stage)

GSE181919@meta.data <- GSE181919@meta.data[ , !colnames(GSE181919@meta.data) %in% c( "hpv")]

## subsite
Idents(GSE181919)="subsite"
GSE181919 <- RenameIdents(GSE181919, 
                          `OC` = "Oral cavity",`HP` = "Larynx/Hypopharynx", `OP` = "Oropharynx")
GSE181919$subsite <- GSE181919@active.ident
table(GSE181919$subsite)


## Gender  
table(GSE181919$Gender)
Idents(GSE181919)="Gender"
GSE181919 <- RenameIdents(GSE181919, 
                          `F` = "Female",`M` = "Male")
GSE181919$Gender <- GSE181919@active.ident
table(GSE181919$Gender)

## group
table(GSE181919$tissue.type)
Idents(GSE181919)="tissue.type"
GSE181919 <- RenameIdents(GSE181919, 
                          `CA` = "HNSCC",`LN` = "LN metastasis", `LP` = "Leukoplakia",`NL` = 'Normal')
GSE181919$group <- GSE181919@active.ident
table(GSE181919$group)

## sample_type
Idents(GSE181919)="tissue.type"
GSE181919 <- RenameIdents(GSE181919, 
                          `CA` = "TIL",`LN` = "LN TIL", `LP` = "Leukoplakia",`NL` = 'Normal')
GSE181919$sample_type <- GSE181919@active.ident
table(GSE181919$sample_type)

head(GSE181919@meta.data)
saveRDS(GSE181919,"~/Project/5Dataset/GSE181919/GSE181919.rds")


# AddMetaData GSE234933 ----------

# dataset

GSE234933@meta.data$dataset =  'GSE234933'

table(GSE234933$orig.ident)
GSE234933$patient.id <- GSE234933$orig.ident
GSE234933$orig.ident <- paste0("GSE234_", GSE234933$orig.ident)

tail(GSE234933@meta.data)
metadata_df <- read.table("~/Project/5Dataset/GSE234933/MGH_HNSCC_sample_annotation.txt", sep = "\t", header = TRUE, row.names = 1,stringsAsFactors = FALSE)

GSE234933$Gender <- metadata_df[GSE234933$patient.id, "Sex"]
GSE234933$Age <- metadata_df[GSE234933$patient.id, "Age"]
GSE234933$HPV <- metadata_df[GSE234933$patient.id, "HPV.Status"]
GSE234933$Smoking.history <- metadata_df[GSE234933$patient.id, "Smoking.history"]
GSE234933$subsite <- metadata_df[GSE234933$patient.id, "Original.anatomic.site"]
GSE234933$tissue.type <- metadata_df[GSE234933$patient.id, "Anatomic.location.of.scRNA.seq.specimen"]

head(GSE234933@meta.data)
table(GSE234933$tissue.type)
table(GSE234933$subsite)


## group
table(GSE234933$tissue.type)
Idents(GSE234933)="tissue.type"
GSE234933 <- RenameIdents(GSE234933, 
                          `Primary` = "HNSCC",`Unknown primary of the head and neck` = "HNSCC", 
                          `Distant metastasis (Liver)` = "Metastasis",
                          `Distant metastasis (Lung)` = "Metastasis",
                          `Distant metastasis (pleura)` = "Metastasis",
                          `Distant metastasis (skin)` = "Metastasis",
                          `Distant metastasis (sternum)` = "Metastasis",
                          `Locoregional recurrence` = "Local recurrence")
GSE234933$group <- GSE234933@active.ident
table(GSE234933$group)

## sample_type
Idents(GSE234933)="tissue.type"
GSE234933 <- RenameIdents(GSE234933, 
                          `Primary` = "TIL",`Unknown primary of the head and neck` = "TIL", 
                          `Distant metastasis (Liver)` = "Metastasis",
                          `Distant metastasis (Lung)` = "Metastasis",
                          `Distant metastasis (pleura)` = "Metastasis",
                          `Distant metastasis (skin)` = "Metastasis",
                          `Distant metastasis (sternum)` = "Metastasis",
                          `Locoregional recurrence` = "Local recurrence")
GSE234933$sample_type <- GSE234933@active.ident
table(GSE234933$sample_type)


## HPV
table(GSE234933$HPV)
Idents(GSE234933)="HPV"
GSE234933 <- RenameIdents(GSE234933, `Negative` = "p16-",`Positive` = "p16+")
GSE234933$HPV <- GSE234933@active.ident
head(GSE234933@meta.data)

saveRDS(GSE234933,"~/Project/5Dataset/GSE234933/GSE234933.rds")

# AddMetaData GSE139324 ----------
GSE139 <- readRDS("~/Project/GSE139324/GSE139324.rds")

## patient.id
table(GSE139$orig.ident)
GSE139$patient.id <- gsub('GSE139*_','',GSE139$orig.ident)


Idents(GSE139)="patient.id"
GSE139 <- RenameIdents(GSE139, 
                          `HD_1_Tonsil` = "HD_Tonsil1_Lymphocytes",
                          `HD_2_Tonsil` = "HD_Tonsil2_Lymphocytes", 
                          `HD_3_Tonsil` = "HD_Tonsil3_Lymphocytes",
                          `HD_4_Tonsil` = "HD_Tonsil4_Lymphocytes",
                          `HD_5_Tonsil` = "HD_Tonsil5_Lymphocytes")
GSE139$patient.id <- GSE139@active.ident
GSE139$patient.id <- sub("(.*)_(.*)$", "\\1", GSE139$patient.id)
table(GSE139$patient.id)

metadata_df <- read.csv("~/Project/GSE139324/GSE139324_ClinicalSheet.csv", header = TRUE, row.names = 1,stringsAsFactors = FALSE)
head(metadata_df)

GSE139$Gender <- metadata_df[GSE139$patient.id, "Sex"]
GSE139$Age <- metadata_df[GSE139$patient.id, "Age"]
GSE139$Race <- metadata_df[GSE139$patient.id, "Race"]
GSE139$HPV <- metadata_df[GSE139$patient.id, "Tumor.p16.status"]
GSE139$Alcohol.history <- metadata_df[GSE139$patient.id, "Alcohol.Use"]
GSE139$Smoking.history <- metadata_df[GSE139$patient.id, "Tobacco.Use"]
GSE139$subsite <- metadata_df[GSE139$patient.id, "Tumor.site"]
GSE139$T_stage <- metadata_df[GSE139$patient.id, "Pathologic.Stage"]
GSE139$N_stage <- metadata_df[GSE139$patient.id, "Pathological.Node"]

head(GSE139@meta.data)
table(GSE139$N_stage)

## N_stage
table(GSE139$N_stage)
Idents(GSE139)="N_stage"
GSE139 <- RenameIdents(GSE139, 
                          `N2B` = "N2b",`N2C` = "N2c", 
                          `N3B` = "N3b",`NX` = "Nx")
GSE139$N_stage <- GSE139@active.ident
table(GSE139$N_stage)

## subsite
table(GSE139$subsite)
Idents(GSE139)="subsite"
GSE139 <- RenameIdents(GSE139, 
                       `Base of tongue` = "Oral cavity",`Buccal mucosa` = "Oral cavity", 
                       `Floor of mouth` = "Oral cavity",`Larynx` = "Larynx/Hypopharynx",
                       `Lower gum` = "Oral cavity",`Mandable` = "Oral cavity", 
                       `Supraglottis` = "Larynx/Hypopharynx",
                       `Tongue` = "Oral cavity",`Tonsil` = "Oropharynx")
GSE139$subsite <- GSE139@active.ident
table(GSE139$subsite)

## group
table(GSE139$group)
table(GSE139$HPV)

Idents(GSE139)="orig.ident"
saveRDS(GSE139,"~/Project/GSE139324/GSE139324.rds")




# AddMetaData GSE164690 ----------
GSE164 <- readRDS("~/Project/GSE164690/GSE164690.rds")


## patient.id
table(GSE164$orig.ident)
GSE164$patient.id <- sub("^.*?_(.*?)_.*$", "\\1", GSE164$orig.ident)
table(GSE164$patient.id)

patient.id <- c("HN01","HN02","HN03","HN04","HN05",
                "HN06","HN07","HN08","HN09",
                "HN10","HN11","HN12","HN13","HN14","HN15",
                "HN16","HN17","HN18") 
Gender <- c("Male", "Female", "Male","Male","Female","Male",
             "Male","Female","Female","Male","Male","Male","Male","Male",
             "Female","Male","Male","Male")  
Age <- c("Male", "Female", "Male","Male","Female","Male",
            "Male","Female","Female","Male","Male","Male","Male","Male",
            "Female","Male","Male","Male")  
Smoking.history <- c("Yes", "No", "No","Yes","Yes","Yes",
            "Yes","Yes","Yes","No","No","Yes","No","Yes",
            "Yes","Yes","Yes","Yes") 
Alcohol.history <- c("No", "No", "No","Yes","Yes","Yes",
            "Yes","Yes","Yes","Yes","No","Yes","No","Yes",
            NA,"Yes","Yes","Yes")  
T_stage <- c("T4A", "T3", "T4A","T3","T3","T3",
             "T3","T1","T3","T3","T2","T2","T2","T1",
             "T2","T2","T1","T2")  
N_stage <- c("N2b", "N2a", "N0","N1","N3b","N0",
             "N0","N0","N2b","N0","N0","N1","N0","N1",
             "N0","N1","N1","N2")  
M_stage <- c("M0", "M0", "M0","M0","M0","M0",
             "M0","M0","M0","M0","M0","M0","M0","M0",
             "M0","M0","M0","M0")  
HPV <- c("p16-", "p16-","p16-","p16-","p16-","p16-","p16-","p16-","p16-",
         "P16-","p16-","p16+","p16+","p16+","p16-","P16+","P16+",
         "p16+") 
subsite <-  c("Oral cavity", "Oral cavity", "Oral cavity","Oral cavity","Oral cavity","Oral cavity",
              "Larynx/Hypopharynx","Oral cavity","Oral cavity","Oral cavity","Oral cavity","Oropharynx","Oropharynx","Oropharynx",
              "Oral cavity","Oropharynx","Oropharynx","Oropharynx")  
metadata_df <- data.frame(patient.id, Gender, Age, Smoking.history, Alcohol.history,
                          T_stage, N_stage, M_stage, HPV, subsite)
rownames(metadata_df) <- metadata_df$patient.id
metadata_df$patient.id <- NULL

GSE164$Gender <- metadata_df[GSE164$patient.id, "Gender"]
GSE164$Age <- metadata_df[GSE164$patient.id, "Age"]
GSE164$HPV <- metadata_df[GSE164$patient.id, "HPV"]
GSE164$Smoking.history <- metadata_df[GSE164$patient.id, "Smoking.history"]
GSE164$Alcohol.history <- metadata_df[GSE164$patient.id, "Alcohol.history"]
GSE164$subsite <- metadata_df[GSE164$patient.id, "subsite"]
GSE164$T_stage <- metadata_df[GSE164$patient.id, "T_stage"]
GSE164$N_stage <- metadata_df[GSE164$patient.id, "N_stage"]
GSE164$M_stage <- metadata_df[GSE164$patient.id, "M_stage"]

head(GSE164@meta.data)
table(GSE164$N_stage)

## N_stage
table(GSE164$N_stage)
Idents(GSE164)="N_stage"
GSE164 <- RenameIdents(GSE164, 
                       `N2B` = "N2b",`N2C` = "N2c", 
                       `N3B` = "N3b",`NX` = "Nx")
GSE164$N_stage <- GSE164@active.ident
table(GSE164$Alcohol.history)

## group
table(GSE164$group)

saveRDS(GSE164,"~/Project/GSE164690/GSE164690.rds")


# Merge ----------
setwd('~/Project/5Dataset')

options(warn=-1) # turn off warning message globally
rm(list=ls());gc()
options(stringsAsFactors = F) 
options(as.is = T)
source('scRNA_scripts/lib.R')


getwd()
set.seed(123)  

GSE181919 <- readRDS("~/Project/5Dataset/GSE181919/GSE181919.rds")
GSE234933 <- readRDS("~/Project/5Dataset/GSE234933/GSE234933.rds")
GSE139 <- readRDS("~/Project/GSE139324/GSE139324.rds")
GSE164 <- readRDS("~/Project/GSE164690/GSE164690.rds")

GSE181 <- subset(GSE181919,subset=cell.type==c("B_Plasma.cells","T.cells","Macrophages",
                                               "Dendritic.cells","Macrophages","Mast.cells"))
GSE234 <- subset(GSE234933,subset=subclustering==c("Lymphocytes","Myeloid cells"))

sce_list <- list(GSE181, GSE234, GSE139, GSE164)
sce.all=merge(x=sce_list[[1]],
              y=sce_list[ -1 ])  

saveRDS(sce.all,"sce.all.rds")

# QC ----------

dir.create("./1-QC")
setwd("./1-QC")

source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
table(sce.all.filt$orig.ident)
saveRDS(sce.all.filt,"sce.all.filt.rds")
setwd('../')


# Harmony ----------

rm(list=ls())
options(stringsAsFactors = F) 
source('scRNA_scripts/lib.R')
sce.all.filt <- readRDS('1-QC/sce.all.filt.rds')

dir.create("2-harmony")
setwd("2-harmony")
source('../scRNA_scripts/harmony.R')

sce.all.int = run_harmony(sce.all.filt)
a <- UMAPPlot(object = sce.all.int, pt.size = 0.3, repel=T, label=T, label.size=5,raster=FALSE);a
setwd('../')

table(Idents(sce.all.int))
table(sce.all.int$seurat_clusters)
table(sce.all.int$RNA_snn_res.0.1) 
table(sce.all.int$RNA_snn_res.0.8) 

getwd()
dir.create('check-by-0.1')
setwd('check-by-0.1')
sel.clust = "RNA_snn_res.0.1"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
sp='human'
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()



dir.create('check-by-0.8')
setwd('check-by-0.8')
sel.clust = "RNA_snn_res.0.8"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
sp='human'
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

last_markers_to_check


# SingleR annotation ----------

library(SingleR)
library(tidyverse)

dir.create("3-CellType")
scRNA = sce.all.int

length(levels(Idents(scRNA)))
table(scRNA@meta.data$seurat_clusters)
#table(scRNA@active.ident) 

metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'3-CellType/cell_cluster.csv',row.names = F)

refdata <- MonacoImmuneData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters

cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellpred$pruned.labels
a <- plotScoreHeatmap(cellpred, clusters=cellpred@rownames, fontsize.row = 9,show_colnames = T)
ggsave("3-CellType/Monaco_heatmap.png", a, width=7 ,height=3)
ggsave("3-CellType/Monaco_heatmap.pdf", a, width=7 ,height=3)


celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"3-CellType/celltype_Monaco.csv",row.names = F)
scRNA@meta.data$celltype_Monaco = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_Monaco'] <- celltype$celltype[i]}

a <- UMAPPlot(object = scRNA, group.by = "celltype_Monaco", pt.size = 0.3, repel=T, label=T, label.size=5)
ggsave("3-CellType/UMAP_celltype_Monaco.png", a, width=7 ,height=6)
ggsave("3-CellType/UMAP_celltype_Monaco.pdf", a, width=7 ,height=6)

ref_DICE <- DatabaseImmuneCellExpressionData()
# saveRDS(ref_DICE,'../SingleR_ref/ref_DICE_1561s.RData')
# saveRDS(ref_DICE,'../SingleR_ref/ref_DICE_1561s.rds')
# load('../SingleR_ref/ref_DICE_1561s.RData')
refdata <- ref_DICE
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                    clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")

cellpred$pruned.labels
b <- plotScoreHeatmap(cellpred, clusters=cellpred@rownames, fontsize.row = 9,show_colnames = T)
ggsave("3-CellType/DICE_heatmap.png", b, width=7 ,height=3)
ggsave("3-CellType/DICE_heatmap.pdf", b, width=7 ,height=3)

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"3-CellType/celltype_DICE.csv",row.names = F)
scRNA@meta.data$celltype_DICE = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_DICE'] <- celltype$celltype[i]}
b <- UMAPPlot(object = scRNA, group.by = "celltype_DICE", pt.size = 0.3, repel=T, label=T, label.size=5)
ggsave("3-CellType/UMAP_celltype_DICE.png", b, width=7 ,height=6)
ggsave("3-CellType/UMAP_celltype_DICE.pdf", b, width=7 ,height=6)

p = a+b
ggsave("3-CellType/Monaco_DICE.png", p, width=12 ,height=5)
saveRDS(scRNA,'3-CellType/sce.all.singleR.rds')

sce.all.int = readRDS('2-harmony/sce.all_int.rds')
sp='human'

DotPlot(sce.all.int , features = c('DCT','PMEL','CD207',
                                   'KRT1','KRT10','KRT14'),
        group.by = 'RNA_snn_res.0.1' )  + 
  coord_flip() + 
  theme(axis.text.x=element_text(angle=45,hjust = 1))

 if(F){
  sce.all.int
  celltype=data.frame(ClusterID=0:13 ,
                      celltype= 0:13)     
  celltype[celltype$ClusterID %in% c( 1,12  ),2]='lym'   
  celltype[celltype$ClusterID %in% c( 11,2 ),2]='Myeloids'  
  celltype[celltype$ClusterID %in% c( 10 ),2]='endo'  
  celltype[celltype$ClusterID %in% c( 7 ),2]='fibro'  
  celltype[celltype$ClusterID %in% c( 9),2]='SMC'
  celltype[celltype$ClusterID %in% c( 8 ),2]='cycle' 
  celltype[celltype$ClusterID %in% c( 13 ),2]='epi'
  
  celltype[celltype$ClusterID %in% c( 5 ),2]='Lang'
  celltype[celltype$ClusterID %in% c( 4 ),2]='Melan'
  celltype[celltype$ClusterID %in% c( 0,3,6 ),2]='Kera'
  
  
  head(celltype)
  celltype
  table(celltype$celltype)
  sce.all.int@meta.data$celltype = "NA"
  
  for(i in 1:nrow(celltype)){
    sce.all.int@meta.data[which(sce.all.int@meta.data$RNA_snn_res.0.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
  Idents(sce.all.int)=sce.all.int$celltype
  
  dir.create('check-by-celltype')
  setwd('check-by-celltype')
  sel.clust = "celltype"
  sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
  table(sce.all.int@active.ident) 
  source('../scRNA_scripts/check-all-markers.R')
  setwd('../') 
  getwd()
}
# phe=sce.all.int@meta.data
# save(phe,file = 'phe.Rdata')


## check integration ## 
sce.all.int = sce
a1 <- DimPlot(object = sce, group.by = "orig.ident", pt.size = 0.3, repel=T, label=F, label.size=5,raster=T)+NoLegend();a1
ggsave("3-CellType/integration_orig.ident.pdf", width = width.ppi*1,height = height.ppi)
a2 <- DimPlot(object = sce, group.by = "dataset", pt.size = 0.3, repel=T, label=F, label.size=5,raster=T)
ggsave("3-CellType/integration_dataset.pdf", width = width.ppi*1,height = height.ppi)
a3 <- DimPlot(object = sce, group.by = "group", pt.size = 0.3, repel=T, label=F, label.size=5,raster=T)
ggsave("3-CellType/integration_group.pdf", width = width.ppi*1,height = height.ppi)
a1|a2|a3
ggsave("3-CellType/integration.pdf", width = width.ppi*3,height = height.ppi)


## check QC metrics ##

# Create histograms of all metrics
p <- (VlnPlot(sce, 'nCount_RNA', group.by = 'dataset', pt.size = 0) + NoLegend() |
        VlnPlot(sce, 'nFeature_RNA', group.by = 'dataset', pt.size = 0) + NoLegend()) /
  (VlnPlot(sce, 'percent_mito', group.by = 'dataset', pt.size = 0) + NoLegend() |
     VlnPlot(sce, 'percent_ribo', group.by = 'dataset', pt.size = 0) + NoLegend())
ggsave('3-CellType/QC_metrics_by_dataset.png',p)

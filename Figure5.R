library(tidyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggpubr)
library(scran)
library(ggforce)
library(gghalves)
library(ggridges)
library(scDblFinder)
library(dplyr)

##set the width and height of plots
width.ppi=6.5
height.ppi=5


# remotes::install_version("matrixStats", version="1.1.0") # restart your session and run previous scripts
# UseName error

#Import data
Bcell <- readRDS("~/Project/5Dataset/6StratifiedAnalysis/Bcell.rds")
table(Bcell$sample_type)
sce=subset(Bcell,subset=sample_type==c("TIL")) #只分析primary HNSCC的tissue samples - i.e. TIL

library(EnhancedVolcano)
mz.de.markers <- FindMarkers(sce, ident.1 = "MZB-2", ident.2 = "MZB-1",  logfc.threshold = 0.25,min.pct = 0.1)

markerslabel <- c("DDX21", "TXNIP", "MIF", "NCL", "NME1", "FOS", "HLA-DRB1", "YBX1", "HLA-DQA1", "PSME2", "HLA-DRA", "PTMA", "FABP5", "HSP90AB1",
                  "TXNIP", "MALAT1", "LINC00926", "LTB", "CD37", "CD27", "AIM2", "CD83",'RPS15A','RPS4X','RPS18','RPS14','EEF1A1','HLA-DRB1','FCER2','MS4A1','LTB','HLA-DRA','HSP90AB1','HLA-DQA1',
                  "ATP8B4","SLC12A8","GNG11")

mz_volcano <- EnhancedVolcano(mz.de.markers, lab = rownames(mz.de.markers), 
                              x = "avg_log2FC", y = "p_val",
                              pCutoff = 10e-6, FCcutoff = 1, 
                              title = 'MZB-1 versus MZB-2', subtitle = NULL,
                              labSize = 4, 
                              # selectLab = markerslabel, 
                              labCol = 'black', labFace = 'bold',
                              boxedLabels = T, drawConnectors = T);mz_volcano


suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
})
gene_list <- c(
  c("FCRL4", "FCRL5", "ITGAX", "TBX21", "CR2"),
  c("CD80", "CD86", "FAS"),
  c("BHLHE40", "RUNX2", "RUNX3", "SPIB", "ZBTB32", "ZEB2"),
  c("IFNGR1", "IL2RB", "IL12RB1", "TNFRSF13B", "TNFRSF17"),
  c("CXCR3", "CCR1", "CCR5", "CCR6"),
  c("HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DQA2", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DOA", "HLA-DOB")
)

row_split <- c(rep("1", 5),
               rep("2", 3),
               rep("3", 6),
               rep("4", 5),
               rep("5", 4),
               rep("6", 12))

if (!"CellID" %in% colnames(sce@meta.data)) {
  sce$CellID <- colnames(sce)
}

meta <- sce@meta.data %>%
  dplyr::select(celltype, CellID) %>%
  dplyr::filter(celltype %in% c("MZB-1", "MZB-2")) %>%

  mutate(anno_order = ifelse(celltype == "MZB-1", 1, 2)) %>%
  arrange(anno_order) %>%
  mutate(col_split = ifelse(celltype == "MZB-1", 1, 2))


cells_use <- meta$CellID


DefaultAssay(sce) <- "RNA"
expr <- sce@assays$RNA@data

genes_present <- intersect(gene_list, rownames(expr))

mat <- as.matrix(expr[genes_present, cells_use, drop = FALSE])


mat_scaled <- t(scale(t(mat)))  # 行标准化


df_rows <- data.frame(
  gene = gene_list,
  row_split = row_split,
  stringsAsFactors = FALSE
)
df_rows <- df_rows[df_rows$gene %in% genes_present, , drop = FALSE]

df_rows <- df_rows[match(rownames(mat_scaled), df_rows$gene), , drop = FALSE]


ha <- HeatmapAnnotation(
  Group = factor(meta$col_split, levels = c(1, 2), labels = c("MZB-1", "MZB-2")),
  col = list(Group = c("MZB-1" = "#009FB9", "MZB-2" = "#E60012"))
)


col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#0F7B9F", "white", "#D83215"))

ht <- ComplexHeatmap::Heatmap(
  mat_scaled,
  name            = "Expression",
  col             = col_fun,
  top_annotation  = ha,
  row_split       = df_rows$row_split,
  column_split    = meta$col_split,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp    = grid::gpar(fontsize = 7),
  width  = ncol(mat_scaled) * unit(0.003, "inch"),
  height = nrow(mat_scaled) * unit(0.15, "inch")
)

draw(ht)





df <- read.csv("/home/data/t210344/Project/5Dataset/11-GSVA/enrichGO_all.csv",
               check.names = FALSE)

if ("...1" %in% names(df)) df <- df[ , setdiff(names(df), "...1")]

topN <- 30
df_show <- head(df, topN)

df_show$Description <- factor(df_show$Description,
                              levels = rev(df_show$Description))  

ggplot2::ggplot(df_show,
                ggplot2::aes(x = Description, y = -log10(p.adjust))) +
  ggplot2::geom_bar(stat = "identity", width = 0.8, fill = "salmon1") +
  ggplot2::coord_flip() +
  ggplot2::labs(x = "GO term", y = expression(-log[10](adj~p)),
                title = "GO BP (use original CSV order)") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_text(size = 6, color = "black"),
    axis.text.y  = ggplot2::element_text(size = 6, color = "black"),
    plot.title   = ggplot2::element_text(hjust = 0.5, size = 8),
    panel.grid   = ggplot2::element_blank()
  )


# 加载包
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})


load("/home/data/t210344/Project/5Dataset/11-GSVA/GSEA_result.RData")  

library(enrichplot)
library(ggplot2)
go_id <- "GO:0050776"

go_title <- GO_kk@result[go_id, "Description"]


gseap1 <- gseaplot2(
  GO_kk,
  geneSetID = go_id,
  title = go_title,
  color = "firebrick",          
  base_size = 14,              
  rel_heights = c(1.5, 0.5, 1), 
  subplots = 1:3,              
  ES_geom = "line",            
  pvalue_table = TRUE          
)

gseap1 <- gseap1 +
  ggplot2::theme_minimal(base_size = 14, base_family = "Arial") +
  ggplot2::theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  ggplot2::labs(
    title = paste0("GSEA: ", go_title),
    x = "Rank in ordered gene list",
    y = "Enrichment Score (ES)"
  )

print(gseap1)















library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tibble)

# 差异基因
markers <- FindMarkers(
  sce, ident.1 = "MZB-2", ident.2 = "MZB-1", 
  group.by = "celltype", 
  only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox"
) %>% rownames_to_column("gene")

markers <- subset(markers, p_val_adj < 0.05 & abs(avg_log2FC) > 0.15)

# GO BP 富集
ego_BP <- enrichGO(
  gene          = markers$gene,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 0.05,   # 建议别设太宽松
  qvalueCutoff  = 0.05
)

# 转成数据框
df_BP <- as.data.frame(ego_BP)
df_BP <- df_BP[order(df_BP$p.adjust), ]  # 按校正 p 排序
df_BP_top30 <- head(df_BP, 30)

# 按 Description 顺序画图
df_BP_top30$Description <- factor(df_BP_top30$Description,
                                  levels = rev(df_BP_top30$Description))

# 画 barplot
ggplot(df_BP_top30, aes(x = Description, y = -log10(p.adjust))) +
  geom_bar(stat="identity", width=0.8, fill='salmon1') +
  coord_flip() +
  xlab("GO term") + ylab("-log10(p.adjust)") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(family = "ArialMT", size = 6, color = "black"),
    axis.text.y  = element_text(family = "ArialMT", size = 6, color = "black"),
    axis.line.x  = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    plot.margin  = unit(c(1, 1, 1, 1), "char"),
    text         = element_text(family = "ArialMT", size = 6),
    plot.title   = element_text(hjust = 0.5, size = 8),
    panel.grid   = element_blank(),
    axis.line    = element_line(linetype = 1, color = "black", size = 0.3),
    axis.ticks   = element_line(linetype = 1, color = "black", size = 0.3)
  ) +
  ggtitle("MZB-2 enriched pathways")


load("/home/data/t210344/Project/5Dataset/11-GSVA/GSEA_result.RData")  

library(enrichplot)

p2 <- gseaplot(GO_kk,
               geneSetID = "GO:0050776",
               title = GO_kk@result["GO:0050776", "Description"],
               by = "runningScore",  
               color = "red")        
print(p2)

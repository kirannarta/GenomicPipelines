
## This script to be run after getting output from Immunedeconv_barlot.R script

## Do: conda activate deconvolution

args <- commandArgs(trailingOnly = TRUE)

library(immunedeconv)
library(tidyverse)
library(data.table)
library(scales)
library(broom)
library(reshape2)
library(ComplexHeatmap)


## rsem_filesPath <- args[1]

## Heatmaps
#data_all <- read.table("BarPlots_includingMuscle/Immune_deconv_allTools_bigTable.txt", header=T)
data_all <- read.table(args[1], header=T, sep="\t") 
cell_types <- data.frame(Var1 = data_all$Var1, CellTypes=chartr(" ", "_", (data_all$Var1))) %>% unique
#cell_annot <- read.table("cell_annot.txt", header=T)
cell_annot <- read.table("cell_types_annot.txt", header=T, sep="\t", row.names=1)
cell_types_annot <- merge(cell_types, cell_annot, by.x="CellTypes", by.y="newCellType", all.x=T)
write.table(cell_types_annot, "cell_types_annot.txt", quote=F, row.names=F, sep="\t")

### Edit cell annotations and create broad categories
##=============================================================== 
## start here #####

data_all <- read.table(args[1], header=T, sep="\t") 
cell_annot <- read.table("cell_types_annot.txt", header=T, row.names=1, sep="\t")
big_mat <- merge(data_all, cell_annot, by="Var1", all.x=T)

big_mat3 <- big_mat %>% .[complete.cases(.), ] %>% group_by(L1) %>% mutate(rescaleMinMaxTool = rescale(as.numeric(as.character(value)), to=c(0,1))) %>% data.frame() %>% group_by(L1, CellTypes) %>% mutate(rescaleMinMaxCellTypeTool = rescale(as.numeric(as.character(value)), to=c(0,1))) %>% data.frame()

big_mat3$CT_Tool <- paste(big_mat3$CellTypes, gsub("res_", "", big_mat3$L1), sep=":")

### Prepare matrix for heatmap
Mat_all <- big_mat3[, c("Var2", "rescaleMinMaxCellTypeTool", "CT_Tool")] %>% acast(., CT_Tool~Var2, value.var= "rescaleMinMaxCellTypeTool")

### Annotation rows and columns
ann_row <- big_mat3[, c("CT_Tool", "BroadCategory")] %>% unique %>% .[order(.$CT_Tool),]
ann_row1 <- rowAnnotation(BroadCategory = ann_row$BroadCategory)

ann_col <- big_mat3[, c("Var2", "condition")] %>% unique() %>% .[order(.$Var2),]
ann_col1 <- HeatmapAnnotation(df= ann_col, which='col')

split_by_col <- factor(big_mat3[, c("Var2", "condition")] %>% unique() %>% .[order(.$Var2), "condition"])
split_by_row <- factor(ann_row[order(ann_row$CT_Tool), "BroadCategory"])

pdf(paste("Heatmap_split.pdf"), height=30, width=15)
h <- Heatmap(Mat_all, name = "MinMaxScale", column_title = "ModelsCompared", row_title = "CellTypes", bottom_annotation= ann_col1, right_annotation = ann_row1, column_split = split_by_col, row_split = split_by_row, col = colorRamp2(c(0,1), c("white", "red")))
print(h)
dev.off()


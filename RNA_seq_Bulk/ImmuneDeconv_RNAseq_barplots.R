## Rscript for immune deconv 
## Do: conda activate deconvolution

## usage: Rscript ImmuneDeconv_RNAseq_barplots.R ../RSEM_files/ ../ImmuneDeconv/ ../metadata_Radiation_Arvind.txt

args <- commandArgs(trailingOnly = TRUE)

library(immunedeconv)
library(tidyverse)
library(data.table)
library("ggpubr")
#--------------------------------------------------------------------------
## Load input data: RSEM files
rsem_filesPath <- args[1]
setwd(rsem_filesPath)  #to where RSEM files are

DF = do.call(cbind,
         lapply( list.files(pattern=".*genes.results"),
                     FUN=function(x) {
            aColumn = read.table(x,header=T)[,c("gene_id", "TPM")];
            colnames(aColumn)[2] = x;
            aColumn;
             }
            )
        )

DF = DF[,!duplicated(colnames(DF))]
colnames(DF) <- gsub("_rsem.genes.results", "", colnames(DF))

deconv_folderPath <- args[2]
setwd(deconv_folderPath)

#--------------------------------------------------------------------------
## Prepare input file: #### Select rows with unique genesymbols -> get rowsums > 0 and if duplicate genes, select ones with max rowsums
DF$GeneSymbol <- map_chr(str_split(DF$gene_id, "_"),2)
DF <- DF %>% mutate(rowsums = rowSums(.[,sapply(., is.numeric)])) %>% .[.$rowsums> 0, ]

DF1 <-  DF %>% group_by(GeneSymbol) %>% slice(which.max(rowsums))  %>% data.frame
rownames(DF1) <- DF1$GeneSymbol
DF1<- subset(DF1, select = -c(GeneSymbol, rowsums, gene_id))
write.table(DF1, "RSEM_table_immunedeconv_inputFiles.txt", sep="\t", quote=F)

#--------------------------------------------------------------------------

## ImmuneDeconv: Mouse based methods: mMCPcounter, seqimmucc, dcq, base
res_mMCPcounter_m38 <- deconvolute_mmcp_counter(DF1, genome="GCRm38")
res_seqimmucc <- deconvolute_seqimmucc(DF1, "LLSR")
res_dcq <- deconvolute_dcq(DF1, n_repeats = 10,combine_cells = TRUE)
res_base <- deconvolute_base_algorithm(DF1, n_permutations = 100,combine_cells = TRUE)


#--------------------------------------------------------------------------
## ImmuneDeconv: Human based methods: quantiseq, XCell, timer, abis, extimate, Cibersort; Convert to human geneids

	# Convert mouse_gene_ids to human_gene_ds
	#tpm_ups2_human <- immunedeconv::mouse_genes_to_human(tpm_ups2)

	mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep = "\t")

	find_corr_gene <- function(gene, mouse_human_genes_df, convert_to = c("human", "mouse")) {
          if (convert_to == "human") {
            orgn.name <- "mouse, laboratory"
            new.orgn <- "human"
          } else {
            orgn.name <- "human"
            new.orgn <- "mouse, laboratory"
          }

          class_key <- (mouse_human_genes_df %>%
            filter(Symbol == gene & Common.Organism.Name == orgn.name))[["DB.Class.Key"]]

          if (!identical(class_key, integer(0))) {
            output <- NULL
            new_genes <- (mouse_human_genes_df %>% filter(DB.Class.Key == class_key & Common.Organism.Name == new.orgn))[, "Symbol"]

            for (new_gene in new_genes) {
              output <- append(output, new_gene)
            }

            if (!is.null(output)) {
              return(data.frame(
                "new_gene" = output,
                "old_gene" = gene
              ))
            }
          }
        }



gene.names <- rownames(DF1)
DF1$gene_name <- gene.names
genes.retrieved <- map_dfr(gene.names, function(x) find_corr_gene(x, mouse_human_genes, "human"))
newGenes.counts <- DF1 %>% left_join(., genes.retrieved, by = c("gene_name" = "old_gene")) %>% select(., -c("gene_name")) %>%select(., c("new_gene", everything()))
colnames(newGenes.counts)[1] <- "gene_name"
newGenes.counts <- newGenes.counts[!(is.na(newGenes.counts$gene_name)), ] %>%group_by(gene_name) %>%summarise_all(median)
newGenes.counts <- as.data.frame(newGenes.counts)
rownames(newGenes.counts) <- newGenes.counts$gene_name
newGenes.counts <- select(newGenes.counts, -c("gene_name"))

DF1_human <- newGenes.counts	
write.table(DF1_human, "RSEM_table_immunedeconv_inputFiles_human.txt", sep="\t", quote=F)


res_xCell <- deconvolute_xcell(DF1_human, arrays = FALSE)
res_timer <- deconvolute_timer(DF1_human, indications = rep("sarc", ncol(DF1_human)))
res_abis <- deconvolute_abis(DF1_human, arrays = FALSE)
res_estimate <- deconvolute_estimate(DF1_human)
res_quantiseq <- deconvolute_quantiseq(DF1_human, tumor = FALSE, arrays = FALSE, scale_mrna=TRUE)
res_epic <- deconvolute_epic(DF1_human, tumor=TRUE, scale_mrna=FALSE)
res_consensus_tme <- deconvolute_consensus_tme(DF1_human, indications = rep("sarc", ncol(DF1_human))  , method= 'ssgsea')

	######## Cibersort
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("/work/morrissy_lab/Kiran/RESOURCES/Software/Cibersort/LM22.txt")
res_cibersort <- deconvolute(DF1_human, "cibersort") 
res_cibersort <- data.frame(res_cibersort)
rownames(res_cibersort) <- res_cibersort$cell_type
res_cibersort <- subset(res_cibersort, select= -(cell_type))
#rownames(res_cibersort) <- chartr(" ", "_", rownames(res_cibersort))
res_cibersort <- res_cibersort %>% as.matrix

#--------------------------------------------------------------------------
## Create Big Table with all values

	## input metadata  ## 2 columns: samplename condition #### rows should be sorted according to case and control so its levels of factors
metadata <- read.table(args[3], header=T)

dfs <- list(res_abis, res_base, res_consensus_tme, res_dcq, res_epic, res_mMCPcounter_m38, res_quantiseq, res_seqimmucc, res_timer, res_xCell, res_cibersort)
names(dfs) <- c("res_abis", "res_base", "res_consensus_tme", "res_dcq", "res_epic", "res_mMCPcounter_m38", "res_quantiseq", "res_seqimmucc", "res_timer", "res_xCell", "res_cibersort")

data_all="";
for (i in 1:length(dfs)){
        data1 <- dfs[i] %>% reshape2::melt(., id.vars = rownames(.), variable.name = "value")
	data1 <- merge(data1, metadata, by.x="Var2", by.y="samplename")
	data1$Var2 <- factor(data1$Var2, levels=metadata$samplename)
	
	## Bar plots
	var_height <- 10+ ((nrow(dfs[[i]])/10)*2)  ##  use [[. nrow(dfs[[i]]) to get content if dfs[i]; do typeof()
	print(var_height)
	pdf(paste(names(dfs[i]), "_barplot.pdf", sep=""), width=15, height=var_height)
	p <- ggplot(data1 %>% .[complete.cases(.), ], aes(Var2, value, fill=condition)) + geom_bar(stat='identity') + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  facet_wrap(~Var1, scales="free_y", ncol=4)
	print(p)
	dev.off()

	## Box plots case and control and pvalues
	comparisons <- combn(metadata$condition %>% unique,2, simplify=F)
	pdf(paste(names(dfs[i]), "_boxplot.pdf", sep=""), width=15, height=var_height)
	p <- ggboxplot(data1, x = "condition", y = "value", color = "condition", add = "jitter") + stat_compare_means(comparisons = comparisons, method = "t.test")
        q <- facet(p, facet.by = "Var1", scales="free_y", ncol=4) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +theme(axis.text.x = element_text(angle = 25, hjust = 1))
	print(q)
	dev.off()

	data_all <- rbind(data_all, data1)
}

write.table(data_all %>% .[complete.cases(.), ] , "Immune_deconv_allTools_bigTable.txt", quote=F, sep="\t", row.names=F)

#--------------------------------------------------------------------------
#### facet grid for diffrent cell types across tools: Bcells, Tcells, Macrophages
  ##Barplots:  between sample comparison Tools: xCell, TIMER, ConsensusTME, ESTIMATE, ABIS, mMCP-counter (mouse based), BASE (mouse based)
  ##Pie  Charts: between-cell-type comparisons: CIBERSORT, DCQ (mouse based)
  ## both EPIC, quanTIseq, CIBERSORT abs. mode, seqImmuCC (mouse based)
  ### T cell CD8+, T cell CD4+, Macrophage, Neutrophil, mDC, Monocyte, NK cell, BCell
#
big_mat3 <- read_table("")
facet_table <- read.table("cell_types_ForFacet_barplot_pieplot.txt", header=T, sep="\t")
facet_tab <- merge(big_mat3, facet_table[, c("CT_Tool", "Barplot", "Piechart", "MajorCellTypes", "Barplot_celltypes")], by="CT_Tool")


pdf("facet_barplot.pdf", height=20, width=20)
ggplot(facet_tab %>% .[.$MajorCellTypes == "yes",], aes(Var2, value, fill = condition)) + geom_bar(stat='identity') + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  facet_grid(L1~Barplot_celltypes, scales="free_y")
dev.off()

## Heatmaps

cell_types <- data.frame(Var1 = data_all$Var1, CellTypes=chartr(" ", "_", (data_all$Var1))) %>% unique
cell_annot <- read.table("cell_annot.txt", header=T)
cell_types_annot <- merge(cell_types, cell_annot, by.x="CellTypes", by.y="newCellType", all.x=T)
write.table(cell_types_annot, "cell_types_annot.txt", quote=F, row.names=F, sep="\t")

### Edit cell annotations and create broad categories
library(scales)
library(broom)
library(reshape2)

cell_annot <- read.table("cell_types_annot.txt", header=T, row.names=1, sep="\t")
big_mat <- merge(data_all, cell_annot, by="Var1", all.x=T)

big_mat3 <- big_mat %>% .[complete.cases(.), ] %>% group_by(L1) %>% mutate(rescaleMinMaxTool = rescale(as.numeric(as.character(value)), to=c(0,1))) %>% data.frame() %>% group_by(L1, CellTypes) %>% mutate(rescaleMinMaxCellTypeTool = rescale(as.numeric(as.character(value)), to=c(0,1))) %>% data.frame() %>% .$rescaleMinMaxTool

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


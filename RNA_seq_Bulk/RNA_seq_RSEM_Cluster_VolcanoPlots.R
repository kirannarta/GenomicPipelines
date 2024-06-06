#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

## code for DE analysis using DEseq2

# load libraries
library(DESeq2)
library(tximport)
library("tidyverse")
library(ggrepel)
library(here) ## to create dir
library(data.table) ## for rbindlist
library(scales)


###Usage: Rscript RNA_seq_RSEM_Cluster_VolcanoPlots.R ../metadata_Kurt_Samples1.txt ../RSEM_files/

#############################3 Import required files : metadata and RSEM files
metafile = args[1]

## input files : metadata (samplename, condition), RSEM files
metadata <- read.table(metafile, sep="\t", header=T)

## RSEM files ## the "samplename column should be same as names before "_rsem.genes.results" in RSEM files"
rsemfiles_path = args[2] ## eg in command line "../InputFiles/STING_RSEMFiles/"
files <- file.path(rsemfiles_path, paste0(metadata$samplename, "_rsem.genes.results"))
names(files) <- metadata$samplename
print(files)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

txi.rsem$abundance <- txi.rsem$abundance[apply(txi.rsem$length,1, function(row) all(row !=0 )),]
txi.rsem$counts <- txi.rsem$counts[apply(txi.rsem$length,1, function(row) all(row !=0 )),]
txi.rsem$length <- txi.rsem$length[apply(txi.rsem$length,1, function(row) all(row !=0 )),]

rownames(metadata) <- metadata$samplename

write.table(txi.rsem, "Combined_RSEM_expected_counts.txt", sep="\t", quote=F, row.names=F)


############################## Preprocessing: Create hclust, PCA 
## http://www.bea.ki.se/documents/Intro2RNAseq.pdf
if (!dir.exists(here("Preprocessing"))) {dir.create(here("Preprocessing"))}

dds = DESeqDataSetFromTximport(txi.rsem, metadata, ~condition)

DESeq.ds <- dds[ rowSums(counts(dds)) > 0, ]
DESeq.ds <- estimateSizeFactors(DESeq.ds)

sizeFactors(DESeq.ds)
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE)
log.norm.counts <- log2(counts.sf_normalized + 1)

DESeq.rlog <- vst(DESeq.ds, blind = TRUE)
rlog.norm.counts <- assay(DESeq.rlog)

	## Normalise Counts plots
pdf("Preprocessing/Normalised.pdf", width=10)
par(mfrow=c(2,1))
boxplot(counts.sf_normalized, notch = TRUE,main = "untransformed read counts", ylab = "read counts")
boxplot(log.norm.counts, notch = TRUE,main = "log2-transformed read counts",ylab = "log2(read counts)", las=2)
dev.off()

	#### Clustering plot

distance.m_rlog <- as.dist(1 - cor(rlog.norm.counts, method = "pearson" ))
pdf("Preprocessing/hclust_dendrogram.pdf")
plot( hclust(distance.m_rlog),labels = colnames(rlog.norm.counts),main = "rlog transformed read counts\ndistance: Pearson correlation")
dev.off()

	### PCA plot
P <- plotPCA(DESeq.rlog)
P <- P + theme_classic() + ggtitle("Rlog transformed counts")

pdf("Preprocessing/PCA_unlabelled.pdf")
print(P)
dev.off()

	### PCA plot Labelled, may change PCs if needed

plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE)
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }

  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
                  intgroup.df, name = rownames(colData(object)))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  p <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + theme_classic() + geom_text_repel(size=3)
  pdf("Preprocessing/PCA_PC1vsPC2_labelled.pdf", width=10)
  print(p)
  dev.off()
}
plotPCA.san(DESeq.rlog, ntop = 500, returnData = FALSE)

	#### Print PC loadings
##### Loadings

pca <- prcomp(t(rlog.norm.counts))
loadings <- as.data.frame(pca$rotation)

PC1_top500 <- loadings %>% mutate(Genes = map_chr(str_split(rownames(loadings), "_"),2))  %>% .[, c("Genes", "PC1")] %>% .[order(-.$PC1),] %>% head(500)  %>% rbind(c(">PC1_top500", "PC1top500") , .)
PC1_bottom500 <- loadings %>% mutate(Genes = map_chr(str_split(rownames(loadings), "_"),2))  %>% .[, c("Genes", "PC1")] %>% .[order(.$PC1),] %>% head(500)  %>% rbind(c(">PC1_bottom500", "PC1bottom500") , .) 
PC2_top500 <- loadings %>% mutate(Genes = map_chr(str_split(rownames(loadings), "_"),2))  %>% .[, c("Genes", "PC2")] %>% .[order(-.$PC2),] %>% head(500)  %>% rbind(c(">PC2_top500", "PC2top500") , .)
PC2_bottom500 <- loadings %>% mutate(Genes = map_chr(str_split(rownames(loadings), "_"),2))  %>% .[, c("Genes", "PC2")] %>% .[order(.$PC2),] %>% head(500)  %>% rbind(c(">PC2_bottom500", "PC2bottom500") , .) 
PC3_top500 <- loadings %>% mutate(Genes = map_chr(str_split(rownames(loadings), "_"),2))  %>% .[, c("Genes", "PC3")] %>% .[order(-.$PC3),] %>% head(500)  %>% rbind(c(">PC3_top500", "PC3top500") , .) 
PC3_bottom500 <- loadings %>% mutate(Genes = map_chr(str_split(rownames(loadings), "_"),2))  %>% .[, c("Genes", "PC3")] %>% .[order(.$PC3),] %>% head(500)  %>% rbind(c(">PC3_bottom500", "PC3bottom500") , .) 

PCA_loading_Gene_list <- rbindlist(list(PC1_top500, PC1_bottom500,PC2_top500 ,  PC2_bottom500,PC3_top500, PC3_bottom500), use.names=FALSE)
write.table(PCA_loading_Gene_list, "Preprocessing/PCA_Loadings_top-bottomGenesPC123.txt", row.names=F, quote=F, sep="\t")


############################## DE: Create Volcano plots DEseq and generate DE tables

dds = DESeqDataSetFromTximport(txi.rsem, metadata, ~condition)
dds <- DESeq(dds)
de_vol <- function(Cond1, Cond2, pval_label, FClabel) {
	res <- results(dds, contrast=c("condition",Cond1,Cond2))
	res_Vol <- res %>% data.frame %>% .[complete.cases(.),] %>% mutate(threshold_OE = abs(.$log2FoldChange) > 1 & .$padj < 0.05)
	res_Vol$Gene_ENSG <- rownames(res_Vol)
	res_Vol$GeneSymbol <- map_chr(str_split(res_Vol$Gene_ENSG, "_"),2)
	### unlabelled volcano
	pdf(paste("res_", Cond1, "-vs-", Cond2, "_DEseq2_allgenes_DE_nolabel.pdf", sep=""), width=10, height=10)
	
	p <- ggplot(res_Vol, aes(x = ifelse(log2FoldChange >5, 5, ifelse(log2FoldChange < -5, -5, log2FoldChange)), y = ifelse(-log10(padj) < 50, -log10(padj), 50), colour = threshold_OE)) + geom_point() +
        ggtitle(paste("DE Volcano: ", Cond1, " -vs- ", Cond2, sep="")) +
        xlab("log2 fold change") +                                                                                                             
        ylab("-log10 adjusted p-value") +
	scale_colour_manual(values = c("TRUE"='red3', "FALSE"='grey'), name="-log10(padj)<=0.05 \n abs|log2FC| >1", labels=c("DE","notDE")) +
        theme_classic() + theme(axis.text=element_text(size=20)) +
        geom_vline(xintercept=c(-1, 1), col="grey", linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="grey")
	
	print(p)
	dev.off()

	### labelled volcano
	label_top50 <- res_Vol %>% .[abs(.$log2FoldChange) > 1,] %>% arrange(-desc(padj)) %>% slice(1:50) %>% .$GeneSymbol

	pdf(paste("res_", Cond1, "-vs-", Cond2, "_DEseq2_allgenes_DE_labelled_ylim100_labelpval30.pdf", sep=""), width=10, height=10)
        
        p <- ggplot(res_Vol, aes(x = ifelse(log2FoldChange >5, 5, ifelse(log2FoldChange < -5, -5, log2FoldChange)), y = ifelse(-log10(padj) < 100, -log10(padj), 100), colour = threshold_OE)) + geom_point() +
        ggtitle(paste("DE Volcano: ", Cond1, " -vs- ", Cond2, sep="")) +
        xlab("log2 fold change") +
        ylab("-log10 adjusted p-value") +
	geom_text_repel(size=3, label = ifelse(res_Vol$GeneSymbol %in% label_top50, res_Vol$GeneSymbol, ""), max.overlaps=1000, min.segment.length = 0)+
	#geom_text_repel(size=3, label = ifelse((-log10(res_Vol$padj) > pval_label & abs(res_Vol$log2FoldChange) > FClabel) , res_Vol$GeneSymbol, ""), max.overlaps=1000, min.segment.length = 0) +
        scale_colour_manual(values = c("TRUE"='red3', "FALSE"='grey'), name="-log10(padj)<=0.05 \n abs|log2FC| >1", labels=c("DE","notDE")) +
        theme_classic() + theme(axis.text=element_text(size=20)) +
        geom_vline(xintercept=c(-1, 1), col="grey", linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="grey")
        
        print(p)
        dev.off()

	write.table(res_Vol, paste0(Cond1,"-vs-", Cond2, "_allgenes.txt"), sep="\t", quote=F, row.names=F)
	write.table(res_Vol  %>% .[.$padj< 0.05 & abs(.$log2FoldChange) > 1,], paste0(Cond1, "-vs-", Cond2, "_DEgenes.txt"), sep="\t", quote=F, row.names=F)
}


comparisons <- combn(metadata$condition %>% unique, 2, simplify=F)
for (i in 1:length(comparisons)){
	de_vol(comparisons[[i]][1], comparisons[[i]][2], 30, 1)
}


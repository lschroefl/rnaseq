# import all necessary libraries
library("DESeq2")
library("readr")
library("tximport")
library("edgeR")
library("httpgd")
library("tximport")
library("regionReport")
library("knitr")
library("printr") # would be needed to print help() in markdown html
library("ggplot2")
library("ggrepel")
library("tidyverse")
library("RColorBrewer")
library("pheatmap")
library("reshape2")


# set up httpgd graphics
hgd()
hgd_browse()
#dev.list()
#dev.off()

# -----------------------------IMPORT & MODEL-----------------------------------------

# GENERAL (conditions for DESeq2 and edgeR)
condition <- factor(c(rep("ctrl", 5), rep("trt", 5), rep("sham", 3)))
edgeR_design <- model.matrix(~condition)

# GENERAL (sample&condition table for salmon)
salmon_samples <- read.table("/projects/neuronet/seq_april/samples_april.csv", sep = ",", header = TRUE, row.names = 1)
salmon_samples$condition <- as.factor(salmon_samples$condition)

# GENERAL (import of countdata from the star mapped data)
ensembl_raw_star_feature_countdata=as.matrix(read.table("/projects/neuronet/seq_april/processing/mapping/ensembl/raw/star/featureCounts_countdata", header = TRUE, row.names=1))
ensembl_trim_star_feature_countdata=as.matrix(read.table("/projects/neuronet/seq_april/processing/mapping/ensembl/trimmomatic/star/featureCounts_countdata", header = TRUE, row.names=1))
ensembl_raw_star_htseq_countdata=as.matrix(read.table("/projects/neuronet/seq_april/processing/mapping/ensembl/raw/star/htseqCount_countdata", header = TRUE, row.names=1))
ensembl_trim_star_htseq_countdata=as.matrix(read.table("/projects/neuronet/seq_april/processing/mapping/ensembl/trimmomatic/star/htseqCount_countdata", header = TRUE, row.names=1))

# GENERAL (import of salmon data via tximport)
ensembl_tx2gene <- read_csv("genome/ensembl/tx2gene_GRCm39.109.csv", col_names = FALSE)
ensembl_raw_salmon_files <- file.path("mapping/ensembl/raw/salmon", rownames(salmon_samples), "quant.sf")
ensembl_trim_salmon_files <- file.path("mapping/ensembl/trimmomatic/salmon", rownames(salmon_samples), "quant.sf")
names(ensembl_raw_salmon_files) <- rownames(salmon_samples)
names(ensembl_trim_salmon_files) <- rownames(salmon_samples)
ensembl_raw_salmon_txi <- tximport(ensembl_raw_salmon_files, type = "salmon", tx2gene = ensembl_tx2gene)
ensembl_trim_salmon_txi <- tximport(ensembl_trim_salmon_files, type = "salmon", tx2gene = ensembl_tx2gene)

# DESeq MAPPED (star)
ensembl_raw_star_feature_DESeq_coldata=data.frame(row.names = colnames(ensembl_raw_star_feature_countdata), condition)
ensembl_trim_star_feature_DESeq_coldata=data.frame(row.names = colnames(ensembl_trim_star_feature_countdata), condition)
ensembl_raw_star_htseq_DESeq_coldata=data.frame(row.names = colnames(ensembl_raw_star_htseq_countdata), condition)
ensembl_trim_star_htseq_DESeq_coldata=data.frame(row.names = colnames(ensembl_trim_star_htseq_countdata), condition)
ensembl_raw_star_feature_dds <- DESeqDataSetFromMatrix( ensembl_raw_star_feature_countdata, ensembl_raw_star_feature_DESeq_coldata, design = ~condition)
ensembl_trim_star_feature_dds <- DESeqDataSetFromMatrix( ensembl_trim_star_feature_countdata, ensembl_trim_star_feature_DESeq_coldata, design = ~condition)
ensembl_raw_star_htseq_dds <- DESeqDataSetFromMatrix( ensembl_raw_star_htseq_countdata, ensembl_raw_star_htseq_DESeq_coldata, design = ~condition)
ensembl_trim_star_htseq_dds <- DESeqDataSetFromMatrix( ensembl_trim_star_htseq_countdata, ensembl_trim_star_htseq_DESeq_coldata, design = ~condition)
ensembl_raw_star_feature_dds$condition <- relevel(ensembl_raw_star_feature_dds$condition, ref = "ctrl")
ensembl_trim_star_feature_dds$condition <- relevel(ensembl_trim_star_feature_dds$condition, ref = "ctrl")
ensembl_raw_star_htseq_dds$condition <- relevel(ensembl_raw_star_feature_dds$condition, ref = "ctrl")
ensembl_trim_star_htseq_dds$condition <- relevel(ensembl_raw_star_feature_dds$condition, ref = "ctrl")
ensembl_raw_star_feature_dds_keep <- rowSums(counts(ensembl_raw_star_feature_dds) >= 5) >= 3
ensembl_trim_star_feature_dds_keep <- rowSums(counts(ensembl_raw_star_feature_dds) >= 5) >= 3
ensembl_raw_star_htseq_dds_keep <- rowSums(counts(ensembl_raw_star_feature_dds) >= 5) >= 3
ensembl_trim_star_htseq_dds_keep <- rowSums(counts(ensembl_raw_star_feature_dds) >= 5) >= 3
ensembl_raw_star_feature_dds <- ensembl_raw_star_feature_dds[ensembl_raw_star_feature_dds_keep,]
ensembl_trim_star_feature_dds <- ensembl_trim_star_feature_dds[ensembl_trim_star_feature_dds_keep,]
ensembl_raw_star_htseq_dds <- ensembl_raw_star_htseq_dds[ensembl_raw_star_htseq_dds_keep,]
ensembl_trim_star_htseq_dds <- ensembl_trim_star_htseq_dds[ensembl_trim_star_htseq_dds_keep,]
ensembl_raw_star_feature_dds <- DESeq(ensembl_raw_star_feature_dds)
ensembl_trim_star_feature_dds <- DESeq(ensembl_trim_star_feature_dds)
ensembl_raw_star_htseq_dds <- DESeq(ensembl_raw_star_htseq_dds)
ensembl_trim_star_htseq_dds <- DESeq(ensembl_trim_star_htseq_dds)

# edgeR MAPPED (star)
ensembl_raw_star_feature_edgeR_y <- DGEList(counts=ensembl_raw_star_feature_countdata, group=condition)
ensembl_trim_star_feature_edgeR_y <- DGEList(counts=ensembl_trim_star_feature_countdata, group=condition)
ensembl_raw_star_htseq_edgeR_y <- DGEList(counts=ensembl_raw_star_htseq_countdata, group=condition)
ensembl_trim_star_htseq_edgeR_y <- DGEList(counts=ensembl_trim_star_htseq_countdata, group=condition)
ensembl_raw_star_feature_edgeR_keep=filterByExpr(ensembl_raw_star_feature_edgeR_y)
ensembl_trim_star_feature_edgeR_keep=filterByExpr(ensembl_trim_star_feature_edgeR_y)
ensembl_raw_star_htseq_edgeR_keep=filterByExpr(ensembl_raw_star_htseq_edgeR_y)
ensembl_trim_star_htseq_edgeR_keep=filterByExpr(ensembl_trim_star_htseq_edgeR_y)
ensembl_raw_star_feature_edgeR_y <- ensembl_raw_star_feature_edgeR_y[ensembl_raw_star_feature_edgeR_keep,, keep.lib.sizes=FALSE]
ensembl_trim_star_feature_edgeR_y <- ensembl_trim_star_feature_edgeR_y[ensembl_trim_star_feature_edgeR_keep,, keep.lib.sizes=FALSE]
ensembl_raw_star_htseq_edgeR_y <- ensembl_raw_star_htseq_edgeR_y[ensembl_raw_star_htseq_edgeR_keep,, keep.lib.sizes=FALSE]
ensembl_trim_star_htseq_edgeR_y <- ensembl_trim_star_htseq_edgeR_y[ensembl_trim_star_htseq_edgeR_keep,, keep.lib.sizes=FALSE]
ensembl_raw_star_feature_edgeR_y <- calcNormFactors(ensembl_raw_star_feature_edgeR_y)
ensembl_trim_star_feature_edgeR_y <- calcNormFactors(ensembl_trim_star_feature_edgeR_y)
ensembl_raw_star_htseq_edgeR_y <- calcNormFactors(ensembl_raw_star_htseq_edgeR_y)
ensembl_trim_star_htseq_edgeR_y <- calcNormFactors(ensembl_trim_star_htseq_edgeR_y)
ensembl_raw_star_feature_edgeR_y <- estimateDisp(ensembl_raw_star_feature_edgeR_y, edgeR_design)
ensembl_trim_star_feature_edgeR_y <- estimateDisp(ensembl_trim_star_feature_edgeR_y, edgeR_design)
ensembl_raw_star_htseq_edgeR_y <- estimateDisp(ensembl_raw_star_htseq_edgeR_y, edgeR_design)
ensembl_trim_star_htseq_edgeR_y <- estimateDisp(ensembl_trim_star_htseq_edgeR_y, edgeR_design)
ensembl_raw_star_feature_edgeR_fit <- glmQLFit(ensembl_raw_star_feature_edgeR_y, edgeR_design)
ensembl_trim_star_feature_edgeR_fit <- glmQLFit(ensembl_trim_star_feature_edgeR_y, edgeR_design)
ensembl_raw_star_htseq_edgeR_fit <- glmQLFit(ensembl_raw_star_htseq_edgeR_y, edgeR_design)
ensembl_trim_star_htseq_edgeR_fit <- glmQLFit(ensembl_trim_star_htseq_edgeR_y, edgeR_design)
ensembl_raw_star_feature_edgeR_qlf <- glmQLFTest(ensembl_raw_star_feature_edgeR_fit, coef="conditiontrt")
ensembl_trim_star_feature_edgeR_qlf <- glmQLFTest(ensembl_trim_star_feature_edgeR_fit, coef="conditiontrt")
ensembl_raw_star_htseq_edgeR_qlf <- glmQLFTest(ensembl_raw_star_htseq_edgeR_fit, coef="conditiontrt")
ensembl_trim_star_htseq_edgeR_qlf <- glmQLFTest(ensembl_trim_star_htseq_edgeR_fit, coef="conditiontrt")


# DESeq Pseudomapped (salmon)
ensembl_raw_salmon_dds <- DESeqDataSetFromTximport( ensembl_raw_salmon_txi, colData = salmon_samples , design = ~condition)
ensembl_trim_salmon_dds <- DESeqDataSetFromTximport( ensembl_trim_salmon_txi, colData = salmon_samples, design = ~condition)
ensembl_raw_salmon_dds$condition <- relevel(ensembl_raw_salmon_dds$condition, ref = "ctrl")
ensembl_trim_salmon_dds$condition <- relevel(ensembl_trim_salmon_dds$condition, ref = "ctrl")
ensembl_raw_salmon_dds_keep <- rowSums(counts(ensembl_raw_salmon_dds) >= 5) >= 3
ensembl_trim_salmon_dds_keep <- rowSums(counts(ensembl_trim_salmon_dds) >= 5) >= 3
ensembl_raw_salmon_dds <- ensembl_raw_salmon_dds[ensembl_raw_salmon_dds_keep,]
ensembl_trim_salmon_dds <- ensembl_trim_salmon_dds[ensembl_trim_salmon_dds_keep,]
ensembl_raw_salmon_dds <- DESeq(ensembl_raw_salmon_dds)
ensembl_trim_salmon_dds <- DESeq(ensembl_trim_salmon_dds)


# edgeR pseudomapped data
ensembl_raw_salmon_txi_cts <- ensembl_raw_salmon_txi$counts
ensembl_raw_salmon_txi_normMat <- ensembl_raw_salmon_txi$length 
ensembl_raw_salmon_txi_normMat <- ensembl_raw_salmon_txi_normMat/exp(rowMeans(log(ensembl_raw_salmon_txi_normMat)))
ensembl_raw_salmon_txi_normCts <- ensembl_raw_salmon_txi_cts/ensembl_raw_salmon_txi_normMat
ensembl_raw_salmon_txi_eff.lib <- calcNormFactors(ensembl_raw_salmon_txi_normCts)*colSums(ensembl_raw_salmon_txi_normCts)
ensembl_raw_salmon_txi_normMat <- sweep(ensembl_raw_salmon_txi_normMat, 2, ensembl_raw_salmon_txi_eff.lib, "*")
ensembl_raw_salmon_txi_normMat <- log(ensembl_raw_salmon_txi_normMat)
ensembl_raw_salmon_edgeR_y <- DGEList(ensembl_raw_salmon_txi_cts)
ensembl_raw_salmon_edgeR_y <- scaleOffset(ensembl_raw_salmon_edgeR_y, ensembl_raw_salmon_txi_normMat)
ensembl_raw_salmon_edgeR_y_keep <- filterByExpr(ensembl_raw_salmon_edgeR_y, edgeR_design)
ensembl_raw_salmon_edgeR_y <- ensembl_raw_salmon_edgeR_y[ensembl_raw_salmon_edgeR_y_keep,]
ensembl_raw_salmon_edgeR_y <- calcNormFactors(ensembl_raw_salmon_edgeR_y)
ensembl_raw_salmon_edgeR_y <- estimateDisp(ensembl_raw_salmon_edgeR_y, edgeR_design)
ensembl_raw_salmon_edgeR_fit <- glmQLFit(ensembl_raw_salmon_edgeR_y, edgeR_design)
ensembl_raw_salmon_edgeR_qlf <- glmQLFTest(ensembl_raw_salmon_edgeR_fit, coef = "conditiontrt")

ensembl_trim_salmon_txi_cts <- ensembl_trim_salmon_txi$counts
ensembl_trim_salmon_txi_normMat <- ensembl_trim_salmon_txi$length 
ensembl_trim_salmon_txi_normMat <- ensembl_trim_salmon_txi_normMat/exp(rowMeans(log(ensembl_trim_salmon_txi_normMat)))
ensembl_trim_salmon_txi_normCts <- ensembl_trim_salmon_txi_cts/ensembl_trim_salmon_txi_normMat
ensembl_trim_salmon_txi_eff.lib <- calcNormFactors(ensembl_trim_salmon_txi_normCts)*colSums(ensembl_trim_salmon_txi_normCts)
ensembl_trim_salmon_txi_normMat <- sweep(ensembl_trim_salmon_txi_normMat, 2, ensembl_trim_salmon_txi_eff.lib, "*")
ensembl_trim_salmon_txi_normMat <- log(ensembl_trim_salmon_txi_normMat)
ensembl_trim_salmon_edgeR_y <- DGEList(ensembl_trim_salmon_txi_cts)
ensembl_trim_salmon_edgeR_y <- scaleOffset(ensembl_trim_salmon_edgeR_y, ensembl_trim_salmon_txi_normMat)
ensembl_trim_salmon_edgeR_y_keep <- filterByExpr(ensembl_trim_salmon_edgeR_y, edgeR_design)
ensembl_trim_salmon_edgeR_y <- ensembl_trim_salmon_edgeR_y[ensembl_trim_salmon_edgeR_y_keep,]
ensembl_trim_salmon_edgeR_y <- calcNormFactors(ensembl_trim_salmon_edgeR_y)
ensembl_trim_salmon_edgeR_y <- estimateDisp(ensembl_trim_salmon_edgeR_y, edgeR_design)
ensembl_trim_salmon_edgeR_fit <- glmQLFit(ensembl_trim_salmon_edgeR_y, edgeR_design)
ensembl_trim_salmon_edgeR_qlf <- glmQLFTest(ensembl_trim_salmon_edgeR_fit, coef = "conditiontrt")



## -------------RESULTS ************** RESULTS -----------

#----DESeq mapped data

ensembl_raw_star_feature_DESeq_Results <- results(ensembl_raw_star_feature_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ensembl_raw_star_feature_DESeq_ResSig <- subset(ensembl_raw_star_feature_DESeq_Results, padj < 0.1)
ensembl_raw_star_feature_DESeq_ResSig_DF <- as.data.frame(ensembl_raw_star_feature_DESeq_ResSig)
#sink("./results/june/sum/ensembl_raw_star_feature_DESeq_Results_summary.txt")
#summary(ensembl_raw_star_feature_DESeq_Results)
#sink()
#write.csv(ensembl_raw_star_feature_DESeq_ResSig_DF, "./results/june/tables/ensembl_raw_star_feature_DESeq_ResSig_DF.csv", row.names = TRUE)

ensembl_trim_star_feature_DESeq_Results <- results(ensembl_trim_star_feature_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ensembl_trim_star_feature_DESeq_ResSig <- subset(ensembl_trim_star_feature_DESeq_Results, padj < 0.1)
ensembl_trim_star_feature_DESeq_ResSig_DF <- as.data.frame(ensembl_trim_star_feature_DESeq_ResSig)
#sink("./results/june/sum/ensembl_trim_star_feature_DESeq_Results_summary.txt")
#summary(ensembl_trim_star_feature_DESeq_Results)
#sink()
#write.csv(ensembl_trim_star_feature_DESeq_ResSig_DF, "./results/june/tables/ensembl_trim_star_feature_DESeq_ResSig_DF.csv", row.names = TRUE)


ensembl_raw_star_htseq_DESeq_Results <- results(ensembl_raw_star_htseq_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ensembl_raw_star_htseq_DESeq_ResSig <- subset(ensembl_raw_star_htseq_DESeq_Results, padj < 0.1)
ensembl_raw_star_htseq_DESeq_ResSig_DF <- as.data.frame(ensembl_raw_star_htseq_DESeq_ResSig)
#sink("./results/june/sum/ensembl_raw_star_htseq_DESeq_Results_summary.txt")
#summary(ensembl_raw_star_htseq_DESeq_Results)
#sink()
#write.csv(ensembl_raw_star_htseq_DESeq_ResSig_DF, "./results/june/tables/ensembl_raw_star_htseq_DESeq_ResSig_DF.csv", row.names = TRUE)



ensembl_trim_star_htseq_DESeq_Results <- results(ensembl_trim_star_htseq_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ensembl_trim_star_htseq_DESeq_ResSig <- subset(ensembl_trim_star_htseq_DESeq_Results, padj < 0.1)
ensembl_trim_star_htseq_DESeq_ResSig_DF <- as.data.frame(ensembl_trim_star_htseq_DESeq_ResSig)
#sink("./results/june/sum/ensembl_trim_star_htseq_DESeq_Results_summary.txt")
#summary(ensembl_trim_star_htseq_DESeq_Results)
#sink()
#write.csv(ensembl_trim_star_htseq_DESeq_ResSig_DF, "./results/june/tables/ensembl_trim_star_htseq_DESeq_ResSig_DF.csv", row.names = TRUE)

# DESeq2 pseudomapped data

ensembl_raw_salmon_DESeq_Results <- results(ensembl_raw_salmon_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ensembl_raw_salmon_DESeq_ResSig <- subset(ensembl_raw_salmon_DESeq_Results, padj < 0.1)
ensembl_raw_salmon_DESeq_ResSig_DF <- as.data.frame(ensembl_raw_salmon_DESeq_ResSig)
#sink("./results/june/sum/ensembl_raw_salmon_DESeq_Results_summary.txt")
#summary(ensembl_raw_salmon_DESeq_Results)
#sink()
#write.csv(ensembl_raw_salmon_DESeq_ResSig_DF, "./results/june/tables/ensembl_raw_salmon_DESeq_ResSig_DF.csv", row.names = TRUE)



ensembl_trim_salmon_DESeq_Results <- results(ensembl_trim_salmon_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ensembl_trim_salmon_DESeq_ResSig <- subset(ensembl_trim_salmon_DESeq_Results, padj < 0.1)
ensembl_trim_salmon_DESeq_ResSig_DF <- as.data.frame(ensembl_trim_salmon_DESeq_ResSig)
#sink("./results/june/sum/ensembl_trim_salmon_DESeq_Results_summary.txt")
#summary(ensembl_trim_salmon_DESeq_Results)
#sink()
#write.csv(ensembl_trim_salmon_DESeq_ResSig_DF, "./results/june/tables/ensembl_trim_salmon_DESeq_ResSig_DF.csv", row.names = TRUE)



#----------RESULTS edgeR mapped data

ensembl_raw_star_feature_edgeR_topTags <- topTags(ensembl_raw_star_feature_edgeR_qlf, p.value = 0.1)
ensembl_raw_star_feature_edgeR_topTags_DF <- as.data.frame(ensembl_raw_star_feature_edgeR_topTags)
#sink("./results/june/sum/ensembl_raw_star_feature_edgeR_qlf_summary.txt")
#summary(decideTests(ensembl_raw_star_feature_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ensembl_raw_star_feature_edgeR_topTags_DF, "./results/june/tables/ensembl_raw_star_feature_edgeR_topTags_DF.csv", row.names = TRUE)



ensembl_trim_star_feature_edgeR_topTags <- topTags(ensembl_trim_star_feature_edgeR_qlf, p.value = 0.1)
ensembl_trim_star_feature_edgeR_topTags_DF <- as.data.frame(ensembl_trim_star_feature_edgeR_topTags)
#sink("./results/june/sum/ensembl_trim_star_feature_edgeR_qlf_summary.txt")
#summary(decideTests(ensembl_trim_star_feature_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ensembl_trim_star_feature_edgeR_topTags_DF, "./results/june/tables/ensembl_trim_star_feature_edgeR_topTags_DF.csv", row.names = TRUE)


ensembl_raw_star_htseq_edgeR_topTags <- topTags(ensembl_raw_star_htseq_edgeR_qlf, p.value = 0.1)
ensembl_raw_star_htseq_edgeR_topTags_DF <- as.data.frame(ensembl_raw_star_htseq_edgeR_topTags)
#sink("./results/june/sum/ensembl_raw_star_htseq_edgeR_qlf_summary.txt")
#summary(decideTests(ensembl_raw_star_htseq_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ensembl_raw_star_htseq_edgeR_topTags_DF, "./results/june/tables/ensembl_raw_star_htseq_edgeR_topTags_DF.csv", row.names = TRUE)


ensembl_trim_star_htseq_edgeR_topTags <- topTags(ensembl_trim_star_htseq_edgeR_qlf, p.value = 0.1)
ensembl_trim_star_htseq_edgeR_topTags_DF <- as.data.frame(ensembl_trim_star_htseq_edgeR_topTags)
#sink("./results/june/sum/ensembl_trim_star_htseq_edgeR_qlf_summary.txt")
#summary(decideTests(ensembl_trim_star_htseq_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ensembl_trim_star_htseq_edgeR_topTags_DF, "./results/june/tables/ensembl_trim_star_htseq_edgeR_topTags_DF.csv", row.names = TRUE)

# edgeR pseudomapped data

ensembl_raw_salmon_edgeR_topTags <- topTags(ensembl_raw_salmon_edgeR_qlf, p.value = 0.1)
ensembl_raw_salmon_edgeR_topTags_DF <- as.data.frame(ensembl_raw_salmon_edgeR_topTags)
#sink("./results/june/sum/ensembl_raw_salmon_edgeR_qlf_summary.txt")
#summary(decideTests(ensembl_raw_salmon_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ensembl_raw_salmon_edgeR_topTags_DF, "./results/june/tables/ensembl_raw_salmon_edgeR_topTags_DF.csv", row.names = TRUE)


ensembl_trim_salmon_edgeR_topTags <- topTags(ensembl_trim_salmon_edgeR_qlf, p.value = 0.1)
ensembl_trim_salmon_edgeR_topTags_DF <- as.data.frame(ensembl_trim_salmon_edgeR_topTags)
#sink("./results/june/sum/ensembl_trim_salmon_edgeR_qlf_summary.txt")
#summary(decideTests(ensembl_trim_salmon_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ensembl_trim_salmon_edgeR_topTags_DF, "./results/june/tables/ensembl_trim_salmon_edgeR_topTags_DF.csv", row.names = TRUE)


###--------------------- visualization ***** visualization -----------------



# COMPARISON TO RAFAELS RESULTS

rafael_genes <- read.table("./results/rafael/rafael_genes.csv", sep = ",", header = TRUE, row.names = 1)


rafael_genes$geneName
genes_ensembl_salmon <- list(
    "rafael_genes" = rafael_genes$geneName,
    "ensembl_raw_salmon_DESeq" = row.names(ensembl_raw_salmon_DESeq_ResSig_DF),
    "ensembl_trim_salmon_DESeq" = row.names(ensembl_trim_salmon_DESeq_ResSig_DF),
    "ensembl_raw_salmon_edgeR" = row.names(ensembl_raw_salmon_edgeR_topTags_DF), 
    "ensembl_trim_salmon_edgeR" = row.names(ensembl_trim_salmon_edgeR_topTags_DF)
)

jpeg(file="./results/june/plots/genes_ensembl_salmon.jpeg", width = 1500, height = 1500)
elements <- ggvenn(genes_ensembl_salmon, c("rafael_genes", "ensembl_raw_salmon_DESeq", "ensembl_trim_salmon_DESeq"), 
    show_elements = TRUE, label_sep = "\n", show_percentage = TRUE, digits = 2)
percentage <- ggvenn(genes_ensembl_salmon, c("rafael_genes", "ensembl_raw_salmon_DESeq", "ensembl_trim_salmon_DESeq"), 
    show_percentage = TRUE, digits = 0, text_size = 6)
ggarrange(percentage, elements, labels = c("Percent", "Elements"),
ncol = 2, nrow = 1)
dev.off()


# data transformation for visualization
ensembl_raw_salmon_dds_rlog <- rlog(ensembl_raw_salmon_dds, blind = TRUE)
ensembl_raw_salmon_dds_lfcShrink <- lfcShrink(ensembl_raw_salmon_dds, coef = "condition_trt_vs_ctrl", type = "apeglm")


# creating a html report

DESeq2Report(
    ensembl_raw_salmon_dds, 
    res = ensembl_raw_salmon_DESeq_Results, 
    intgroup = "condition", 
    nBest = sum(ensembl_raw_salmon_DESeq_Results$padj < 0.1, na.rm = TRUE), 
    outdir = "./results/june/html", 
    output = "genes_ensembl_raw_salmon", 
    customCode = "/projects/neuronet/seq_april/scripts/genes_ensembl_raw_salmon.Rmd"
    )

knit2html(
    "/projects/neuronet/seq_april/scripts/genes_ensembl_raw_salmon.Rmd", 
    output = "/projects/neuronet/seq_april/processing/results/june/html/test.html")


#-------------------------



rmarkdown::render("/projects/neuronet/seq_april/scripts/genes_ensembl_raw_salmon.md",
    output_format = "html_document",
    output_dir = "/projects/neuronet/seq_april/processing/results/june/html/",
    clean = FALSE)

# my plots  
# change the axis title of "number of assigned reads"
# na.omit to not have the warnings 

nrows(ensembl_raw_salmon_dds)
ensembl_raw_salmon_dds



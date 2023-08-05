# import all necessary libraries
library("DESeq2")
library("readr")
library("tximport")
library("edgeR")
library("httpgd")
library("tximport")
library("ggvenn")
library("ggpubr")


# set up graphics
hgd()
hgd_browse()
dev.list()
dev.off()

# -----------------------------IMPORT & MODEL-----------------------------------------

# GENERAL (conditions for DESeq2)
condition <- factor(c(rep("ctrl", 5), rep("trt", 5), rep("sham", 3)))
edgeR_design <- model.matrix(~condition)

# GENERAL (sample&condition table for salmon)
salmon_samples <- read.table("../samples_april.csv", sep = ",", header = TRUE, row.names = 1)
salmon_samples$condition <- as.factor(salmon_samples$condition)

# GENERAL (import of countdata from the star mapped data)
ucsc_raw_star_feature_countdata=as.matrix(read.table("./mapping/ucsc/raw/star/featureCounts_countdata", header = TRUE, row.names=1))
ucsc_trim_star_feature_countdata=as.matrix(read.table("./mapping/ucsc/trimmomatic/star/featureCounts_countdata", header = TRUE, row.names=1))
ucsc_raw_star_htseq_countdata=as.matrix(read.table("./mapping/ucsc/raw/star/htseqCount_countdata", header = TRUE, row.names=1))
ucsc_trim_star_htseq_countdata=as.matrix(read.table("./mapping/ucsc/trimmomatic/star/htseqCount_countdata", header = TRUE, row.names=1))

# GENERAL (import of salmon data via tximport)
ucsc_tx2gene <- read_csv("genome/ucsc/tx2gene_mm39.csv", col_names = FALSE)
ucsc_raw_salmon_files <- file.path("mapping/ucsc/raw/salmon", rownames(salmon_samples), "quant.sf")
ucsc_trim_salmon_files <- file.path("mapping/ucsc/trimmomatic/salmon", rownames(salmon_samples), "quant.sf")
names(ucsc_raw_salmon_files) <- rownames(salmon_samples)
names(ucsc_trim_salmon_files) <- rownames(salmon_samples)
ucsc_raw_salmon_txi <- tximport(ucsc_raw_salmon_files, type = "salmon", tx2gene = ucsc_tx2gene)
ucsc_trim_salmon_txi <- tximport(ucsc_trim_salmon_files, type = "salmon", tx2gene = ucsc_tx2gene)



#---------------- MODEL ********* MODEL ------------------

# DESeq MAPPED (star)
ucsc_raw_star_feature_DESeq_coldata=data.frame(row.names = colnames(ucsc_raw_star_feature_countdata), condition)
ucsc_trim_star_feature_DESeq_coldata=data.frame(row.names = colnames(ucsc_trim_star_feature_countdata), condition)
ucsc_raw_star_htseq_DESeq_coldata=data.frame(row.names = colnames(ucsc_raw_star_htseq_countdata), condition)
ucsc_trim_star_htseq_DESeq_coldata=data.frame(row.names = colnames(ucsc_trim_star_htseq_countdata), condition)
ucsc_raw_star_feature_dds <- DESeqDataSetFromMatrix( ucsc_raw_star_feature_countdata, ucsc_raw_star_feature_DESeq_coldata, design = ~condition)
ucsc_trim_star_feature_dds <- DESeqDataSetFromMatrix( ucsc_trim_star_feature_countdata, ucsc_trim_star_feature_DESeq_coldata, design = ~condition)
ucsc_raw_star_htseq_dds <- DESeqDataSetFromMatrix( ucsc_raw_star_htseq_countdata, ucsc_raw_star_htseq_DESeq_coldata, design = ~condition)
ucsc_trim_star_htseq_dds <- DESeqDataSetFromMatrix( ucsc_trim_star_htseq_countdata, ucsc_trim_star_htseq_DESeq_coldata, design = ~condition)
ucsc_raw_star_feature_dds$condition <- relevel(ucsc_raw_star_feature_dds$condition, ref = "ctrl")
ucsc_trim_star_feature_dds$condition <- relevel(ucsc_trim_star_feature_dds$condition, ref = "ctrl")
ucsc_raw_star_htseq_dds$condition <- relevel(ucsc_raw_star_feature_dds$condition, ref = "ctrl")
ucsc_trim_star_htseq_dds$condition <- relevel(ucsc_raw_star_feature_dds$condition, ref = "ctrl")
ucsc_raw_star_feature_dds_keep <- rowSums(counts(ucsc_raw_star_feature_dds) >= 5) >= 3
ucsc_trim_star_feature_dds_keep <- rowSums(counts(ucsc_raw_star_feature_dds) >= 5) >= 3
ucsc_raw_star_htseq_dds_keep <- rowSums(counts(ucsc_raw_star_feature_dds) >= 5) >= 3
ucsc_trim_star_htseq_dds_keep <- rowSums(counts(ucsc_raw_star_feature_dds) >= 5) >= 3
ucsc_raw_star_feature_dds <- ucsc_raw_star_feature_dds[ucsc_raw_star_feature_dds_keep,]
ucsc_trim_star_feature_dds <- ucsc_trim_star_feature_dds[ucsc_trim_star_feature_dds_keep,]
ucsc_raw_star_htseq_dds <- ucsc_raw_star_htseq_dds[ucsc_raw_star_htseq_dds_keep,]
ucsc_trim_star_htseq_dds <- ucsc_trim_star_htseq_dds[ucsc_trim_star_htseq_dds_keep,]
ucsc_raw_star_feature_dds <- DESeq(ucsc_raw_star_feature_dds)
ucsc_trim_star_feature_dds <- DESeq(ucsc_trim_star_feature_dds)
ucsc_raw_star_htseq_dds <- DESeq(ucsc_raw_star_htseq_dds)
ucsc_trim_star_htseq_dds <- DESeq(ucsc_trim_star_htseq_dds)

# edgeR MAPPED (star)
ucsc_raw_star_feature_edgeR_y <- DGEList(counts=ucsc_raw_star_feature_countdata, group=condition)
ucsc_trim_star_feature_edgeR_y <- DGEList(counts=ucsc_trim_star_feature_countdata, group=condition)
ucsc_raw_star_htseq_edgeR_y <- DGEList(counts=ucsc_raw_star_htseq_countdata, group=condition)
ucsc_trim_star_htseq_edgeR_y <- DGEList(counts=ucsc_trim_star_htseq_countdata, group=condition)
ucsc_raw_star_feature_edgeR_keep=filterByExpr(ucsc_raw_star_feature_edgeR_y)
ucsc_trim_star_feature_edgeR_keep=filterByExpr(ucsc_trim_star_feature_edgeR_y)
ucsc_raw_star_htseq_edgeR_keep=filterByExpr(ucsc_raw_star_htseq_edgeR_y)
ucsc_trim_star_htseq_edgeR_keep=filterByExpr(ucsc_trim_star_htseq_edgeR_y)
ucsc_raw_star_feature_edgeR_y <- ucsc_raw_star_feature_edgeR_y[ucsc_raw_star_feature_edgeR_keep,, keep.lib.sizes=FALSE]
ucsc_trim_star_feature_edgeR_y <- ucsc_trim_star_feature_edgeR_y[ucsc_trim_star_feature_edgeR_keep,, keep.lib.sizes=FALSE]
ucsc_raw_star_htseq_edgeR_y <- ucsc_raw_star_htseq_edgeR_y[ucsc_raw_star_htseq_edgeR_keep,, keep.lib.sizes=FALSE]
ucsc_trim_star_htseq_edgeR_y <- ucsc_trim_star_htseq_edgeR_y[ucsc_trim_star_htseq_edgeR_keep,, keep.lib.sizes=FALSE]
ucsc_raw_star_feature_edgeR_y <- calcNormFactors(ucsc_raw_star_feature_edgeR_y)
ucsc_trim_star_feature_edgeR_y <- calcNormFactors(ucsc_trim_star_feature_edgeR_y)
ucsc_raw_star_htseq_edgeR_y <- calcNormFactors(ucsc_raw_star_htseq_edgeR_y)
ucsc_trim_star_htseq_edgeR_y <- calcNormFactors(ucsc_trim_star_htseq_edgeR_y)
ucsc_raw_star_feature_edgeR_y <- estimateDisp(ucsc_raw_star_feature_edgeR_y, edgeR_design)
ucsc_trim_star_feature_edgeR_y <- estimateDisp(ucsc_trim_star_feature_edgeR_y, edgeR_design)
ucsc_raw_star_htseq_edgeR_y <- estimateDisp(ucsc_raw_star_htseq_edgeR_y, edgeR_design)
ucsc_trim_star_htseq_edgeR_y <- estimateDisp(ucsc_trim_star_htseq_edgeR_y, edgeR_design)
ucsc_raw_star_feature_edgeR_fit <- glmQLFit(ucsc_raw_star_feature_edgeR_y, edgeR_design)
ucsc_trim_star_feature_edgeR_fit <- glmQLFit(ucsc_trim_star_feature_edgeR_y, edgeR_design)
ucsc_raw_star_htseq_edgeR_fit <- glmQLFit(ucsc_raw_star_htseq_edgeR_y, edgeR_design)
ucsc_trim_star_htseq_edgeR_fit <- glmQLFit(ucsc_trim_star_htseq_edgeR_y, edgeR_design)
ucsc_raw_star_feature_edgeR_qlf <- glmQLFTest(ucsc_raw_star_feature_edgeR_fit, coef="conditiontrt")
ucsc_trim_star_feature_edgeR_qlf <- glmQLFTest(ucsc_trim_star_feature_edgeR_fit, coef="conditiontrt")
ucsc_raw_star_htseq_edgeR_qlf <- glmQLFTest(ucsc_raw_star_htseq_edgeR_fit, coef="conditiontrt")
ucsc_trim_star_htseq_edgeR_qlf <- glmQLFTest(ucsc_trim_star_htseq_edgeR_fit, coef="conditiontrt")


# DESeq Pseudomapped (salmon)
ucsc_raw_salmon_dds <- DESeqDataSetFromTximport( ucsc_raw_salmon_txi, colData = salmon_samples , design = ~condition)
ucsc_trim_salmon_dds <- DESeqDataSetFromTximport( ucsc_trim_salmon_txi, colData = salmon_samples, design = ~condition)
ucsc_raw_salmon_dds$condition <- relevel(ucsc_raw_salmon_dds$condition, ref = "ctrl")
ucsc_trim_salmon_dds$condition <- relevel(ucsc_trim_salmon_dds$condition, ref = "ctrl")
ucsc_raw_salmon_dds_keep <- rowSums(counts(ucsc_raw_salmon_dds) >= 5) >= 3
ucsc_trim_salmon_dds_keep <- rowSums(counts(ucsc_trim_salmon_dds) >= 5) >= 3
ucsc_raw_salmon_dds <- ucsc_raw_salmon_dds[ucsc_raw_salmon_dds_keep,]
ucsc_trim_salmon_dds <- ucsc_trim_salmon_dds[ucsc_trim_salmon_dds_keep,]
ucsc_raw_salmon_dds <- DESeq(ucsc_raw_salmon_dds)
ucsc_trim_salmon_dds <- DESeq(ucsc_trim_salmon_dds)


# edgeR pseudomapped data
ucsc_raw_salmon_txi_cts <- ucsc_raw_salmon_txi$counts
ucsc_raw_salmon_txi_normMat <- ucsc_raw_salmon_txi$length 
ucsc_raw_salmon_txi_normMat <- ucsc_raw_salmon_txi_normMat/exp(rowMeans(log(ucsc_raw_salmon_txi_normMat)))
ucsc_raw_salmon_txi_normCts <- ucsc_raw_salmon_txi_cts/ucsc_raw_salmon_txi_normMat
ucsc_raw_salmon_txi_eff.lib <- calcNormFactors(ucsc_raw_salmon_txi_normCts)*colSums(ucsc_raw_salmon_txi_normCts)
ucsc_raw_salmon_txi_normMat <- sweep(ucsc_raw_salmon_txi_normMat, 2, ucsc_raw_salmon_txi_eff.lib, "*")
ucsc_raw_salmon_txi_normMat <- log(ucsc_raw_salmon_txi_normMat)
ucsc_raw_salmon_edgeR_y <- DGEList(ucsc_raw_salmon_txi_cts)
ucsc_raw_salmon_edgeR_y <- scaleOffset(ucsc_raw_salmon_edgeR_y, ucsc_raw_salmon_txi_normMat)
ucsc_raw_salmon_edgeR_y_keep <- filterByExpr(ucsc_raw_salmon_edgeR_y, edgeR_design)
ucsc_raw_salmon_edgeR_y <- ucsc_raw_salmon_edgeR_y[ucsc_raw_salmon_edgeR_y_keep,]
ucsc_raw_salmon_edgeR_y <- calcNormFactors(ucsc_raw_salmon_edgeR_y)
ucsc_raw_salmon_edgeR_y <- estimateDisp(ucsc_raw_salmon_edgeR_y, edgeR_design)
ucsc_raw_salmon_edgeR_fit <- glmQLFit(ucsc_raw_salmon_edgeR_y, edgeR_design)
ucsc_raw_salmon_edgeR_qlf <- glmQLFTest(ucsc_raw_salmon_edgeR_fit, coef = "conditiontrt")

ucsc_trim_salmon_txi_cts <- ucsc_trim_salmon_txi$counts
ucsc_trim_salmon_txi_normMat <- ucsc_trim_salmon_txi$length 
ucsc_trim_salmon_txi_normMat <- ucsc_trim_salmon_txi_normMat/exp(rowMeans(log(ucsc_trim_salmon_txi_normMat)))
ucsc_trim_salmon_txi_normCts <- ucsc_trim_salmon_txi_cts/ucsc_trim_salmon_txi_normMat
ucsc_trim_salmon_txi_eff.lib <- calcNormFactors(ucsc_trim_salmon_txi_normCts)*colSums(ucsc_trim_salmon_txi_normCts)
ucsc_trim_salmon_txi_normMat <- sweep(ucsc_trim_salmon_txi_normMat, 2, ucsc_trim_salmon_txi_eff.lib, "*")
ucsc_trim_salmon_txi_normMat <- log(ucsc_trim_salmon_txi_normMat)
ucsc_trim_salmon_edgeR_y <- DGEList(ucsc_trim_salmon_txi_cts)
ucsc_trim_salmon_edgeR_y <- scaleOffset(ucsc_trim_salmon_edgeR_y, ucsc_trim_salmon_txi_normMat)
ucsc_trim_salmon_edgeR_y_keep <- filterByExpr(ucsc_trim_salmon_edgeR_y, edgeR_design)
ucsc_trim_salmon_edgeR_y <- ucsc_trim_salmon_edgeR_y[ucsc_trim_salmon_edgeR_y_keep,]
ucsc_trim_salmon_edgeR_y <- calcNormFactors(ucsc_trim_salmon_edgeR_y)
ucsc_trim_salmon_edgeR_y <- estimateDisp(ucsc_trim_salmon_edgeR_y, edgeR_design)
ucsc_trim_salmon_edgeR_fit <- glmQLFit(ucsc_trim_salmon_edgeR_y, edgeR_design)
ucsc_trim_salmon_edgeR_qlf <- glmQLFTest(ucsc_trim_salmon_edgeR_fit, coef = "conditiontrt")







## -------------RESULTS ************** RESULTS -----------

#----DESeq mapped data

ucsc_raw_star_feature_DESeq_Results <- results(ucsc_raw_star_feature_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ucsc_raw_star_feature_DESeq_ResSig <- subset(ucsc_raw_star_feature_DESeq_Results, padj < 0.1)
ucsc_raw_star_feature_DESeq_ResSig_DF <- as.data.frame(ucsc_raw_star_feature_DESeq_ResSig)
#sink("./results/june/ucsc_raw_star_feature_DESeq_Results_summary.txt")
#summary(ucsc_raw_star_feature_DESeq_Results)
#sink()
#write.csv(ucsc_raw_star_feature_DESeq_ResSig_DF, "./results/june/ucsc_raw_star_feature_DESeq_ResSig_DF.csv", row.names = TRUE)

ucsc_trim_star_feature_DESeq_Results <- results(ucsc_trim_star_feature_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ucsc_trim_star_feature_DESeq_ResSig <- subset(ucsc_trim_star_feature_DESeq_Results, padj < 0.1)
ucsc_trim_star_feature_DESeq_ResSig_DF <- as.data.frame(ucsc_trim_star_feature_DESeq_ResSig)
#sink("./results/june/ucsc_trim_star_feature_DESeq_Results_summary.txt")
#summary(ucsc_trim_star_feature_DESeq_Results)
#sink()
#write.csv(ucsc_trim_star_feature_DESeq_ResSig_DF, "./results/june/ucsc_trim_star_feature_DESeq_ResSig_DF.csv", row.names = TRUE)


ucsc_raw_star_htseq_DESeq_Results <- results(ucsc_raw_star_htseq_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ucsc_raw_star_htseq_DESeq_ResSig <- subset(ucsc_raw_star_htseq_DESeq_Results, padj < 0.1)
ucsc_raw_star_htseq_DESeq_ResSig_DF <- as.data.frame(ucsc_raw_star_htseq_DESeq_ResSig)
#sink("./results/june/ucsc_raw_star_htseq_DESeq_Results_summary.txt")
#summary(ucsc_raw_star_htseq_DESeq_Results)
#sink()
#write.csv(ucsc_raw_star_htseq_DESeq_ResSig_DF, "./results/june/ucsc_raw_star_htseq_DESeq_ResSig_DF.csv", row.names = TRUE)



ucsc_trim_star_htseq_DESeq_Results <- results(ucsc_trim_star_htseq_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ucsc_trim_star_htseq_DESeq_ResSig <- subset(ucsc_trim_star_htseq_DESeq_Results, padj < 0.1)
ucsc_trim_star_htseq_DESeq_ResSig_DF <- as.data.frame(ucsc_trim_star_htseq_DESeq_ResSig)
#sink("./results/june/ucsc_trim_star_htseq_DESeq_Results_summary.txt")
#summary(ucsc_trim_star_htseq_DESeq_Results)
#sink()
#write.csv(ucsc_trim_star_htseq_DESeq_ResSig_DF, "./results/june/ucsc_trim_star_htseq_DESeq_ResSig_DF.csv", row.names = TRUE)

# DESeq2 pseudomapped data

ucsc_raw_salmon_DESeq_Results <- results(ucsc_raw_salmon_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ucsc_raw_salmon_DESeq_ResSig <- subset(ucsc_raw_salmon_DESeq_Results, padj < 0.1)
ucsc_raw_salmon_DESeq_ResSig_DF <- as.data.frame(ucsc_raw_salmon_DESeq_ResSig)
#sink("./results/june/ucsc_raw_salmon_DESeq_Results_summary.txt")
#summary(ucsc_raw_salmon_DESeq_Results)
#sink()
#write.csv(ucsc_raw_salmon_DESeq_ResSig_DF, "./results/june/ucsc_raw_salmon_DESeq_ResSig_DF.csv", row.names = TRUE)



ucsc_trim_salmon_DESeq_Results <- results(ucsc_trim_salmon_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
ucsc_trim_salmon_DESeq_ResSig <- subset(ucsc_trim_salmon_DESeq_Results, padj < 0.1)
ucsc_trim_salmon_DESeq_ResSig_DF <- as.data.frame(ucsc_trim_salmon_DESeq_ResSig)
#sink("./results/june/ucsc_trim_salmon_DESeq_Results_summary.txt")
#summary(ucsc_trim_salmon_DESeq_Results)
#sink()
#write.csv(ucsc_trim_salmon_DESeq_ResSig_DF, "./results/june/ucsc_trim_salmon_DESeq_ResSig_DF.csv", row.names = TRUE)



#----------RESULTS edgeR mapped data

ucsc_raw_star_feature_edgeR_topTags <- topTags(ucsc_raw_star_feature_edgeR_qlf, p.value = 0.1)
ucsc_raw_star_feature_edgeR_topTags_DF <- as.data.frame(ucsc_raw_star_feature_edgeR_topTags)
#sink("./results/june/ucsc_raw_star_feature_edgeR_qlf_summary.txt")
#summary(decideTests(ucsc_raw_star_feature_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ucsc_raw_star_feature_edgeR_topTags_DF, "./results/june/ucsc_raw_star_feature_edgeR_topTags_DF.csv", row.names = TRUE)



ucsc_trim_star_feature_edgeR_topTags <- topTags(ucsc_trim_star_feature_edgeR_qlf, p.value = 0.1)
ucsc_trim_star_feature_edgeR_topTags_DF <- as.data.frame(ucsc_trim_star_feature_edgeR_topTags)
#sink("./results/june/ucsc_trim_star_feature_edgeR_qlf_summary.txt")
#summary(decideTests(ucsc_trim_star_feature_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ucsc_trim_star_feature_edgeR_topTags_DF, "./results/june/ucsc_trim_star_feature_edgeR_topTags_DF.csv", row.names = TRUE)


ucsc_raw_star_htseq_edgeR_topTags <- topTags(ucsc_raw_star_htseq_edgeR_qlf, p.value = 0.1)
ucsc_raw_star_htseq_edgeR_topTags_DF <- as.data.frame(ucsc_raw_star_htseq_edgeR_topTags)
sink("./results/june/ucsc_raw_star_htseq_edgeR_qlf_summary.txt")
#summary(decideTests(ucsc_raw_star_htseq_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ucsc_raw_star_htseq_edgeR_topTags_DF, "./results/june/ucsc_raw_star_htseq_edgeR_topTags_DF.csv", row.names = TRUE)


ucsc_trim_star_htseq_edgeR_topTags <- topTags(ucsc_trim_star_htseq_edgeR_qlf, p.value = 0.1)
ucsc_trim_star_htseq_edgeR_topTags_DF <- as.data.frame(ucsc_trim_star_htseq_edgeR_topTags)
#sink("./results/june/ucsc_trim_star_htseq_edgeR_qlf_summary.txt")
#summary(decideTests(ucsc_trim_star_htseq_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ucsc_trim_star_htseq_edgeR_topTags_DF, "./results/june/ucsc_trim_star_htseq_edgeR_topTags_DF.csv", row.names = TRUE)

# edgeR pseudomapped data

ucsc_raw_salmon_edgeR_topTags <- topTags(ucsc_raw_salmon_edgeR_qlf, p.value = 0.1)
ucsc_raw_salmon_edgeR_topTags_DF <- as.data.frame(ucsc_raw_salmon_edgeR_topTags)
#sink("./results/june/ucsc_raw_salmon_edgeR_qlf_summary.txt")
#summary(decideTests(ucsc_raw_salmon_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ucsc_raw_salmon_edgeR_topTags_DF, "./results/june/ucsc_raw_salmon_edgeR_topTags_DF.csv", row.names = TRUE)


ucsc_trim_salmon_edgeR_topTags <- topTags(ucsc_trim_salmon_edgeR_qlf, p.value = 0.1)
ucsc_trim_salmon_edgeR_topTags_DF <- as.data.frame(ucsc_trim_salmon_edgeR_topTags)
#sink("./results/june/ucsc_trim_salmon_edgeR_qlf_summary.txt")
#summary(decideTests(ucsc_trim_salmon_edgeR_qlf, p.value = 0.1))
#sink()
#write.csv(ucsc_trim_salmon_edgeR_topTags_DF, "./results/june/ucsc_trim_salmon_edgeR_topTags_DF.csv", row.names = TRUE)


###---------------------------------------


ucsc_trim_salmon_edgeR_topTags_DF







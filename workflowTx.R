# import all necessary libraries
library("DESeq2")
library("readr")
library("tximport")
library("edgeR")
library("httpgd")
library("tximport")

#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn")
# or just via cran install.packages("ggvenn")
library("ggvenn")
library("ggpubr")
library("vsn")


# set up graphics
hgd()



# GENERAL conditions
condition <- factor(c(rep("ctrl", 5), rep("trt", 5), rep("sham", 3)))
edgeR_design <- model.matrix(~condition)

# GENERAL (sample&condition table for salmon)
salmon_samples <- read.table("../samples_april.csv", sep = ",", header = TRUE, row.names = 1)
salmon_samples$condition <- as.factor(salmon_samples$condition)

# GENERAL import of SALMON mapped to ensembl
ensembl_tx2gene <- read_csv("genome/ensembl/tx2gene_GRCm39.109.csv", col_names = FALSE)
ensembl_raw_salmon_files <- file.path("mapping/ensembl/raw/salmon", rownames(salmon_samples), "quant.sf")
ensembl_trim_salmon_files <- file.path("mapping/ensembl/trimmomatic/salmon", rownames(salmon_samples), "quant.sf")
names(ensembl_raw_salmon_files) <- rownames(salmon_samples)
names(ensembl_trim_salmon_files) <- rownames(salmon_samples)
tx_ensembl_raw_salmon_txi <- tximport(ensembl_raw_salmon_files, type = "salmon", txIn = TRUE, txOut = TRUE)
tx_ensembl_trim_salmon_txi <- tximport(ensembl_trim_salmon_files, type = "salmon", txIn = TRUE, txOut = TRUE)

# DESeq2 model, ensembl
tx_ensembl_raw_salmon_dds <- DESeqDataSetFromTximport(tx_ensembl_raw_salmon_txi, colData = salmon_samples , design = ~condition)
tx_ensembl_trim_salmon_dds <- DESeqDataSetFromTximport(tx_ensembl_trim_salmon_txi, colData = salmon_samples , design = ~condition)
tx_ensembl_raw_salmon_dds$condition <- relevel(tx_ensembl_raw_salmon_dds$condition, ref = "ctrl")
tx_ensembl_trim_salmon_dds$condition <- relevel(tx_ensembl_trim_salmon_dds$condition, ref = "ctrl")
tx_ensembl_raw_salmon_dds_keep <- rowSums(counts(tx_ensembl_raw_salmon_dds) >= 5) >=3
tx_ensembl_trim_salmon_dds_keep <- rowSums(counts(tx_ensembl_trim_salmon_dds) >= 5) >=3
tx_ensembl_raw_salmon_dds <- tx_ensembl_raw_salmon_dds[tx_ensembl_raw_salmon_dds_keep,]
tx_ensembl_trim_salmon_dds <- tx_ensembl_trim_salmon_dds[tx_ensembl_trim_salmon_dds_keep,]
tx_ensembl_raw_salmon_dds <- DESeq(tx_ensembl_raw_salmon_dds)
tx_ensembl_trim_salmon_dds <- DESeq(tx_ensembl_trim_salmon_dds)

# edgeR Model, ensembl
tx_ensembl_raw_salmon_txi_cts <- tx_ensembl_raw_salmon_txi$counts
tx_ensembl_raw_salmon_txi_normMat <- tx_ensembl_raw_salmon_txi$length 
tx_ensembl_raw_salmon_txi_normMat <- tx_ensembl_raw_salmon_txi_normMat/exp(rowMeans(log(tx_ensembl_raw_salmon_txi_normMat)))
tx_ensembl_raw_salmon_txi_normCts <- tx_ensembl_raw_salmon_txi_cts/tx_ensembl_raw_salmon_txi_normMat
tx_ensembl_raw_salmon_txi_eff.lib <- calcNormFactors(tx_ensembl_raw_salmon_txi_normCts)*colSums(tx_ensembl_raw_salmon_txi_normCts)
tx_ensembl_raw_salmon_txi_normMat <- sweep(tx_ensembl_raw_salmon_txi_normMat, 2, tx_ensembl_raw_salmon_txi_eff.lib, "*")
tx_ensembl_raw_salmon_txi_normMat <- log(tx_ensembl_raw_salmon_txi_normMat)
tx_ensembl_raw_salmon_edgeR_y <- DGEList(tx_ensembl_raw_salmon_txi_cts)
tx_ensembl_raw_salmon_edgeR_y <- scaleOffset(tx_ensembl_raw_salmon_edgeR_y, tx_ensembl_raw_salmon_txi_normMat)
tx_ensembl_raw_salmon_edgeR_y_keep <- filterByExpr(tx_ensembl_raw_salmon_edgeR_y, edgeR_design)
tx_ensembl_raw_salmon_edgeR_y <- tx_ensembl_raw_salmon_edgeR_y[tx_ensembl_raw_salmon_edgeR_y_keep,]
tx_ensembl_raw_salmon_edgeR_y <- calcNormFactors(tx_ensembl_raw_salmon_edgeR_y)
tx_ensembl_raw_salmon_edgeR_y <- estimateDisp(tx_ensembl_raw_salmon_edgeR_y, edgeR_design)
tx_ensembl_raw_salmon_edgeR_fit <- glmQLFit(tx_ensembl_raw_salmon_edgeR_y, edgeR_design)
tx_ensembl_raw_salmon_edgeR_qlf <- glmQLFTest(tx_ensembl_raw_salmon_edgeR_fit, coef = "conditiontrt")

tx_ensembl_trim_salmon_txi_cts <- tx_ensembl_trim_salmon_txi$counts
tx_ensembl_trim_salmon_txi_normMat <- tx_ensembl_trim_salmon_txi$length 
tx_ensembl_trim_salmon_txi_normMat <- tx_ensembl_trim_salmon_txi_normMat/exp(rowMeans(log(tx_ensembl_trim_salmon_txi_normMat)))
tx_ensembl_trim_salmon_txi_normCts <- tx_ensembl_trim_salmon_txi_cts/tx_ensembl_trim_salmon_txi_normMat
tx_ensembl_trim_salmon_txi_eff.lib <- calcNormFactors(tx_ensembl_trim_salmon_txi_normCts)*colSums(tx_ensembl_trim_salmon_txi_normCts)
tx_ensembl_trim_salmon_txi_normMat <- sweep(tx_ensembl_trim_salmon_txi_normMat, 2, tx_ensembl_trim_salmon_txi_eff.lib, "*")
tx_ensembl_trim_salmon_txi_normMat <- log(tx_ensembl_trim_salmon_txi_normMat)
tx_ensembl_trim_salmon_edgeR_y <- DGEList(tx_ensembl_trim_salmon_txi_cts)
tx_ensembl_trim_salmon_edgeR_y <- scaleOffset(tx_ensembl_trim_salmon_edgeR_y, tx_ensembl_trim_salmon_txi_normMat)
tx_ensembl_trim_salmon_edgeR_y_keep <- filterByExpr(tx_ensembl_trim_salmon_edgeR_y, edgeR_design)
tx_ensembl_trim_salmon_edgeR_y <- tx_ensembl_trim_salmon_edgeR_y[tx_ensembl_trim_salmon_edgeR_y_keep,]
tx_ensembl_trim_salmon_edgeR_y <- calcNormFactors(tx_ensembl_trim_salmon_edgeR_y)
tx_ensembl_trim_salmon_edgeR_y <- estimateDisp(tx_ensembl_trim_salmon_edgeR_y, edgeR_design)
tx_ensembl_trim_salmon_edgeR_fit <- glmQLFit(tx_ensembl_trim_salmon_edgeR_y, edgeR_design)
tx_ensembl_trim_salmon_edgeR_qlf <- glmQLFTest(tx_ensembl_trim_salmon_edgeR_fit, coef = "conditiontrt")



## -------------RESULTS ************** RESULTS -----------

# DESeq2 pseudomapped data
tx_ensembl_raw_salmon_DESeq_Results <- results(tx_ensembl_raw_salmon_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
tx_ensembl_raw_salmon_DESeq_ResSig <- subset(tx_ensembl_raw_salmon_DESeq_Results, padj < 0.1)
tx_ensembl_raw_salmon_DESeq_ResSig_DF <- as.data.frame(tx_ensembl_raw_salmon_DESeq_ResSig)
sink("./results/june/sum/tx_ensembl_raw_salmon_DESeq_Results_summary.txt")
summary(tx_ensembl_raw_salmon_DESeq_Results)
sink()
write.csv(tx_ensembl_raw_salmon_DESeq_ResSig_DF, "./results/june/tables/tx_ensembl_raw_salmon_DESeq_ResSig_DF.csv", row.names = TRUE)

summary(subset(tx_ensembl_raw_salmon_DESeq_Results, abs(log2FoldChange) > 1))
summary(tx_ensembl_raw_salmon_DESeq_Results)

tx_ensembl_trim_salmon_DESeq_Results <- results(tx_ensembl_trim_salmon_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
tx_ensembl_trim_salmon_DESeq_ResSig <- subset(tx_ensembl_trim_salmon_DESeq_Results, padj < 0.1)
tx_ensembl_trim_salmon_DESeq_ResSig_DF <- as.data.frame(tx_ensembl_trim_salmon_DESeq_ResSig)
sink("./results/june/sum/tx_ensembl_trim_salmon_DESeq_Results_summary.txt")
summary(tx_ensembl_trim_salmon_DESeq_Results)
sink()
write.csv(tx_ensembl_trim_salmon_DESeq_ResSig_DF, "./results/june/tables/tx_ensembl_trim_salmon_DESeq_ResSig_DF.csv", row.names = TRUE)


# edgeR pseudomapped data
tx_ensembl_raw_salmon_edgeR_topTags <- topTags(tx_ensembl_raw_salmon_edgeR_qlf, p.value = 0.1)
tx_ensembl_raw_salmon_edgeR_topTags_DF <- as.data.frame(tx_ensembl_raw_salmon_edgeR_topTags)
sink("./results/june/sum/tx_ensembl_raw_salmon_edgeR_qlf_summary.txt")
summary(decideTests(tx_ensembl_raw_salmon_edgeR_qlf, p.value = 0.1))
sink()
write.csv(tx_ensembl_raw_salmon_edgeR_topTags_DF, "./results/june/tables/tx_ensembl_raw_salmon_edgeR_topTags_DF.csv", row.names = TRUE)


tx_ensembl_trim_salmon_edgeR_topTags <- topTags(tx_ensembl_trim_salmon_edgeR_qlf, p.value = 0.1)
tx_ensembl_trim_salmon_edgeR_topTags_DF <- as.data.frame(tx_ensembl_trim_salmon_edgeR_topTags)
sink("./results/june/sum/tx_ensembl_trim_salmon_edgeR_qlf_summary.txt")
summary(decideTests(tx_ensembl_trim_salmon_edgeR_qlf, p.value = 0.1))
sink()
write.csv(tx_ensembl_trim_salmon_edgeR_topTags_DF, "./results/june/tables/tx_ensembl_trim_salmon_edgeR_topTags_DF.csv", row.names = TRUE)







####------------------ VISUALIZATION ***** VISUALIZATION ---------------------

# PLOTS compare to rafaels results
rafael_tx <- read.table("./results/rafael/rafael_Tx.csv", sep = ",", header = TRUE, row.names = 1)

significant_tx <- list(
    "rafael" = row.names(rafael_tx),
    "tx_ensembl_raw_salmon_DESeq" = row.names(tx_ensembl_raw_salmon_DESeq_ResSig_DF),
    "tx_ensembl_trim_salmon_DESeq" = row.names(tx_ensembl_trim_salmon_DESeq_ResSig_DF),
    "tx_ensembl_raw_salmon_edgeR" = row.names(tx_ensembl_raw_salmon_edgeR_topTags_DF), 
    "tx_ensembl_trim_salmon_edgeR" = row.names(tx_ensembl_trim_salmon_edgeR_topTags_DF)
)

jpeg(file="./results/june/plots/tx_ensembl_salmon.jpeg", width = 1500, height = 1500)
elements <- ggvenn(significant_tx, c("rafael", "tx_ensembl_raw_salmon_DESeq", "tx_ensembl_trim_salmon_DESeq"), 
    show_elements = TRUE, label_sep = "\n", show_percentage = TRUE, digits = 2)
percentage <- ggvenn(significant_tx, c("rafael", "tx_ensembl_raw_salmon_DESeq", "tx_ensembl_trim_salmon_DESeq"), 
    show_percentage = TRUE, digits = 0, text_size = 6)
ggarrange(percentage, elements, labels = c("Percent", "Elements"),
ncol = 2, nrow = 1)
dev.off()

















#----IF I USE THE STAR ALIGNMENT TO TRANSCRIPTS

# DESeq2 
tx_ensembl_raw_star_feature_DESeq_Results <- results(tx_ensembl_raw_star_feature_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
tx_ensembl_raw_star_feature_DESeq_ResSig <- subset(tx_ensembl_raw_star_feature_DESeq_Results, padj < 0.1)
tx_ensembl_raw_star_feature_DESeq_ResSig_DF <- as.data.frame(tx_ensembl_raw_star_feature_DESeq_ResSig)
sink("./results/june/sum/tx_ensembl_raw_star_feature_DESeq_Results_summary.txt")
summary(tx_ensembl_raw_star_feature_DESeq_Results)
sink()
write.csv(tx_ensembl_raw_star_feature_DESeq_ResSig_DF, "./results/june/tables/tx_ensembl_raw_star_feature_DESeq_ResSig_DF.csv", row.names = TRUE)

tx_ensembl_trim_star_feature_DESeq_Results <- results(tx_ensembl_trim_star_feature_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
tx_ensembl_trim_star_feature_DESeq_ResSig <- subset(tx_ensembl_trim_star_feature_DESeq_Results, padj < 0.1)
tx_ensembl_trim_star_feature_DESeq_ResSig_DF <- as.data.frame(tx_ensembl_trim_star_feature_DESeq_ResSig)
sink("./results/june/sum/tx_ensembl_trim_star_feature_DESeq_Results_summary.txt")
summary(tx_ensembl_trim_star_feature_DESeq_Results)
sink()
write.csv(tx_ensembl_trim_star_feature_DESeq_ResSig_DF, "./results/june/tables/tx_ensembl_trim_star_feature_DESeq_ResSig_DF.csv", row.names = TRUE)


tx_ensembl_raw_star_htseq_DESeq_Results <- results(tx_ensembl_raw_star_htseq_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
tx_ensembl_raw_star_htseq_DESeq_ResSig <- subset(tx_ensembl_raw_star_htseq_DESeq_Results, padj < 0.1)
tx_ensembl_raw_star_htseq_DESeq_ResSig_DF <- as.data.frame(tx_ensembl_raw_star_htseq_DESeq_ResSig)
sink("./results/june/sum/tx_ensembl_raw_star_htseq_DESeq_Results_summary.txt")
summary(tx_ensembl_raw_star_htseq_DESeq_Results)
sink()
write.csv(tx_ensembl_raw_star_htseq_DESeq_ResSig_DF, "./results/june/tables/tx_ensembl_raw_star_htseq_DESeq_ResSig_DF.csv", row.names = TRUE)






tx_ensembl_trim_star_htseq_DESeq_Results <- results(tx_ensembl_trim_star_htseq_dds, alpha = 0.1, contrast = c("condition", "trt", "ctrl"))
tx_ensembl_trim_star_htseq_DESeq_ResSig <- subset(tx_ensembl_trim_star_htseq_DESeq_Results, padj < 0.1)
tx_ensembl_trim_star_htseq_DESeq_ResSig_DF <- as.data.frame(tx_ensembl_trim_star_htseq_DESeq_ResSig)
sink("./results/june/sum/tx_ensembl_trim_star_htseq_DESeq_Results_summary.txt")
summary(tx_ensembl_trim_star_htseq_DESeq_Results)
sink()
write.csv(tx_ensembl_trim_star_htseq_DESeq_ResSig_DF, "./results/june/tables/tx_ensembl_trim_star_htseq_DESeq_ResSig_DF.csv", row.names = TRUE)

#----------RESULTS edgeR mapped data

tx_ensembl_raw_star_feature_edgeR_topTags <- topTags(tx_ensembl_raw_star_feature_edgeR_qlf, p.value = 0.1)
tx_ensembl_raw_star_feature_edgeR_topTags_DF <- as.data.frame(tx_ensembl_raw_star_feature_edgeR_topTags)
sink("./results/june/sum/tx_ensembl_raw_star_feature_edgeR_qlf_summary.txt")
summary(decideTests(tx_ensembl_raw_star_feature_edgeR_qlf, p.value = 0.1))
sink()
write.csv(tx_ensembl_raw_star_feature_edgeR_topTags_DF, "./results/june/tables/tx_ensembl_raw_star_feature_edgeR_topTags_DF.csv", row.names = TRUE)


tx_ensembl_trim_star_feature_edgeR_topTags <- topTags(tx_ensembl_trim_star_feature_edgeR_qlf, p.value = 0.1)
tx_ensembl_trim_star_feature_edgeR_topTags_DF <- as.data.frame(tx_ensembl_trim_star_feature_edgeR_topTags)
sink("./results/june/sum/tx_ensembl_trim_star_feature_edgeR_qlf_summary.txt")
summary(decideTests(tx_ensembl_trim_star_feature_edgeR_qlf, p.value = 0.1))
sink()
write.csv(tx_ensembl_trim_star_feature_edgeR_topTags_DF, "./results/june/tables/tx_ensembl_trim_star_feature_edgeR_topTags_DF.csv", row.names = TRUE)


tx_ensembl_raw_star_htseq_edgeR_topTags <- topTags(tx_ensembl_raw_star_htseq_edgeR_qlf, p.value = 0.1)
tx_ensembl_raw_star_htseq_edgeR_topTags_DF <- as.data.frame(tx_ensembl_raw_star_htseq_edgeR_topTags)
sink("./results/june/sum/tx_ensembl_raw_star_htseq_edgeR_qlf_summary.txt")
summary(decideTests(tx_ensembl_raw_star_htseq_edgeR_qlf, p.value = 0.1))
sink()
write.csv(tx_ensembl_raw_star_htseq_edgeR_topTags_DF, "./results/june/tables/tx_ensembl_raw_star_htseq_edgeR_topTags_DF.csv", row.names = TRUE)


tx_ensembl_trim_star_htseq_edgeR_topTags <- topTags(tx_ensembl_trim_star_htseq_edgeR_qlf, p.value = 0.1)
tx_ensembl_trim_star_htseq_edgeR_topTags_DF <- as.data.frame(tx_ensembl_trim_star_htseq_edgeR_topTags)
sink("./results/june/sum/tx_ensembl_trim_star_htseq_edgeR_qlf_summary.txt")
summary(decideTests(tx_ensembl_trim_star_htseq_edgeR_qlf, p.value = 0.1))
sink()
write.csv(tx_ensembl_trim_star_htseq_edgeR_topTags_DF, "./results/june/tables/tx_ensembl_trim_star_htseq_edgeR_topTags_DF.csv", row.names = TRUE)





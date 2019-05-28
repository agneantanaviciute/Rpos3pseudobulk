#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2", version = "3.8")
#
#if (!require("ggplot2")){
#  install.packahes("ggplot2")
#}

##required libraries
options(stringsAsFactors = FALSE)
library(DESeq2)
library(ggplot2)

#read in cell meta data
meta.data.HC <- read.delim("meta.data.HC.txt")
meta.data.DSS <- read.delim("meta.data.DSS.txt")

##read in raw UMI counts matrices
readRDS("counts.HC.Rds")
readRDS("counts.DSS.Rds")


##pseudobulk counts per sample for stromal cells only in HC:
##discarding non-fibroblast cell clusters
stromal.cells.HC <- which(meta.data.HC$Cluster %in% c("Str1.1", "Str1.3", 
                                           "Str1.2", "Str2", "Str3",
                                           "Myofibroblast"))
counts.HC <- counts.HC[, stromal.cells.HC]
meta.data.HC <- meta.data.HC[stromal.cells.HC, ]

pseudobulk.HC <- do.call(cbind, lapply(unique(meta.data.HC$Sample), function(x){
  Matrix::rowSums(counts.HC [, meta.data.HC$Sample == x])
}))
colnames(pseudobulk.HC) <- paste0("HC_", unique(meta.data.HC$Sample))


##pseudobulk counts per sample for stromal cells only in DSS:
##discarding non-fiboblast cell clusters
stromal.cells.DSS <- which(meta.data.DSS$Cluster %in% c("Str1.1", "Str1.3", "Str1.2",
                                            "Str2", "Str3.1",  "Str3.2", 
                                            "Str3.3", "Str4", "Myofibroblast"))

counts.DSS <- counts.DSS[, stromal.cells.DSS]
meta.data.DSS <- meta.data.DSS[stromal.cells.DSS, ]

pseudobulk.DSS <- do.call(cbind, lapply(unique(meta.data.DSS$Sample), function(x){
  Matrix::rowSums(counts.DSS [, meta.data.DSS$Sample == x])
}))
colnames(pseudobulk.DSS) <- paste0("DSS_", unique(meta.data.DSS$Sample))



##pseudobulk format for DESeq2:
pseudobulk <- cbind(pseudobulk.HC, pseudobulk.DSS)
pseudobulk.meta.data <- data.frame(Treatment=c("HC", "HC", "HC", "DSS", "DSS", "DSS"))
rownames(pseudobulk.meta.data) <- colnames(pseudobulk)

##DESeq2 normalisation and DE test
dds <- DESeqDataSetFromMatrix(pseudobulk, pseudobulk.meta.data, design= ~Treatment)
dds <- DESeq(dds)
rl <- rlog(dds)
res <- results(dds)

##gather plot data- normalised Rspo3 counts and means per sample:
plot.data <- data.frame(Rspo3=counts(dds, normalized=T)["Rspo3", ], pseudobulk.meta.data)
plot.data$Mean <- c(rep(mean(plot.data$Rspo3[1:3]), 3),
                    rep(mean(plot.data$Rspo3[4:6]), 3))

##fetch and format log fold change and FDR values:
log2.FC <- res["Rspo3", ]$log2FoldChange
log2.FC <- round(log2.FC, digits = 2)

FDR <- res["Rspo3", ]$padj
FDR <- formatC(FDR)

##plot and export:
png("Rspo3.png", res=300, height=2000, width=2000)

ggplot(plot.data, aes(Treatment, Rspo3, fill=Treatment)
       ) + geom_jitter(shape=21, color="black", size=7, height = 0) + geom_crossbar(aes(y=Mean,
        ymin=Mean, ymax=Mean)) + ylim(0, 1200) + theme_light(base_size = 24) + labs(y="Rspo3 Expression",
title=paste0("log2 FC: ", log2.FC,
             "\nFDR: ", FDR))

dev.off()





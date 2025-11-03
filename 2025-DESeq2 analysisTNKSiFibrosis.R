### DESeq2 analysis of selected contrasts 
### data = mRNA sequencing dataset from treated NHLFs (Control, OM-153, IPF-RC, IPF-RC + Nintedanib, and IPF-RC + OM-153) 

# Installing packages (if needed):
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ggplot2")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("genefilter")
BiocManager::install("regionReport")
BiocManager::install("apeglm")
BiocManager::install(c("DT", "sessioninfo"))
BiocManager::install("InteractiveComplexHeatmap")
install.packages("pheatmap")
install.packages("gplots")

# Loading packages: 
library("gplots")
library("RColorBrewer")
library("ggplot2")
library("pheatmap")
library("genefilter")
library("DESeq2")
library("regionReport")
library("apeglm")
library("edgeR")
library(dplyr)
library(InteractiveComplexHeatmap)
library(ggrepel)

# Load and prepare dataset
setwd("~/Desktop/ANA")
data<-read.table("genes.txt",header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=FALSE)
rownames(data) <- make.unique(data[,2])
dataset <- data[,-c(1:2)]
dataset <- as.matrix(dataset)
dim(dataset)
dataset <- round(dataset, 0)

# Setting factors/metadata from file
mfactors<-read.table("Meta.txt", header=TRUE, sep="\t", comment="", quote="", stringsAsFactors=TRUE)
rownames(mfactors) <- mfactors$sample
Meta <- mfactors[,c(2:4)]
Meta$batch <- as.factor(Meta$batch)
Meta$group <- as.factor(Meta$group)
colnames(dataset) <- mfactors$sample

#Filtering data into contrasts for analysis
OvC <- dataset[,c(1:6)]
OvCm <- Meta[c(1:6),]
IvC <- dataset[,c(1:3,7:9)]
IvCm <- Meta[c(1:3,7:9),]
IOvI <- dataset[,c(7:12)]
IOvIm <- Meta[c(7:12),]
IOvO <- dataset[,c(4:6,10:12)]
IOvOm <- Meta[c(4:6,10:12),]
INvI <- dataset[,c(7:9,13:15)]
INvIm <- Meta[c(7:9,13:15),]
BvI <- dataset[,c(7:9,16:18)]
BvIm <- Meta[c(7:9,16:18),]
BvIO <- dataset[,c(10:12,16:18)]
BvIOm <- Meta[c(10:12,16:18),]
BvIN <- dataset[,c(13:18)]
BvINm <- Meta[c(13:18),]



### Contrast1 - OM-153 vs IPF-RC
IOvIm$group <- c(1,1,1,2,2,2)

# Generating a DESeq2 object:
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = IOvI,
  colData = IOvIm,
  design = ~ group + batch)

# DESeq2 analysis
dds <- ddsFullCountTable
as.data.frame(colData(dds)) #Checking factors
dds <- DESeq(dds)
rld <- rlog(dds, blind=FALSE)

pdf(file="DESeq2-PCAIOvI.pdf", width = 14, height = 8)
z <- plotPCA(rld, intgroup=c("sample"), ntop = 500)
z + geom_label(aes(label = name))
dev.off()

# Result
res <- results(dds, list(c("group")))
summary(res)
res <- lfcShrink(dds, coef="group", type="apeglm")
resOrdered <- res[order(res$padj),]
write.table(as.data.frame(resOrdered), file=paste("DESeq-IPF-153vsIPF.txt", sep="."), quote=FALSE, sep="\t")

pdf(file="DESeq2-IPF-153vsIPF.pdf")

DESeq2::plotMA(res, ylim = (c(-3,3)))

##Result vulcanoplot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res$log2FoldChange) > 1.5 & res$padj < alpha 
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.6)

plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()


# Bach ajustment and data export
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mm <- model.matrix(~group, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
colData(vsd)

write.table(as.data.frame(mat), file=paste("DESeq-normalized-LIMMA-IOvI.txt", sep="."), quote=FALSE, sep="\t", col.names = NA, row.names = TRUE)




### Contrast2 - OM-153 vs Control
IvC <- dataset[,c(1:3,7:9)]
IvCm <- Meta[c(1:3,7:9),]
IvCm$group <- c(1,1,1,2,2,2)


# Generating a DESeq2 object:
ddsFullCountTableIvC <- DESeqDataSetFromMatrix(
  countData = IvC,
  colData = IvCm,
  design = ~ group + batch)

# DESeq2 analysis
ddsIvC <- ddsFullCountTableIvC
as.data.frame(colData(ddsIvC)) #Checking factors
ddsIvC <- DESeq(ddsIvC)
rldIvC <- rlog(ddsIvC, blind=FALSE)

pdf(file="DESeq2-PCAIvC.pdf", width = 14, height = 8)
z <- plotPCA(rldIvC, intgroup=c("sample"), ntop = 500)
z + geom_label(aes(label = name))
dev.off()

resultsNames(ddsIvC)

# Result
resIvC <- results(ddsIvC, list(c("group")))
summary(resIvC)
resIvC <- lfcShrink(ddsIvC, coef="group", type="apeglm")
resOrderedIvC <- resIvC[order(resIvC$padj),]
write.table(as.data.frame(resOrderedIvC), file=paste("DESeq-IvC.txt", sep="."), quote=FALSE, sep="\t")

pdf(file="DESeq2-IvC.pdf")

DESeq2::plotMA(resIvC, ylim = (c(-3,3)))

##Result vulcanoplot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(resIvC$log2FoldChange, -log10(resIvC$pvalue))
plot(resIvC$log2FoldChange, -log10(resIvC$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(resIvC$log2FoldChange) > 1.5 & resIvC$padj < alpha 
text(resIvC$log2FoldChange[gn.selected],
     -log10(resIvC$padj)[gn.selected],
     lab=rownames(resIvC)[gn.selected ], cex=0.6)

plot(resIvC$log2FoldChange, -log10(resIvC$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()


# Bach ajustment and data export
vsdIvC <- vst(ddsIvC, blind=FALSE)
mat <- assay(vsdIvC)
mm <- model.matrix(~group, colData(vsdIvC))
mat <- limma::removeBatchEffect(mat, batch=vsdIvC$batch, design=mm)
assay(vsdIvC) <- mat
colData(vsdIvC)

write.table(as.data.frame(mat), file=paste("DESeq-normalized-LIMMA-IvC.txt", sep="."), quote=FALSE, sep="\t", col.names = NA, row.names = TRUE)




### Contrast3
IOvO <- dataset[,c(4:6,10:12)]
IOvOm <- Meta[c(4:6,10:12),]
IOvOm$group <- c(1,1,1,2,2,2)


# Generating a DESeq2 object:
ddsFullCountTableIOvO <- DESeqDataSetFromMatrix(
  countData = IOvO,
  colData = IOvOm,
  design = ~ group + batch)

# DESeq2 analysis
ddsIOvO <- ddsFullCountTableIOvO
as.data.frame(colData(ddsIOvO)) #Checking factors
ddsIOvO <- DESeq(ddsIOvO)
rldIOvO <- rlog(ddsIOvO, blind=FALSE)

pdf(file="DESeq2-PCAIOvO.pdf", width = 14, height = 8)
z <- plotPCA(rldIOvO, intgroup=c("sample"), ntop = 500)
z + geom_label(aes(label = name))
dev.off()

resultsNames(ddsIOvO)

# Result
resIOvO <- results(ddsIOvO, list(c("group")))
summary(resIOvO)
resIOvO <- lfcShrink(ddsIOvO, coef="group", type="apeglm")
resOrderedIOvO <- resIOvO[order(resIOvO$padj),]
write.table(as.data.frame(resOrderedIOvO), file=paste("DESeq-IOvO.txt", sep="."), quote=FALSE, sep="\t")

pdf(file="DESeq2-IOvO.pdf")

DESeq2::plotMA(resIOvO, ylim = (c(-3,3)))

##Result vulcanoplot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(resIOvO$log2FoldChange, -log10(resIOvO$pvalue))
plot(resIOvO$log2FoldChange, -log10(resIOvO$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(resIOvO$log2FoldChange) > 1.5 & resIOvO$padj < alpha 
text(resIOvO$log2FoldChange[gn.selected],
     -log10(resIOvO$padj)[gn.selected],
     lab=rownames(resIOvO)[gn.selected ], cex=0.6)

plot(resIOvO$log2FoldChange, -log10(resIOvO$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()


# Bach ajustment and data export
vsdIOvO <- vst(ddsIOvO, blind=FALSE)
mat <- assay(vsdIOvO)
mm <- model.matrix(~group, colData(vsdIOvO))
mat <- limma::removeBatchEffect(mat, batch=vsdIOvO$batch, design=mm)
assay(vsdIOvO) <- mat
colData(vsdIOvO)

write.table(as.data.frame(mat), file=paste("DESeq-normalized-LIMMA-IOvO.txt", sep="."), quote=FALSE, sep="\t", col.names = NA, row.names = TRUE)



### Contrast4
INvI <- dataset[,c(7:9,13:15)]
INvIm <- Meta[c(7:9,13:15),]
INvIm$group <- c(1,1,1,2,2,2)


# Generating a DESeq2 object:
ddsFullCountTableINvI <- DESeqDataSetFromMatrix(
  countData = INvI,
  colData = INvIm,
  design = ~ group + batch)

# DESeq2 analysis
ddsINvI <- ddsFullCountTableINvI
as.data.frame(colData(ddsINvI)) #Checking factors
ddsINvI <- DESeq(ddsINvI)
rldINvI <- rlog(ddsINvI, blind=FALSE)

pdf(file="DESeq2-PCAINvI.pdf", width = 14, height = 8)
z <- plotPCA(rldINvI, intgroup=c("sample"), ntop = 500)
z + geom_label(aes(label = name))
dev.off()


# Result
resINvI <- results(ddsINvI, list(c("group")))
summary(resINvI)
resINvI <- lfcShrink(ddsINvI, coef="group", type="apeglm")
resOrderedINvI <- resINvI[order(resINvI$padj),]
write.table(as.data.frame(resOrderedINvI), file=paste("DESeq-INvI.txt", sep="."), quote=FALSE, sep="\t")

pdf(file="DESeq2-INvI.pdf")

DESeq2::plotMA(resINvI, ylim = (c(-3,3)))

##Result vulcanoplot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(resINvI$log2FoldChange, -log10(resINvI$pvalue))
plot(resINvI$log2FoldChange, -log10(resINvI$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(resINvI$log2FoldChange) > 1.5 & resINvI$padj < alpha 
text(resINvI$log2FoldChange[gn.selected],
     -log10(resINvI$padj)[gn.selected],
     lab=rownames(resINvI)[gn.selected ], cex=0.6)

plot(resINvI$log2FoldChange, -log10(resINvI$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()


# Bach ajustment and data export
vsdINvI <- vst(ddsINvI, blind=FALSE)
mat <- assay(vsdINvI)
mm <- model.matrix(~group, colData(vsdINvI))
mat <- limma::removeBatchEffect(mat, batch=vsdINvI$batch, design=mm)
assay(vsdINvI) <- mat
colData(vsdINvI)

write.table(as.data.frame(mat), file=paste("DESeq-normalized-LIMMA-INvI.txt", sep="."), quote=FALSE, sep="\t", col.names = NA, row.names = TRUE)




### Contrast5
BvI <- dataset[,c(7:9,16:18)]
BvIm <- Meta[c(7:9,16:18),]
BvIm$group <- c(1,1,1,2,2,2)


# Generating a DESeq2 object:
ddsFullCountTableBvI <- DESeqDataSetFromMatrix(
  countData = BvI,
  colData = BvIm,
  design = ~ group + batch)

# DESeq2 analysis
ddsBvI <- ddsFullCountTableBvI
as.data.frame(colData(ddsBvI)) #Checking factors
ddsBvI <- DESeq(ddsBvI)
rldBvI <- rlog(ddsBvI, blind=FALSE)

pdf(file="DESeq2-PCABvI.pdf", width = 14, height = 8)
z <- plotPCA(rldBvI, intgroup=c("sample"), ntop = 500)
z + geom_label(aes(label = name))
dev.off()



# Result
resBvI <- results(ddsBvI, list(c("group")))
summary(resBvI)
resBvI <- lfcShrink(ddsBvI, coef="group", type="apeglm")
resOrderedBvI <- resBvI[order(resBvI$padj),]
write.table(as.data.frame(resOrderedBvI), file=paste("DESeq-BvI.txt", sep="."), quote=FALSE, sep="\t")

pdf(file="DESeq2-BvI.pdf")

DESeq2::plotMA(resBvI, ylim = (c(-3,3)))

##Result vulcanoplot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(resBvI$log2FoldChange, -log10(resBvI$pvalue))
plot(resBvI$log2FoldChange, -log10(resBvI$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(resBvI$log2FoldChange) > 1.5 & resBvI$padj < alpha 
text(resBvI$log2FoldChange[gn.selected],
     -log10(resBvI$padj)[gn.selected],
     lab=rownames(resBvI)[gn.selected ], cex=0.6)

plot(resBvI$log2FoldChange, -log10(resBvI$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()


# Bach ajustment and data export
vsdBvI <- vst(ddsBvI, blind=FALSE)
mat <- assay(vsdBvI)
mm <- model.matrix(~group, colData(vsdBvI))
mat <- limma::removeBatchEffect(mat, batch=vsdBvI$batch, design=mm)
assay(vsdBvI) <- mat
colData(vsdBvI)

write.table(as.data.frame(mat), file=paste("DESeq-normalized-LIMMA-BvI.txt", sep="."), quote=FALSE, sep="\t", col.names = NA, row.names = TRUE)




### Contrast6
BvIO <- dataset[,c(10:12,16:18)]
BvIOm <- Meta[c(10:12,16:18),]
BvIOm$group <- c(1,1,1,2,2,2)


# Generating a DESeq2 object:
ddsFullCountTableBvIO <- DESeqDataSetFromMatrix(
  countData = BvIO,
  colData = BvIOm,
  design = ~ group + batch)

# DESeq2 analysis
ddsBvIO <- ddsFullCountTableBvIO
as.data.frame(colData(ddsBvIO)) #Checking factors
ddsBvIO <- DESeq(ddsBvIO)
rldBvIO <- rlog(ddsBvIO, blind=FALSE)

pdf(file="DESeq2-PCABvIO.pdf", width = 14, height = 8)
z <- plotPCA(rldBvIO, intgroup=c("sample"), ntop = 500)
z + geom_label(aes(label = name))
dev.off()


# Result
resBvIO <- results(ddsBvIO, list(c("group")))
summary(resBvIO)
resBvIO <- lfcShrink(ddsBvIO, coef="group", type="apeglm")
resOrderedBvIO <- resBvIO[order(resBvIO$padj),]
write.table(as.data.frame(resOrderedBvIO), file=paste("DESeq-BvIO.txt", sep="."), quote=FALSE, sep="\t")

pdf(file="DESeq2-BvIO.pdf")

DESeq2::plotMA(resBvIO, ylim = (c(-3,3)))

##Result vulcanoplot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(resBvIO$log2FoldChange, -log10(resBvIO$pvalue))
plot(resBvIO$log2FoldChange, -log10(resBvIO$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(resBvIO$log2FoldChange) > 1.5 & resBvIO$padj < alpha 
text(resBvIO$log2FoldChange[gn.selected],
     -log10(resBvIO$padj)[gn.selected],
     lab=rownames(resBvIO)[gn.selected ], cex=0.6)

plot(resBvIO$log2FoldChange, -log10(resBvIO$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()


# Bach ajustment and data export
vsdBvIO <- vst(ddsBvIO, blind=FALSE)
mat <- assay(vsdBvIO)
mm <- model.matrix(~group, colData(vsdBvIO))
mat <- limma::removeBatchEffect(mat, batch=vsdBvIO$batch, design=mm)
assay(vsdBvIO) <- mat
colData(vsdBvIO)

write.table(as.data.frame(mat), file=paste("DESeq-normalized-LIMMA-BvIO.txt", sep="."), quote=FALSE, sep="\t", col.names = NA, row.names = TRUE)



### Contrast7
BvIN <- dataset[,c(13:18)]
BvINm <- Meta[c(13:18),]
BvINm$group <- c(1,1,1,2,2,2)


# Generating a DESeq2 object:
ddsFullCountTableBvIN <- DESeqDataSetFromMatrix(
  countData = BvIN,
  colData = BvINm,
  design = ~ group + batch)

# DESeq2 analysis
ddsBvIN <- ddsFullCountTableBvIN
as.data.frame(colData(ddsBvIN)) #Checking factors
ddsBvIN <- DESeq(ddsBvIN)
rldBvIN <- rlog(ddsBvIN, blind=FALSE)

pdf(file="DESeq2-PCABvIN.pdf", width = 14, height = 8)
z <- plotPCA(rldBvIN, intgroup=c("sample"), ntop = 500)
z + geom_label(aes(label = name))
dev.off()


# Result
resBvIN <- results(ddsBvIN, list(c("group")))
summary(resBvIN)
resBvIN <- lfcShrink(ddsBvIN, coef="group", type="apeglm")
resOrderedBvIN <- resBvIN[order(resBvIN$padj),]
write.table(as.data.frame(resOrderedBvIN), file=paste("DESeq-BvIN.txt", sep="."), quote=FALSE, sep="\t")

pdf(file="DESeq2-BvIN.pdf")

DESeq2::plotMA(resBvIN, ylim = (c(-3,3)))

##Result vulcanoplot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(resBvIN$log2FoldChange, -log10(resBvIN$pvalue))
plot(resBvIN$log2FoldChange, -log10(resBvIN$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(resBvIN$log2FoldChange) > 1.5 & resBvIN$padj < alpha 
text(resBvIN$log2FoldChange[gn.selected],
     -log10(resBvIN$padj)[gn.selected],
     lab=rownames(resBvIN)[gn.selected ], cex=0.6)

plot(resBvIN$log2FoldChange, -log10(resBvIN$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()


# Bach ajustment and data export
vsdBvIN <- vst(ddsBvIN, blind=FALSE)
mat <- assay(vsdBvIN)
mm <- model.matrix(~group, colData(vsdBvIN))
mat <- limma::removeBatchEffect(mat, batch=vsdBvIN$batch, design=mm)
assay(vsdBvIN) <- mat
colData(vsdBvIN)

write.table(as.data.frame(mat), file=paste("DESeq-normalized-LIMMA-BvIN.txt", sep="."), quote=FALSE, sep="\t", col.names = NA, row.names = TRUE)




### Contrast8
OvC <- dataset[,c(1:6)]
OvCm <- Meta[c(1:6),]
OvCm$group <- c(1,1,1,2,2,2)


# Generating a DESeq2 object:
ddsFullCountTableOvC <- DESeqDataSetFromMatrix(
  countData = OvC,
  colData = OvCm,
  design = ~ group + batch)

# DESeq2 analysis
ddsOvC <- ddsFullCountTableOvC
as.data.frame(colData(ddsOvC)) #Checking factors
ddsOvC <- DESeq(ddsOvC)
rldOvC <- rlog(ddsOvC, blind=FALSE)

pdf(file="DESeq2-PCAOvC.pdf", width = 14, height = 8)
z <- plotPCA(rldOvC, intgroup=c("sample"), ntop = 500)
z + geom_label(aes(label = name))
dev.off()


# Result
resOvC <- results(ddsOvC, list(c("group")))
summary(resOvC)
resOvC <- lfcShrink(ddsOvC, coef="group", type="apeglm")
resOrderedOvC <- resOvC[order(resOvC$padj),]
write.table(as.data.frame(resOrderedOvC), file=paste("DESeq-OvC.txt", sep="."), quote=FALSE, sep="\t")

pdf(file="DESeq2-OvC.pdf")
DESeq2::plotMA(resOvC, ylim = (c(-3,3)))
##Result vulcanoplot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(resOvC$log2FoldChange, -log10(resOvC$pvalue))
plot(resOvC$log2FoldChange, -log10(resOvC$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(resOvC$log2FoldChange) > 1.5 & resOvC$padj < alpha 
text(resOvC$log2FoldChange[gn.selected],
     -log10(resOvC$padj)[gn.selected],
     lab=rownames(resOvC)[gn.selected ], cex=0.6)

plot(resOvC$log2FoldChange, -log10(resOvC$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()


# Bach ajustment and data export
vsdOvC <- vst(ddsOvC, blind=FALSE)
mat <- assay(vsdOvC)
mm <- model.matrix(~group, colData(vsdOvC))
mat <- limma::removeBatchEffect(mat, batch=vsdOvC$batch, design=mm)
assay(vsdOvC) <- mat
colData(vsdOvC)

write.table(as.data.frame(mat), file=paste("DESeq-normalized-LIMMA-OvC.txt", sep="."), quote=FALSE, sep="\t", col.names = NA, row.names = TRUE)

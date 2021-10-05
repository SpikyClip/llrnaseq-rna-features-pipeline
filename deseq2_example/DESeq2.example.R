
# Setting up the R session ------------------------------------------------
# https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html
# https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#performing-differential-expression-testing-with-edger

#clear memory
rm(list=ls())

# if packages need to be installed
#chooseCRANmirror()
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install(c("statmod"))
#BiocManager::install(c("fastcluster"), force = TRUE)

# Loading the packages
library("repmis")
library("ggplot2")
library("reshape")
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("dendextend")
library("vsn")
library("hexbin")
library("DCGL")
library("ggrepel")


#load image if you have a saved version already (this does not apply the first time you run this script)
# load("test.Rdata")

#1. Read in the phenotypic data
pheno <- read.csv("test.pheno.csv", header=T)

# Have a look at the pheno object, list the variables
head(pheno)
table(pheno$SampleID)
table(pheno$Treatment)
table(pheno$Time)
table(pheno$Replicate)

#Read in the RNA-Seq gene counts
counts <- read.delim("test.counts.txt", header=T, row.names=1, stringsAsFactors=FALSE, na.strings="")
dim(counts) #33630 genes, 54 samples
colnames(counts)

# Boxplot of raw read counts
# Check the distribution of the raw read counts to identify potentially outlying samples
png("boxplot.png",width = 1000, height=500)
boxplot(counts[, !(colnames(counts) %in% c("Description"))]^0.2, las = 2, main = "Boxplots of raw counts to the power of 0.2", ylim = c(0,15), 
        names = colnames(counts[, !(colnames(counts) %in% c("Description"))]), cex.axis = 0.6)
abline(h = median(data.matrix(counts[, !(colnames(counts) %in% c("Description"))]))^0.2, col = "red")
legend("topleft", legend = "Median", lty = 1, col = "red")
dev.off()

# Check how many genes have evidence of expression (at least 1 read)
# At least 10,000 genes should have detectable expression in a high quality RNA-Seq data set and 
# this is the case for all thresholds
# a. all samples
table(rowSums(counts > 0) >= 54) #18771 TRUE
# b. one tissue
table(rowSums(counts > 0) >= 3) #25439 TRUE

# Check that order of count data matches pheno data 
# Should return TRUE
identical(colnames(counts[, !(colnames(counts) %in% c("Description"))]), as.character(pheno$SampleID))

## Prepare data for DESeq2
dds <- DESeqDataSetFromMatrix(countData=counts, colData=pheno, design=~1)

# Add phenotypic variables to the DESeqDataSet object (dds)
dds@colData$SampleID <- pheno$SampleID
dds@colData$Treatment <- factor(pheno$Treatment[match(dds@colData$SampleID, pheno$SampleID)], labels = c("Control", "JA"))
dds@colData$Time <- pheno$Time[match(dds@colData$SampleID, pheno$SampleID)]
dds@colData$Replicate <- factor(pheno$Replicate[match(dds@colData$SampleID, pheno$SampleID)], labels = c("rep.1","rep.2","rep.3"))

##For DESeq2: Apply pre-filtering to remove rows that have no reads or only 1 read.

min.filt <- rowSums(counts(dds)) >1
table(min.filt)
#FALSE  TRUE 
#7500 26130
dds <- dds[min.filt, ]


# NORMALISATION
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
png("vsd.png",width = 1000, height=500)
meanSdPlot(assay(vsd))
dev.off()

#Unsupervised clustering
# generate PCA plot
# colour palettes: http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
# ggrepel: https://stackoverflow.com/questions/15624656/label-points-in-geom-point
#treatment and time
# pc1 and pc2
pcadata <- plotPCA(vsd2, intgroup = c("Treatment", "Time"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))
png("pca.all.png",width = 1000, height=500)
p <- ggplot(pcadata, aes(PC1, PC2, color=Time, shape=Treatment)) +   geom_point(size=3) + scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

p + geom_label_repel(aes(label = dds@colData$SampleID),
                   box.padding   = 0.2, 
                   point.padding = 0.2,
                   segment.color = 'grey50', label.size = 0.1, max.overlaps = Inf) + theme_classic()
dev.off()

#control only
## pc1 and pc2
pcadata <- plotPCA(vsd2[,c(1:25)], intgroup = c("Time", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))
png("pca.control.png",width = 1000, height=500)
p <- ggplot(pcadata, aes(PC1, PC2, color=Time, shape=Replicate)) +   geom_point(size=3) + scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

p + geom_label_repel(aes(label = pheno$SampleID[c(1:25)]),
                     box.padding   = 0.2, 
                     point.padding = 0.2,
                     segment.color = 'grey50', max.overlaps = Inf) + theme_classic()
dev.off()

#JA only
## pc1 and pc2
pcadata <- plotPCA(vsd2[,c(28:54)], intgroup = c("Time", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))
png("pca.JA.png",width = 1000, height=500)
p <- ggplot(pcadata, aes(PC1, PC2, color=Time, shape=Replicate)) + geom_point(size=3) + scale_color_brewer(palette="Paired") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

p + geom_label_repel(aes(label = pheno$SampleID[c(28:54)]),
                     box.padding   = 0.2, 
                     point.padding = 0.2,
                     segment.color = 'grey50', max.overlaps = Inf) + theme_classic()

dev.off()

###DEG analysis
# aim of the analysis is to identify genes that are DE over time in response to JA treatment

# design 1
pheno$group <- factor(paste(pheno$Time, pheno$Treatment, sep="_"))
table(pheno$group)

# Re-create a DESeqDataSet object to include group- this is needed for the contrasts for DESeq2
# the results are similar to using ~Time+Treatment:Time
dds.2.deg <- DESeqDataSetFromMatrix(countData=counts[ , !(colnames(counts) %in% c("Description"))], colData=pheno, design = ~0+group)

# adding phenotypic information
dds.2.deg@colData$SampleID <- pheno$SampleID
dds.2.deg@colData$Treatment <- factor(pheno$Treatment[match(dds@colData$SampleID, pheno$SampleID)], labels = c("Control", "JA"))
dds.2.deg@colData$Time <- pheno$Time[match(dds@colData$SampleID, pheno$SampleID)]
dds.2.deg@colData$Replicate <- factor(pheno$Replicate[match(dds@colData$SampleID, pheno$SampleID)], labels = c("rep.1","rep.2","rep.3"))


min.filt2.deg <- rowSums(counts(dds.2.deg)) > 10
table(min.filt2.deg)
#FALSE  TRUE 
#9493 24137 

dds.2.deg <- dds.2.deg[min.filt2.deg, ]
dds.2.deg

# First, estimate library sizes - this adds an extra column in ColData called "Group sizeFactor"
dds.2.deg <- estimateSizeFactors(dds.2.deg)

#estimate dispersions
dds.2.deg <- estimateDispersions(dds.2.deg)

#Fit the negative binomial model:
DESeq2.fit <- nbinomWaldTest(dds.2.deg, betaPrior = FALSE)

# Check if any rows did not converge in beta, labelled in mcols(object)$betaConv.
table(mcols(DESeq2.fit)$betaConv) #all  converged
dds.2.deg <- dds.2.deg[which(mcols(DESeq2.fit)$betaConv), ]
dds.2.deg #all 24137 genes converged

resultsNames(DESeq2.fit)
#returns the following:
#[1] "group0h_Control"  "group0h_JA"       "group12h_Control" "group12h_JA"      "group14d_Control" "group14d_JA"      "group24h_Control" "group24h_JA"      "group2h_Control" 
#[10] "group2h_JA"       "group3d_Control"  "group3d_JA"       "group4h_Control"  "group4h_JA"       "group7d_Control"  "group7d_JA"       "group8h_Control"  "group8h_JA" 

#test contrasts
# JA vs control at 0h
JA.vs.control.0h <- results(DESeq2.fit, contrast=c("group","0h_JA","0h_Control"), alpha = 0.05)
#summary(JA.vs.control.0h)

# JA vs control at 2h
JA.vs.control.2h <- results(DESeq2.fit, contrast=c("group","2h_JA","2h_Control"), alpha = 0.05)
#summary(JA.vs.control.2h)

# JA vs control at 4h
JA.vs.control.4h <- results(DESeq2.fit, contrast=c("group","4h_JA","4h_Control"), alpha = 0.05)
#summary(JA.vs.control.4h)

# JA vs control at 8h
JA.vs.control.8h <- results(DESeq2.fit, contrast=c("group","8h_JA","8h_Control"), alpha = 0.05)
#summary(JA.vs.control.8h)

# JA vs control at 12h
JA.vs.control.12h <- results(DESeq2.fit, contrast=c("group","12h_JA","12h_Control"), alpha = 0.05)
#summary(JA.vs.control.12h)

# JA vs control at 24h
JA.vs.control.24h <- results(DESeq2.fit, contrast=c("group","24h_JA","24h_Control"), alpha = 0.05)
#summary(JA.vs.control.24h)

# JA vs control at 3d
JA.vs.control.3d <- results(DESeq2.fit, contrast=c("group","3d_JA","3d_Control"), alpha = 0.05)
#summary(JA.vs.control.3d)

# JA vs control at 7d
JA.vs.control.7d <- results(DESeq2.fit, contrast=c("group","7d_JA","7d_Control"), alpha = 0.05)
#summary(JA.vs.control.7d)

# JA vs control at 14d
JA.vs.control.14d <- results(DESeq2.fit, contrast=c("group","14d_JA","14d_Control"), alpha = 0.05)
#summary(JA.vs.control.14d)

##extract results
JA.vs.control.0h.results <- as.data.frame(cbind(JA.vs.control.0h, mcols(DESeq2.fit)[, 1:3]))
JA.vs.control.0h.sig <- JA.vs.control.0h[which(JA.vs.control.0h$padj<0.05 & JA.vs.control.0h$pvalue<0.05),]
write.csv(JA.vs.control.0h.sig, file = "JA.vs.control.0h.sig.csv", row.names = TRUE)

JA.vs.control.2h.results <- as.data.frame(cbind(JA.vs.control.2h, mcols(DESeq2.fit)[, 1:3]))
JA.vs.control.2h.sig <- JA.vs.control.2h[which(JA.vs.control.2h$padj<0.05 & JA.vs.control.2h$pvalue<0.05),]
write.csv(JA.vs.control.2h.sig, file = "JA.vs.control.2h.sig.csv", row.names = TRUE)

JA.vs.control.4h.results <- as.data.frame(cbind(JA.vs.control.4h, mcols(DESeq2.fit)[, 1:3]))
JA.vs.control.4h.sig <- JA.vs.control.4h[which(JA.vs.control.4h$padj<0.05 & JA.vs.control.4h$pvalue<0.05),]
write.csv(JA.vs.control.4h.sig, file = "JA.vs.control.4h.sig.csv", row.names = TRUE)

JA.vs.control.8h.results <- as.data.frame(cbind(JA.vs.control.8h, mcols(DESeq2.fit)[, 1:3]))
JA.vs.control.8h.sig <- JA.vs.control.8h[which(JA.vs.control.8h$padj<0.05 & JA.vs.control.8h$pvalue<0.05),]
write.csv(JA.vs.control.8h.sig, file = "JA.vs.control.8h.sig.csv", row.names = TRUE)

JA.vs.control.12h.results <- as.data.frame(cbind(JA.vs.control.12h, mcols(DESeq2.fit)[, 1:3]))
JA.vs.control.12h.sig <- JA.vs.control.12h[which(JA.vs.control.12h$padj<0.05 & JA.vs.control.12h$pvalue<0.05),]
write.csv(JA.vs.control.12h.sig, file = "JA.vs.control.12h.sig.csv", row.names = TRUE)

JA.vs.control.24h.results <- as.data.frame(cbind(JA.vs.control.24h, mcols(DESeq2.fit)[, 1:3]))
JA.vs.control.24h.sig <- JA.vs.control.24h[which(JA.vs.control.24h$padj<0.05 & JA.vs.control.24h$pvalue<0.05),]
write.csv(JA.vs.control.24h.sig, file = "JA.vs.control.24h.sig.csv", row.names = TRUE)

JA.vs.control.3d.results <- as.data.frame(cbind(JA.vs.control.3d, mcols(DESeq2.fit)[, 1:3]))
JA.vs.control.3d.sig <- JA.vs.control.3d[which(JA.vs.control.3d$padj<0.05 & JA.vs.control.3d$pvalue<0.05),]
write.csv(JA.vs.control.3d.sig, file = "JA.vs.control.3d.sig.csv", row.names = TRUE)

JA.vs.control.7d.results <- as.data.frame(cbind(JA.vs.control.7d, mcols(DESeq2.fit)[, 1:3]))
JA.vs.control.7d.sig <- JA.vs.control.7d[which(JA.vs.control.7d$padj<0.05 & JA.vs.control.7d$pvalue<0.05),]
write.csv(JA.vs.control.7d.sig, file = "JA.vs.control.7d.sig.csv", row.names = TRUE)

JA.vs.control.14d.results <- as.data.frame(cbind(JA.vs.control.14d, mcols(DESeq2.fit)[, 1:3]))
JA.vs.control.14d.sig <- JA.vs.control.14d[which(JA.vs.control.14d$padj<0.05 & JA.vs.control.14d$pvalue<0.05),]
write.csv(JA.vs.control.14d.sig, file = "JA.vs.control.14d.sig.csv", row.names = TRUE)

# save R image that can be loaded later on
save.image("test.Rdata")


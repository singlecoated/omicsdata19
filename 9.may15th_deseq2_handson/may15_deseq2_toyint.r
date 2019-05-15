######################################################
# Interaction test
######################################################

# Scenario: two mice groups (A=WT, B=KO) and two diet types (ND=Normal Diet, HFD=High Fat Diet)

library(DESeq2)
library(ggplot2)
library(pheatmap)

## Example 3: two conditions (genotype), two groups (diet), with interaction terms
set.seed(149) 
args(makeExampleDESeqDataSet)
dds <- makeExampleDESeqDataSet(n=1000,m=16,betaSD=2,interceptSD=4) # 1000 genes, 4 samples per condition/group
#save(dds,file='dds.RData')

# Assign real mus musculus names
#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
#set.seed(249)
#genes.1000 <- genes[sample(1:length(genes$gene_id),1000,rep=FALSE),] # retrieve 1000 random entrez ids
#save(genes.1000,file='genes1000mm10.RData')
load('genes1000mm10.RData')
rownames(dds) <- genes.1000$gene_id

# Explore
dim(dds)
head(assay(dds))
as.data.frame(colData(dds)) # Your sampleinfo

# Now let's add a second factor, in this case, diet
dds$group <- factor(rep(rep(c("ND","HFD"),each=4),2))
x <- as.data.frame(colData(dds))
x
table(x$condition,x$group) # Balanced, not paired

# Explore visually
boxplot(log2(assay(dds)+1),col=as.numeric(factor(x$condition))+1)
dds <- DESeq(dds)
colData(dds)
boxplot(log2(counts(dds,normalize=TRUE)+1),col=as.numeric(factor(x$condition))+1)
plot(log2(assay(dds)+1)[,1],log2(counts(dds,normalize=TRUE)+1)[,1],col=densCols(log2(assay(dds)+1)[,1],log2(counts(dds,normalize=TRUE)+1)[,1])) # Identical
rl <- rlog(dds)
plot(log2(counts(dds,normalize=TRUE)+1)[,1],assay(rl)[,1],col=densCols(log2(assay(dds)+1)[,1],assay(rl)[,1]))

# PCA basic
plotPCA(rl,intgroup=c('condition'))
plotPCA(rl,intgroup=c('group'))
plotPCA(rl,intgroup=c('condition','group'))

# Return the PCA object
pcaData <- plotPCA(rl, intgroup=c("condition", "group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot nicely

# Option 1
ggplot(pcaData, aes(PC1, PC2, color=x$condition, shape=x$group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text(aes(label=paste(rownames(x),x$condition,x$group,sep='_')),hjust=.5, vjust=-.8, size=3) +
    stat_ellipse()

# Option 2, colors only, generate dummy variable combining genotype and diet
dds$cg <- factor(paste(dds$condition,dds$group,sep='_'))
x <- colData(dds)
head(x)
rl <- rlog(dds) # recalculate to include new dummy var
pcaData <- plotPCA(rl, intgroup=c("condition", "group"), returnData=TRUE)
ggplot(pcaData, aes(PC1, PC2, color=x$cg)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text(aes(label=paste(rownames(x),x$condition,x$group,sep='_')),hjust=.5, vjust=-.8, size=3) +
    stat_ellipse()

# Model without interaction term
design(dds) <- ~ group + condition
dds <- DESeq(dds)

# Explore main effects
resultsNames(dds)
results(dds,contrast=c('condition','B','A')) # untransformed foldchanges, we will go this way
## lfcShrink(dds,contrast=c('condition','B','A'),type='apeglm') # lfcShrink for foldchanges in designs with interactions has to be of type apeglm or ashr, and need the respective packages
results(dds,contrast=c('group','HFD','ND'))
summary(results(dds,contrast=c('condition','B','A'))) # conditions explain a lot of variance
summary(results(dds,contrast=c('group','HFD','ND'))) # group explains almost nothing

# What if we want to compare individual diet effects on a certain genotype ?
# We can either add interaction terms, or more practical, our introduced dummy variable
design(dds) <- ~ group + condition + group:condition
dds <- DESeq(dds)
resultsNames(dds) # We have the interaction, but only in the direction of the factor levels, not very user friendly for nested effects more on this later

# Include dummy var
design(dds) <- ~ cg
dds <- DESeq(dds)
resultsNames(dds)

# Now we can obtain all possible combinations
results(dds,contrast=c('cg','A_HFD','A_ND'))
results(dds,contrast=c('cg','B_HFD','B_ND'))
results(dds,contrast=c('cg','A_HFD','B_ND'))
results(dds,contrast=c('cg','B_HFD','A_ND'))
summary(results(dds,contrast=c('cg','A_HFD','A_ND')))
summary(results(dds,contrast=c('cg','B_HFD','B_ND')))

# Not let's get back to the interaction term
# We are interested in knowing if the diet has a differential effect depending on the genotype (according to the PCA, seems so)
# Intuitively, we are basically interested in 'the difference of the differences', that is, the differences between foldchanges, or the slopes
design(dds) <- ~ condition + group + condition:group # full model
dds <- DESeq(dds)
resultsNames(dds) # We have the interaction, but as it is specified, it is measured using Normal Diet as the basal reference, and we want it the other way

# In order to change the reference, we have to relevel the factors (and not the factor order in the interaction, as A:B will be the same as B:A)

# check factor levels
x$condition # Reference is A, WT, that's usually ok
x$group # Reference is HFD, let's change that
dds$group <- relevel(dds$group,'ND')
dds <- DESeq(dds)
resultsNames(dds) # We have the interaction, and in the sense we want.
res <- results(dds,pAdjustMethod='BH',name='conditionB.groupHFD')
## This returns the differential effect of the diet between genotypes
## Or in other words, (B_HFD - B_ND) - (A_HFD - A_ND)

# lfcShrink(dds,coef='conditionB.groupHFD',type='apeglm') # lfcShrink for foldchanges in designs with interactions has to be of type apeglm or ashr, and need the respective packages
summary(res) # 2 genes up, 2 down, not a lot, but let's see

# Some plots and results, these are for the interactions, but you can extract main effects results tables and play with those too
# You can also check the main effects results foldchanges for those contrasts, and see what are the foldchanges for these two genes, etc etc.
res.sel <- res[!is.na(res$padj) & res$padj<0.1,]
res.sel
plotCounts(dds,gene=rownames(res.sel)[1],intgroup=c('condition','group'),col=x$cg,pch=19)
plotCounts(dds,gene=rownames(res.sel)[2],intgroup=c('condition','group'),col=x$cg,pch=19)

# Heatmap for genes selected based on unadj pvalues
res.sel <- res[!is.na(res$pvalue) & res$pvalue<0.05,]
library("pheatmap")
df <- as.data.frame(colData(dds)[,c("condition","group")])
library(pheatmap)
pheatmap(assay(dds)[rownames(res.sel),], cluster_rows=TRUE, show_rownames=TRUE,cluster_cols=FALSE, annotation_col=df,scale='row')

# To do...
# export CSV/XLSX tables with results for pairwise and interaction comparisons
# generate count plots for significant genes
# establish a 'rej' column in output tables indicating -1,0,1 if gene is differentially expressed with a foldchange / pvalue-padj of your criteria
# integrate the genes.1000 table in the output tables, and include additional information using for instance Biomart for gene name, description, mgi symbol, GO associated terms, gene biotype, etc
# explore the effect of performing lfcShrinkage on the results we find and the fold change estimations
# try to export HTML tables using the Glimma package...
# ...








######################################################
# DESeq2 template code for toy dataset analysis
# 2 factors + interaction effect
# Exploration and QC, modelling, DE, results, Viz
# UPC - MESIO. May 15th 2019
######################################################

#################################
# 1. GENERATING OUR DATA SET
# 2. WORKING WITH MOUSE TXDB
# 3. UPDATING SAMPLEINFO
# 4. COUNTS AND TRANSFORMATIONS
# 5. PRINCIPAL COMPONENT ANALYSIS
# 6. DIFFERENTIAL EXPRESSION
# 6.1 MAIN EFFECTS
# 6.2 NESTED EFFECTS
# 6.3 INTERACTION TERM
# 7. REPORTS AND VISUALIZATION
# 8. ADDITIONAL WORK...
#################################

#################################
# 1. GENERATING OUR DATA SET
#################################

# Introducing the biological scenario: 
# two mice groups (A=WT, B=KO) 
# and two diet types (ND=Normal Diet, HFD=High Fat Diet)

# We will need these libraries, install them from Bioconductor and CRAN
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Using the provided function in DESeq2 to generate an artificial read count dataset
# It will reflect our two conditions (genotype), and two groups (diet), allowing us to explore interaction terms

# Random seed
set.seed(149) # Use 149 to replicate the example using the parameters below, or change for generating a new random dataset

# You can explore a bit the function help
?makeExampleDESeqDataSet

# Generating the data set per se, you can try playing with the parameters
dds <- makeExampleDESeqDataSet(n=1000,m=16,betaSD=2,interceptSD=4) # 1000 genes, 4 samples per condition/group, betaSD=2, interceptSD=4
#save(dds,file='dds.RData') # This is the object provided with the materials

#################################
# 2. WORKING WITH MOUSE TXDB
#################################

# Since the toy example will be generated with gene1...geneN names, let's assign real mus musculus names
# For that purpose, we will take advantage of the TxDB genome objects in Bioconductor, in our case, the one for mouse
# You can find more information about this extensive package and its uses in 
# http://bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm10.knownGene.html

# For speeding up things, I already extracted 1000 random genes from this dataset, stored in the genes.1000.RData file
# The code needed to install, load and do this is provided commented for your reference or if you want to use it.

#################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
#set.seed(249)
#genes.1000 <- genes[sample(1:length(genes$gene_id),1000,rep=FALSE),] # retrieve 1000 random entrez ids
#save(genes.1000,file='genes1000mm10.RData')
##################################################################

# Loading the saved gene information, in GRanges format, remember to give the correct path for loading it
load('genes1000mm10.RData')
class(genes.1000)
head(genes.1000)

# Assign these genes names to our toy dataset
rownames(dds) <- genes.1000$gene_id

# Explore the resulting object
dim(dds)
head(assay(dds))

#################################
# 3. UPDATING SAMPLEINFO
#################################

# By creating the toy dataset, we have our samples and our main experimental factor
as.data.frame(colData(dds)) # Your sampleinfo

# So let's add a second factor, in this case, diet, in this case, we can do it directly to our dds object
dds$group <- factor(rep(rep(c("ND","HFD"),each=4),2))

# For practical purposes, let's store this sample information file in an auxiliary variable x (not a very informative name, but it will do)
x <- as.data.frame(colData(dds))
x
table(x$condition,x$group) # We have a balanced design

#################################
# 4. COUNTS AND TRANSFORMATIONS
#################################

# Now that we have our count matrix and experimental information, we can start to explore visually our log2 pseudocounts
boxplot(log2(assay(dds)+1),col=as.numeric(factor(x$condition))+1)

# In order to try some normalization methods, we have to run DESeq over our DESeq dataset
dds <- DESeq(dds)

# We can see that we have now estimated size factors that will be used for normalization
colData(dds)

# Now, we can play with DESeq normalized counts (using size factors, rlog, etc)
boxplot(log2(counts(dds,normalize=TRUE)+1),col=as.numeric(factor(x$condition))+1)

# Let's explore how log2 pseudocounts and log2 normalized counts relate
plot(log2(assay(dds)+1)[,1],log2(counts(dds,normalize=TRUE)+1)[,1],col=densCols(log2(assay(dds)+1)[,1],log2(counts(dds,normalize=TRUE)+1)[,1])) # Identical

# And now, let's see how they look compared to the rlog transformation
rl <- rlog(dds)
plot(log2(counts(dds,normalize=TRUE)+1)[,1],assay(rl)[,1],col=densCols(log2(assay(dds)+1)[,1],assay(rl)[,1]))

#################################
# 5. PRINCIPAL COMPONENT ANALYSIS
#################################

# As we know, the PCA is a basic technique for data exploration and quality control, sample assessment
# PCA basic, remember, there may be additional variables (biological or technical, known or not) introducing other sources of variability
# A PCA is also a good method to explore these potencial batch effects, which could involve additional actions in downstream analysis
plotPCA(rl,intgroup=c('condition'))
plotPCA(rl,intgroup=c('group'))
plotPCA(rl,intgroup=c('condition','group'))

# In order to generate nicer plots, we can tell the function to return the values
# Return the PCA object and compute percentage of variance explained by the components
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
# This second option reflects clearly what's happening in our data, the different variability between genotypes,
# and between diets within genotypes, let's create a dummy variable combining our genotypes and diets
# for now we will use that to assign 4 different colors
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

#################################
# 6.1 DIFFERENTIAL EXPRESSION
# MAIN EFFECTS
#################################

# Let's create a model without interaction term, just our two variables
design(dds) <- ~ group + condition
dds <- DESeq(dds)

# This will allow us to explore main effects, global effects between genotypes, or diets
resultsNames(dds)

# To obtain the results tables for the different contrasts, we do as specified in the DESeq2 manual
results(dds,contrast=c('condition','B','A')) # untransformed foldchanges (not shrinked), for now we will go this way

# In case we want to perform log2 fold change shrinkage to adjust high foldchanges for low count genes
## lfcShrink(dds,contrast=c('condition','B','A'),type='apeglm')
# Remember, in designs with interactions has to be of type apeglm or ashr, and need the respective packages to be installed

# On this way we can obtain main effects results, and assess how the number of DE genes assess to the variance we see in the PCA
results(dds,contrast=c('group','HFD','ND'))
summary(results(dds,contrast=c('condition','B','A'))) # conditions explain a lot of variance
summary(results(dds,contrast=c('group','HFD','ND'))) # group explains almost nothing

#################################
# 6.2 DIFFERENTIAL EXPRESSION
# NESTED EFFECTS
#################################

# But what if we want to compare individual diet effects within a certain genotype ?
# We can either add interaction terms, or more practical, our introduced dummy variable
# Let's try first with the interaction term
design(dds) <- ~ group + condition + group:condition
dds <- DESeq(dds)
resultsNames(dds) 
# We have the interaction, but only in the direction of the factor levels, not very user friendly for nested effects, more on this later

# Now let's instead Include dummy var we created before using the line commented below
#dds$cg <- factor(paste(dds$condition,dds$group,sep='_'))
# And let's build a model based on that new dummy variable
design(dds) <- ~ cg
dds <- DESeq(dds)
resultsNames(dds) # Now we have results for contrasts also within groups, so we can assess all combinations

# So, to obtain all possible combinations
results(dds,contrast=c('cg','A_HFD','A_ND'))
results(dds,contrast=c('cg','B_HFD','B_ND'))
results(dds,contrast=c('cg','A_HFD','B_ND'))
results(dds,contrast=c('cg','B_HFD','A_ND'))
summary(results(dds,contrast=c('cg','A_HFD','A_ND')))
summary(results(dds,contrast=c('cg','B_HFD','B_ND')))

#################################
# 6.3 DIFFERENTIAL EXPRESSION
# INTERACTION TERM
#################################

# This is already pretty interesting from a biological point of view as a resuls for an RNASeq D.E. report
# But probably we are also interested in the interaction term
# We are interested in knowing if the diet has a differential effect depending on the genotype (according to the PCA, seems so)
# Intuitively, we are basically interested in 'the difference of the differences', that is, the differences between foldchanges, or the slopes
# For that, we just build the complete model including main effects and interaction, just as you would do with limma
design(dds) <- ~ condition + group + condition:group # full model
dds <- DESeq(dds)
resultsNames(dds) 
# We have the interaction, but as it is specified, it is measured using High Fat Diet as the basal reference, and we want it the other way
# In order to change the reference, we have to relevel the factors (and not the factor order in the interaction, as A:B will be the same as B:A)
# Retrieving the correct interaction sign always needs the factor levels to be in the adequate order for our definition of 'reference'
# And remember, in DE analysis, X vs Y always means that Y is your reference sample 

# So let's check factor levels and rename accordingly
dds$condition # Reference is A, WT, that's usually ok
dds$group # Reference is HFD, as it's done by default in alphabetical order, let's change that
dds$group <- relevel(dds$group,'ND')
dds <- DESeq(dds)
resultsNames(dds) # We have the interaction, and in the sense we want.

# Let's extract and see what we have
res <- results(dds,pAdjustMethod='BH',name='conditionB.groupHFD')
## This returns the differential effect of the diet between genotypes
## Or in other words, (B_HFD - B_ND) - (A_HFD - A_ND)

# In case you'd like to play with shrinked fold changes
# lfcShrink(dds,coef='conditionB.groupHFD',type='apeglm') # lfcShrink for foldchanges in designs with interactions has to be of type apeglm or ashr, and need the respective packages
summary(res) # 2 genes up, 2 outliers, not a lot, but let's see

#################################
# 7. REPORTS AND VISUALIZATION
#################################

# Some plots and results, these are for the interactions, but you can extract main effects results tables and play with those too
# You can also check the main effects results foldchanges for those contrasts, and see what are the foldchanges for these two genes, etc etc.
# Also, by returning the values (as we did in the PCA), you can later make nicer and more complete representations using ggplot
res.sel <- res[!is.na(res$padj) & res$padj<0.1,] # BH adjusted pvalues by default
res.sel
plotCounts(dds,gene=rownames(res.sel)[1],intgroup=c('condition','group'),col=x$cg,pch=19)
plotCounts(dds,gene=rownames(res.sel)[2],intgroup=c('condition','group'),col=x$cg,pch=19)

# These plots summarize in a simple way what we can expect from an interaction effect
# Remember though, that reporting only interaction effects is not usually very informative,
# So, it is also usually important to report individual effects tables (between genotypes/conditions)

# Additionally, it is usually a good approach to use heatmaps to visualize how things change and behave betweeen samples
# In our case we will generate a pheatmap for genes selected based on unadj pvalues (to get a few more than just two)
res.sel <- res[!is.na(res$pvalue) & res$pvalue<0.05,]
library("pheatmap")
df <- as.data.frame(colData(dds)[,c("condition","group")])
library(pheatmap)
pheatmap(assay(dds)[rownames(res.sel),], cluster_rows=TRUE, show_rownames=TRUE,cluster_cols=FALSE, annotation_col=df,scale='row')

#################################
# 8. ADDITIONAL WORK...
#################################

# Now, some open questions for you to do, explore, complete the code, etc
# For instance...
# Use this code template and complete it for your own Rmarkdown document
# export CSV/XLSX tables with results for pairwise and interaction comparisons (see library openxlsx)
# generate count plots for significant genes
# establish a 'rej' column in output tables indicating -1,0,1 if gene is differentially expressed with a foldchange / pvalue-padj of your criteria
# integrate the genes.1000 table in the output tables
# Include additional information using for instance Biomart for gene name, description, mgi symbol, GO associated terms, gene biotype, etc
# Generate additional quality control or exploratory plots such as MA and others
# explore the effect of performing lfcShrinkage on the results we find and the fold change estimations
# try to export HTML tables using the Glimma package...
# perform gene ontology related enrichment analysis (using the online DAVID tool, GSEA using the camera limma function, etc)
# References
# The huge and complete DESeq2 vignette: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# openxlsx library: https://cran.r-project.org/web/packages/openxlsx/index.html
# TxDB mmusculus http://bioconductor.org/packages/release/data/annotation/html/TxDb.Mmusculus.UCSC.mm10.knownGene.html
# biomaRt: https://bioconductor.org/packages/release/bioc/html/biomaRt.html
# Glimma interactive HTML reports: https://bioconductor.org/packages/release/bioc/html/Glimma.html
# DAVID gene set enrichment tests: https://david.ncifcrf.gov/
# GSEA enrichment analysis using the limma camera function: https://bioconductor.org/packages/release/bioc/html/limma.html 




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Execute DEseq2 on the a given featurecounts table
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# list all lib paths availabe for R
#libpaths <- .libPaths()
wd <- getwd()
# select the lib path specific to the R installed in snakemake conda
#conda_R_libpath <- libpaths[grep(paste(wd, ".snakemake/conda/", sep = "/"), libpaths)]

#print("R library path is set to: ")
#print(conda_R_libpath)

# forget all other available R lib paths except conda_R_libpath
#assign(".lib.loc", conda_R_libpath, envir = environment(.libPaths))

# load the installed packages
library(DESeq2)#, lib.loc = conda_R_libpath)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

## Differential expression analysis with DESeq2
## 12.1_deseq2.R <workdir> <samplesheet> <featurecount>  <design_formula>

args <- commandArgs(trailingOnly = TRUE)
workdir <- as.character(args[1])
samplesheet <- as.character(args[2])
featurecount <- as.character(args[3])
formula <- as.character(args[4]) # put testing condition last
print(formula)
# Import & pre-process ----------------------------------------------------
setwd(workdir)
#column data input
coldata <- read.csv(samplesheet, sep = '\t', row.names = "Sample_ID")
#count matrix input
counts <- as.matrix(read.table(featurecount, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE))
# Analysis with DESeq2 ----------------------------------------------------
print("Launching DESeq2")
#Creating the DESeq2 DESeqDataSet (dds) object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = as.formula(formula))
# Run the DESeq pipeline
dds <- DESeq(dds)
# Get differential expression results
res <- results(dds, contrast = c("Condition", "Edited", "Original"))
# Order by adjusted p-value
res <- res[order(res$padj),]
# Merge with normalized count data
resdata <- merge(res, counts(dds, normalized = TRUE), by = "row.names")
rownames(resdata) <- resdata$Row.names
resdata$Row.names <- NULL
write.table(resdata, file = "diffexpr-results.tsv", sep = '\t', col.names = NA)

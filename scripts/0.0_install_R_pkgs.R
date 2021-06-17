#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# All the required R packages with specific versions are installed through this
# script inside conda env R library directory. 
# Note: These packages are not installed using conda because it creats some conflicts
# other envs. Also Conda currently do not have full support for all the R packages
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

args <- commandArgs(trailingOnly = TRUE)
threads <- as.numeric(args[1])
outcheckfile <- as.character(args[2])

libpaths <- .libPaths()
wd <- getwd()
# select the lib path specific to the R installed in snakemake conda
conda_R_libpath <- libpaths[grep(paste(wd, ".snakemake/conda/", sep = "/"), libpaths)]

# forget all other available R lib paths except conda_R_libpath
assign(".lib.loc", conda_R_libpath, envir = environment(.libPaths))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

options(Ncpus = threads)

# install required packages, if not already installed in conda_R_libpath
if (!any(installed.packages(priority = "NA", lib = conda_R_libpath)[, 1] %in% "versions")) {
  install.packages("https://cran.r-project.org/src/contrib/versions_0.3.tar.gz", repos = NULL, type = "source", lib = conda_R_libpath)
  library("versions", lib.loc = conda_R_libpath)
}else {
  library("versions", lib.loc = conda_R_libpath)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# required packages for deseq2 discussed here: https://bioconductor.riken.jp/packages/3.8/bioc/manuals/DESeq2/man/DESeq2.pdf
# we need the exact same version of R packages as hard-coded below (they are available from bioconductor 3.8 version)
# other way to get the right package is to type in the package name in the "Search table:" text box in the following webpage -
# https://bioconductor.org/packages/3.8/BiocViews.html#___Software
# if we change the version they might not be available in the source to download

cran_pkgs = data.frame(pkg = c("Rcpp", "foreign", "Hmisc", "ggplot2", "locfit", "RCurl", "matrixStats", "xtable", "XML", "DBI", "RSQLite", "futile.logger", "snow", "RcppArmadillo"),
                       version = c("1.0.5", "0.8-75", "4.4-2", "3.3.3", "1.5-9.2", "1.98-1.2", "0.57.0", "1.8-4", "3.99-0.1", "1.1.1", "2.2.3", "1.4.3", "0.4-3", "0.10.2.2.0"))

bioc_pkgs <- data.frame(pkg = c("BiocGenerics", "S4Vectors", "IRanges", "zlibbioc", "XVector", "GenomeInfoDbData", "GenomeInfoDb", "GenomicRanges", "BiocParallel",
                                "DelayedArray", "Biobase", "AnnotationDbi", "annotate", "SummarizedExperiment", "genefilter", "geneplotter", "DESeq2"),
                        source = c("https://bioconductor.org/packages/3.8/bioc/src/contrib/BiocGenerics_0.28.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/S4Vectors_0.20.1.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/IRanges_2.16.0.tar.gz",
                                   "https://bioconductor.org/packages/3.5/bioc/src/contrib/zlibbioc_1.22.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/XVector_0.22.0.tar.gz",
                                   "http://bioconductor.riken.jp/packages/3.8/data/annotation/src/contrib/GenomeInfoDbData_1.2.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/GenomeInfoDb_1.18.2.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/GenomicRanges_1.34.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/BiocParallel_1.16.6.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/DelayedArray_0.8.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/Biobase_2.42.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/AnnotationDbi_1.44.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/annotate_1.60.1.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/SummarizedExperiment_1.12.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/genefilter_1.64.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/geneplotter_1.60.0.tar.gz",
                                   "https://bioconductor.org/packages/3.8/bioc/src/contrib/DESeq2_1.22.2.tar.gz"))

# check and install required cran packages if not already installed at conda_R_libpath
lapply(as.character(cran_pkgs$pkg),
       function(x) {
         if (!any(installed.packages(priority = "NA", lib = conda_R_libpath)[, 1] %in% x)) {
           install.versions(x, as.character(subset(cran_pkgs, pkg == x)$version), lib = conda_R_libpath)
         }
       })

# check and install required bioconductor packages if not already installed at conda_R_libpath
lapply(as.character(bioc_pkgs$pkg),
       function(x) {
         if (!any(installed.packages(priority = "NA", lib = conda_R_libpath)[, 1] %in% x)) {
           install.packages(as.character(subset(bioc_pkgs, pkg == x)$source), repos = NULL, type = "source", lib = conda_R_libpath)
         }
       })

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#required for eSNPKaryotype
pkgs_esnp = data.frame(pkg = c("usethis", "devtools", "zoo", "gplots"), version = c("1.6.3", "2.3.2", "1.8-7", "3.1.1"))

lapply(as.character(pkgs_esnp$pkg),
       function(x) {
         if (!any(installed.packages(priority = "NA", lib = conda_R_libpath)[, 1] %in% x)) {
           install.versions(x, as.character(subset(pkgs_esnp, pkg == x)$version), lib = conda_R_libpath)
         }
       })


library(devtools, lib.loc = conda_R_libpath)
Sys.setenv(TAR = "/bin/tar")
install_github("BenvenLab/eSNPKaryotyping/eSNPKaryotyping", lib = conda_R_libpath)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# required packages step 7 concatenate vcf
pkgs_concat = data.frame(pkg = c("data.table", "gtools"), version = c("1.13.6", "3.8.2"))
#pkgs_concat = data.frame(pkg = c("data.table"),  version = c("1.13.6"))

lapply(as.character(pkgs_concat$pkg),
       function(x) {
         if (!any(installed.packages(priority = "NA", lib = conda_R_libpath)[, 1] %in% x)) {
           install.versions(x, as.character(subset(pkgs_concat, pkg == x)$version), lib = conda_R_libpath)
         }
       })

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

write(x = 'All packages are installed', file = outcheckfile)

# to check all the installed packages, do
# installed.packages( priority = "NA" , lib = conda_R_libpath)[,1]

# ro remove the installed packages, do
# remove.packages(installed.packages( priority = "NA" , lib=conda_R_libpath)[,1], lib=conda_R_libpath)

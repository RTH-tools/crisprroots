#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# install/load required libraries
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# list all lib paths availabe for R
#libpaths = .libPaths()
wd = getwd()
# select the lib path specific to the R installed in snakemake conda
#conda_R_libpath = libpaths[grep(paste(wd, ".snakemake/conda/", sep = "/"), libpaths)]

#print("R library path is set to: ")
#print(conda_R_libpath)

# forget all other available R lib paths except conda_R_libpath
#assign(".lib.loc", conda_R_libpath, envir = environment(.libPaths))

# load the installed packages
library(data.table)#, lib.loc = conda_R_libpath)
library(gtools)#, lib.loc = conda_R_libpath)

# user inputs
args = commandArgs(TRUE)
dir = args[1]
outfile = args[2]

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# list .vcf file to be concatenated
f = Sys.glob(paste(dir, "*_filtered.vcf", sep = "/"))
f = gsub("chrM", "chrZM", f)
f = gsub("chrZM", "chrM", f[mixedorder(f)])

# get the header lines for the vcf from the .vcf containing all sample columns
data_header_list = lapply(f, function(x) { o = readLines(x); o[grep("##", o)] }); names(data_header_list) = f
# read the content of all .vcf files
data_list = lapply(f, function(x) fread(x, skip = "#CHROM", header = F)); names(data_list) = f
data_dim = rbindlist(lapply(names(data_list), function(x) data.table(file = x, t(dim(data_list[[x]])))))
header_file = data_dim[which.max(data_dim$V2)]$file
data_concat = unique(as.data.frame(rbindlist(data_list, fill = TRUE)))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# generate a merged .vcf file
sink(outfile)
cat(paste(paste(data_header_list[[header_file]], collapse = "\n"), "\n", sep = ""))
cat(paste(unlist(lapply(1:nrow(data_concat), function(x) paste(unlist(data_concat[x,]), collapse = "\t"))), collapse = "\n"))
sink()
closeAllConnections()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
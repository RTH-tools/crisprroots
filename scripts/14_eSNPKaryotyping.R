#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# eSNPKaryotyping v.1.0 pipeline from https://github.com/BenvenLab/eSNPKaryotyping.
# Weissbein et al, 2015, modified to use with CRISPRroots output files
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# list all lib paths availabe for R
libpaths = .libPaths()
wd = getwd()
# select the lib path specific to the R installed in snakemake conda
conda_R_libpath = libpaths[grep(paste(wd, ".snakemake/conda/", sep = "/"), libpaths)]

print("R library path is set to: ")
print(conda_R_libpath)

# forget all other available R lib paths except conda_R_libpath
assign(".lib.loc", conda_R_libpath, envir = environment(.libPaths))

# load the installed packages
library(devtools, lib.loc = conda_R_libpath)
library(zoo, lib.loc = conda_R_libpath)
library(gplots, lib.loc = conda_R_libpath)
#library(eSNPKaryotyping, lib.loc = conda_R_libpath)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

args <- commandArgs(trailingOnly = TRUE)
vcf_edited_dir <- as.character(args[1])
organism <- as.character(args[2])
dbSNP_dir <- as.character(args[3])
ref_picard_dic <- as.character(args[4])
bams_path <- as.character(args[5])
out_file_prefix <- as.character(args[6])

DeletionTable2 <- function(Directory, Table, dbSNP_Data_Directory, Genome_Fa_dict, Organism, bam_Directory) {
  print("Reading SNPs table")
  i = 1

  if (Organism == "Human") { mx = 25 }
  if (Organism == "Mouse") { mx = 22 }
  while (i < mx) {
    print(paste("Chromosome:", i))
    chrTable = read.delim(paste(dbSNP_Data_Directory, "/", "Edited_Common_chr", i, sep = ""))
    if (i == 1) { snpTable = chrTable }
    if (i > 1) { snpTable = rbind(snpTable, chrTable) }
    i = i + 1
  }

  table2 = Table
  colnames(table2)[2] = "start"
  print("Merging Tables")
  x = merge(snpTable, table2, by = c("chr", "start"), all.x = T)
  x = x[order(x$chr, x$start),]


  setwd(Directory)
  dict = read.csv(Genome_Fa_dict, as.is = T)
  dict_type = grep(pattern = "chr", x = dict[1, 1])
  if (length(dict_type) == 0) {
    dict_type = 0
  }
  i = 1
  while (i < mx) {
    print(paste("Chromosome ", i, "| ", Sys.time(), sep = ""))
    x1 = x[x$chr == i,]

    if (dict_type == 0) {
      loc = paste(x1$chr[1], ":", x1$start[1], "-", x1$start[dim(x1)[1]],
                  sep = "")
      if (Organism == "Human") {
        if (i == 23) {
          loc = paste("X:", x1$start[1], "-", x1$start[dim(x1)[1]],
                      sep = "")
        }
        if (i == 24) {
          loc = paste("Y:", x1$start[1], "-", x1$start[dim(x1)[1]],
                      sep = "")
        }
      }
      if (Organism == "Mouse") {
        if (i == 20) {
          loc = paste("X:", x1$start[1], "-", x1$start[dim(x1)[1]],
                      sep = "")
        }
        if (i == 21) {
          loc = paste("Y:", x1$start[1], "-", x1$start[dim(x1)[1]],
                      sep = "")
        }
      }
    }
    if (dict_type == 1) {
      loc = paste("chr", x1$chr[1], ":", x1$start[1],
                  "-", x1$start[dim(x1)[1]], sep = "")
      if (Organism == "Human") {
        if (i == 23) {
          loc = paste("chrX:", x1$start[1], "-", x1$start[dim(x1)[1]],
                      sep = "")
        }
        if (i == 24) {
          loc = paste("chrY:", x1$start[1], "-", x1$start[dim(x1)[1]],
                      sep = "")
        }
      }
      if (Organism == "Mouse") {
        if (i == 20) {
          loc = paste("chrX:", x1$start[1], "-", x1$start[dim(x1)[1]],
                      sep = "")
        }
        if (i == 21) {
          loc = paste("chrY:", x1$start[1], "-", x1$start[dim(x1)[1]],
                      sep = "")
        }
      }

      command = paste("samtools depth -r ", loc, " ", bam_Directory, "/Aligned.Sorted.bam > reads-per-position.txt", sep = "")
      print(command)
      system(command)
      system("awk -F \" \" '($3 >20){print $0}' reads-per-position.txt >reads-per-position2.txt")
      chr = read.delim("reads-per-position2.txt", header = F)
      colnames(chr) = c("chr", "start", "Depth")
      x2 = merge(x1, chr, by = "start")
      x2 = cbind(x2, Depth_group = x2$Depth)
      x2$Depth_group[x2$Depth < 50] = "20-50"
      x2$Depth_group[x2$Depth > 49 & x2$Depth < 100] = "50-100"
      x2$Depth_group[x2$Depth > 99 & x2$Depth < 200] = "100-200"
      x2$Depth_group[x2$Depth > 199 & x2$Depth < 500] = "200-500"
      x2$Depth_group[x2$Depth > 499] = ">500"
      x2$Depth_group = as.factor(x2$Depth_group)
      if (i == 1) { tbl = x2 }
      if (i > 1) { tbl = rbind(tbl, x2) }
      i = i + 1
    }
  }
  print("Writing Table")
  tbl[is.na(tbl)] = 0
  write.table(tbl, "Deletions.txt", sep = "\t", row.names = F, quote = F)
  return(tbl)
}

EditVCF2 <- function(Directory, Organism) {
  print("Editing VCF File")
  Dir = Directory
  file = "variants_filtered.vcf"
  path = paste(Dir, file, sep = "/")
  readData = read.delim(path, as.is = T)
  readData = as.character(readData[-c(1:which(readData == "#CHROM") - 1), 1])
  jump = 10
  startChr = 1 + jump
  startPos = 2 + jump
  startInfo = 10 + jump
  len = length(readData)
  chrRegex = "^chr(\\w+)$"
  infoRegex = "^([01])\\/([01]):(\\d+)\\,(\\d+):(\\d+):\\d+:\\d+\\,\\d+\\,\\d+$"
  chrVector = readData[startChr]
  posVector = readData[startPos]
  infoVector = readData[startInfo]
  while (startInfo + jump < len) {
    startChr = startChr + jump
    startPos = startPos + jump
    startInfo = startInfo + jump
    chrVector = append(chrVector, readData[startChr])
    posVector = append(posVector, readData[startPos])
    infoVector = append(infoVector, readData[startInfo])
  }
  chrNum = gsub(chrRegex, "\\1", chrVector)
  if (Organism == "Human") {
    chrNum[chrNum == "X"] = "23"
    chrNum[chrNum == "Y"] = "24"
  }
  if (Organism == "Mouse") {
    chrNum[chrNum == "X"] = "20"
    chrNum[chrNum == "Y"] = "21"
  }
  chrNum = as.numeric(chrNum)
  Karyotape = 10 * abs(as.numeric(gsub(infoRegex, "\\1", infoVector)) -
                         as.numeric(gsub(infoRegex, "\\2", infoVector)))
  AD1 = as.numeric(gsub(infoRegex, "\\3", infoVector))
  AD2 = as.numeric(gsub(infoRegex, "\\4", infoVector))
  DP = as.numeric(gsub(infoRegex, "\\5", infoVector))
  posVector = as.numeric(posVector)
  table = data.frame(chr = chrNum, position = posVector, AD1 = AD1,
                     AD2 = AD2, DP = DP, Karyotape = Karyotape)
  fileName = "variantTable.csv"
  pathToSave = paste(Dir, fileName, sep = "/")
  table[is.na(table)] = 0
  write.table(table, pathToSave, sep = "\t", row.names = F,
              quote = F)
  return(table)
}

Edit_dbSNP_Files2 <- function(Directory, File_Name, Organism) {
  if (Organism == "Human") {
    for (i in 1:24) {
      output = paste(Directory, "/Edited_Common_chr", i, sep = "")
      if (!file.exists(output)) {
        print(paste("Chromosme:", i))
        p = paste(Directory, File_Name, i, ".tsv", sep = "")
        readData = read.delim(p, as.is = T, header = F,
        )
        chrRegex = "^chr(\\w+)$"
        chrNum = gsub(chrRegex, "\\1", readData$V1)
        chrNum[chrNum == "X"] = "23"
        chrNum[chrNum == "Y"] = "24"
        snpTable = data.frame(chr = chrNum, start = readData$V4,
                              end = readData$V5, snp = array(0, length(readData$V1)))
        snpTable = snpTable[order(as.numeric(snpTable[,
                                               1]), as.numeric(snpTable[, 2])),]
        d = duplicated(snpTable, MARGIN = 1)
        snpTable = snpTable[d == FALSE,]
        size = snpTable$end - snpTable$start
        snpTable = snpTable[size == 0,]
        write.table(snpTable, output, sep = "\t", row.names = F,
                    quote = F)
      }
    }
  }
  if (Organism == "Mouse") {
    for (i in 1:21) {
      output = paste(Directory, "Edited_Common_chr", i, sep = "")
      if (!file.exists(output)) {
        print(paste("Chromosme:", i))
        p = paste(Directory, File_Name, i, sep = "")
        readData = read.delim(p, as.is = T, header = F,
        )
        chrRegex = "^chr(\\w+)$"
        chrNum = gsub(chrRegex, "\\1", readData$V1)
        chrNum[chrNum == "X"] = "20"
        chrNum[chrNum == "Y"] = "21"
        snpTable = data.frame(chr = chrNum, start = readData$V4,
                              end = readData$V5, snp = array(0, length(readData$V1)))
        snpTable = snpTable[order(as.numeric(snpTable[,
                                               1]), as.numeric(snpTable[, 2])),]
        d = duplicated(snpTable, MARGIN = 1)
        snpTable = snpTable[d == FALSE,]
        size = snpTable$end - snpTable$start
        snpTable = snpTable[size == 0,]
        write.table(snpTable, output, sep = "\t", row.names = F,
                    quote = F)
      }
    }
  }
}

PlotGenome2 <- function(orderedTable, Window, Ylim, PValue, Organism) {
  par(mar = c(5, 4, 4, 2))
  window = Window
  if (Organism == "Human") {
    centromere_pos = c(125, 93.3, 91, 50.4, 48.4, 61, 59.9,
                       45.6, 49, 40.2, 53.7, 35.8, 17.9, 17.6, 19, 36.6,
                       24, 17.2, 26.5, 27.5, 13.2, 14.7, 60.6, 12.5)
    centromere_pos = centromere_pos * 1e+06
    chr_size = c(248956422, 242193529, 198295559, 190214555,
                 181538259, 170805979, 159345973, 145138636, 138394717,
                 133797422, 135086622, 133275309, 114364328, 107043718,
                 101991189, 90338345, 83257441, 80373285, 58617616,
                 64444167, 46709983, 50818468, 156040895, 57227415)
  }
  if (Organism == "Mouse") {
    chr_size = c(195471971, 182113224, 160039680, 156508116,
                 151834684, 149736546, 145441459, 129401213, 124595110,
                 130694993, 122082543, 120129022, 120421639, 124902244,
                 104043685, 98207768, 94987271, 90702639, 61431566,
                 171031299, 91744698)
  }
  chr_total = 0
  for (i in 1:(length(chr_size) - 1)) {
    chr_total = c(chr_total, sum(chr_size[1:i]))
  }
  genome_size = sum(chr_size)
  orderedTable$position = orderedTable$position + chr_total[orderedTable$chr]
  col = rollmedian(orderedTable$chr, window) %% 2 + 1
  col[col == 1] = "dodgerblue4"
  col[col == 2] = "dodgerblue1"
  plot(rollmedian(orderedTable$position, window), rollmedian(orderedTable$MajorMinor,
                                                             window), col = "dimgrey", pch = 15, cex = 0.4, ylim = c(1,
                                                                                                                     Ylim), ylab = "Allelic Ratio", typ = "l", xlab = "",
       xaxt = "n", xlim = c(1, genome_size))
  points(rollmedian(orderedTable$position, window), rollmedian(orderedTable$MajorMinor,
                                                               window), col = col, pch = 15, cex = 0.4, ylim = c(1,
                                                                                                                 Ylim), , xlab = "Chromosomal Position")
  if (PValue == TRUE) {
    g = orderedTable$MajorMinor

    ttest = function(x) {
      ttt = t.test(x, g, alternative = "greater")$p.value
      return(ttt)
    }

    tt = rollapply(g, width = window, FUN = ttest)
    tt = p.adjust(tt, "fdr")
    tt = -1 * log(tt, 10)
    cl = colorpanel(50000, "white", "grey", "red")
    pmax = ceiling(max(tt))
    if (pmax < 2) {
      pmax = 2
    }
    l = length(rollmedian(orderedTable$MajorMinor, window))
    xx = rollmedian(orderedTable$position, window)
    for (i in 1:l) {
      lines(c(xx[i], xx[i]), c(2.5, 2.7), col = cl[round(tt[i] *
                                                           (50000 / pmax))])
      if (tt[i] > 2) {
        lines(c(xx[i], xx[i]), c(2.85, 2.75), col = "gray34")
      }
    }
    left = 5e+07
    for (i in 1:100) {
      rect(left, 3.1, left + 5e+06, 3.3, xpd = T, col = cl[i *
                                                             500], border = NA)
      left = left + 5e+06
    }
    left = 5e+07
    lines(c(left, 5.5e+08), c(3.1, 3.1), xpd = T)
    lines(c(left, 5.5e+08), c(3.3, 3.3), xpd = T)
    lines(c(left, left), c(3.1, 3.3), xpd = T)
    lines(c(5.5e+08, 5.5e+08), c(3.1, 3.3), xpd = T)
    mtext(side = 3, at = 0, "0", line = 0.8, cex = 0.8)
    mtext(side = 3, at = 6e+08, pmax, line = 0.8, cex = 0.8)
    mtext(side = 3, at = 6e+08 / 2, "-Log10(P-Value)", line = 2.5,
          cex = 0.8)
  }
  if (Organism == "Human") {
    for (i in 1:24) {
      if (i > 1) {
        lines(c(chr_total[i], chr_total[i]), c(1, Ylim),
              col = "gray48")
      }
      lines(c(chr_total[i] + centromere_pos[i], chr_total[i] +
        centromere_pos[i]), c(1, Ylim), col = "gray55",
            lty = 4)
      text(chr_total[i] + chr_size[i] / 2, Ylim - 0.1, i,
           cex = 0.8)
    }
  }
  if (Organism == "Mouse") {
    for (i in 1:21) {
      if (i > 1) {
        lines(c(chr_total[i], chr_total[i]), c(1, Ylim),
              col = "gray48")
      }
      text(chr_total[i] + chr_size[i] / 2, Ylim - 0.1, i,
           cex = 0.8)
    }
  }
  sum = 0
  for (i in 1:(length(chr_size) - 1)) {
    lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1,
                                                       1), lwd = 10, col = "black")
    lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1,
                                                       1), lwd = 8, col = "gray50")
    lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1,
                                                       1), lwd = 7, col = "gray53")
    lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1,
                                                       1), lwd = 6, col = "gray59")
    lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1,
                                                       1), lwd = 4, col = "gray75")
    lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1,
                                                       1), lwd = 2, col = "gray85")
    lines(c(sum + 2e+07, sum + chr_size[i] - 2e+07), c(1,
                                                       1), lwd = 1, col = "gray90")
    if (Organism == "Human") {
      lines(c(sum + centromere_pos[i], sum + centromere_pos[i]),
            c(1.01, 0.99), col = "grey13", lwd = 2)
    }
    sum = sum + chr_size[i]
  }
  lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 10,
        col = "black")
  lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 8,
        col = "gray50")
  lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 7,
        col = "gray53")
  lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 6,
        col = "gray59")
  lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 4,
        col = "gray75")
  lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 2,
        col = "gray85")
  lines(c(sum + 2e+07, genome_size - 2e+07), c(1, 1), lwd = 1,
        col = "gray90")
  if (Organism == "Human") {
    lines(c(sum + centromere_pos[24], sum + centromere_pos[24]),
          c(1.01, 0.99), col = "grey13", lwd = 2)
  }
}

PlotZygositySinle2 <- function(Table, Organism) {
  tbl = Table
  if (Organism == "Human") {
    centromere_pos = c(125, 93.3, 91, 50.4, 48.4, 61, 59.9,
                       45.6, 49, 40.2, 53.7, 35.8, 17.9, 17.6, 19, 36.6,
                       24, 17.2, 26.5, 27.5, 13.2, 14.7, 60.6, 12.5)
    centromere_pos = centromere_pos * 1e+06
    chr_size = c(248956422, 242193529, 198295559, 190214555,
                 181538259, 170805979, 159345973, 145138636, 138394717,
                 133797422, 135086622, 133275309, 114364328, 107043718,
                 101991189, 90338345, 83257441, 80373285, 58617616,
                 64444167, 46709983, 50818468, 156040895, 57227415)
    lb = c(seq(from = 150, to = 10, by = -20), seq(from = 10,
                                                   to = 150, by = 20))
    plot(c(1:24), rep(0, 24), ylim = c(-1.5e+08, 1.5e+08),
         cex.axis = 0.7, xlab = "Chromsome", ylab = "Position (Mbp)",
         xaxt = "n", yaxt = "n")
    axis(2, seq(from = 1.5e+08, to = -1.5e+08, by = -2e+07),
         labels = lb, cex.axis = 0.7, las = 1)
    mx = 24
  }
  if (Organism == "Mouse") {
    chr_size = c(195471971, 182113224, 160039680, 156508116,
                 151834684, 149736546, 145441459, 129401213, 124595110,
                 130694993, 122082543, 120129022, 120421639, 124902244,
                 104043685, 98207768, 94987271, 90702639, 61431566,
                 171031299, 91744698)
    centromere_pos = rep(2e+08, 21)
    lb = c(seq(from = 0, to = 200, by = 20))
    plot(c(1:21), rep(2e+08, 21), ylim = c(0, 2.15e+08),
         cex.axis = 0.7, xlab = "Chromsome", ylab = "Position (Mbp)",
         xaxt = "n", yaxt = "n")
    axis(2, seq(from = 2e+08, to = 0, by = -2e+07), labels = lb,
         cex.axis = 0.7, las = 1)
    mx = 21
  }
  for (i in 1:mx) {
    chr = i
    tbl2 = tbl[tbl$chr.x == chr,]
    tbl2$snp = tbl2$AF1 > 0
    tr = tbl2[tbl2$snp == T,]
    fl = tbl2[tbl2$snp == F,]
    if (dim(fl)[1] > 0) {
      for (j in 1:dim(fl)[1]) {
        x = centromere_pos[i] - fl$start[j]
        lines(c(i, i - 0.4), c(x, x), col = "blue",
              lwd = 1)
      }
    }
    if (dim(tr)[1] > 0) {
      for (j in 1:dim(tr)[1]) {
        x = centromere_pos[i] - tr$start[j]
        lines(c(i, i + 0.4), c(x, x), col = "red", cex = 0.05,
              lwd = 1)
      }
    }
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 10, col = "black")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 8, col = "gray50")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 7, col = "gray53")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 6, col = "gray59")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 4, col = "gray75")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 2, col = "gray85")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 1, col = "gray90")
    if (Organism == "Human") {
      points(i, 0, pch = 16, col = "grey13")
      if (i < 23) {
        text(i, centromere_pos[i] + 1.8e+07, i, cex = 0.8)
      }
      if (i == 23) {
        text(i, centromere_pos[i] + 1.8e+07, "X", cex = 0.8)
      }
      if (i == 24) {
        text(i, centromere_pos[i] + 1.8e+07, "Y", cex = 0.8)
      }
    }
    if (Organism == "Mouse") {
      points(i, 2e+08, pch = 16, col = "grey13")
      if (i < 20) {
        text(i, 2.12e+08, i, cex = 0.8)
      }
      if (i == 20) {
        text(i, 2.12e+08, "X", cex = 0.8)
      }
      if (i == 21) {
        text(i, 2.12e+08, "Y", cex = 0.8)
      }
    }
  }
}

PlotZygosityBlocks2 <- function(Table, Window, Max, Max2, Organism) {
  max = Max
  max2 = Max2
  window = Window
  tbl = Table
  if (Organism == "Human") {
    centromere_pos = c(125, 93.3, 91, 50.4, 48.4, 61, 59.9,
                       45.6, 49, 40.2, 53.7, 35.8, 17.9, 17.6, 19, 36.6,
                       24, 17.2, 26.5, 27.5, 13.2, 14.7, 60.6, 12.5)
    centromere_pos = centromere_pos * 1e+06
    chr_size = c(248956422, 242193529, 198295559, 190214555,
                 181538259, 170805979, 159345973, 145138636, 138394717,
                 133797422, 135086622, 133275309, 114364328, 107043718,
                 101991189, 90338345, 83257441, 80373285, 58617616,
                 64444167, 46709983, 50818468, 156040895, 57227415)
    lb = c(seq(from = 150, to = 10, by = -20), seq(from = 10,
                                                   to = 150, by = 20))
    plot(c(1:24), rep(0, 24), ylim = c(-1.5e+08, 1.5e+08),
         cex.axis = 0.7, xlab = "Chromsome", ylab = "Position (Mbp)",
         xaxt = "n", yaxt = "n")
    axis(2, seq(from = 1.5e+08, to = -1.5e+08, by = -2e+07),
         labels = lb, cex.axis = 0.7, las = 1)
    mx = 24
  }
  if (Organism == "Mouse") {
    chr_size = c(195471971, 182113224, 160039680, 156508116,
                 151834684, 149736546, 145441459, 129401213, 124595110,
                 130694993, 122082543, 120129022, 120421639, 124902244,
                 104043685, 98207768, 94987271, 90702639, 61431566,
                 171031299, 91744698)
    centromere_pos = rep(2e+08, 21)
    lb = c(seq(from = 0, to = 200, by = 20))
    plot(c(1:21), rep(2e+08, 21), ylim = c(0, 2.15e+08),
         cex.axis = 0.7, xlab = "Chromsome", ylab = "Position (Mbp)",
         xaxt = "n", yaxt = "n")
    axis(2, seq(from = 2e+08, to = 0, by = -2e+07), labels = lb,
         cex.axis = 0.7, las = 1)
    mx = 21
  }
  tbl$snp = tbl$AF1 > 0
  total_homo = sum(tbl$snp == FALSE)
  total_hetro = sum(tbl$snp == TRUE)
  pval = NULL
  rat = NULL
  for (i in 1:mx) {
    if (i < (mx - 1)) {
      chr = i
      tbl2 = tbl[tbl$chr.x == chr,]
      d = tbl2[tbl2$start < centromere_pos[i],]
      homo = sum(d$snp == FALSE)
      hetro = sum(d$snp == TRUE)
      rat = c(rat, homo / hetro)
      d = tbl2[tbl2$start > centromere_pos[i],]
      homo = sum(d$snp == FALSE)
      hetro = sum(d$snp == TRUE)
      rat = c(rat, homo / hetro)
    }
  }
  rat[is.infinite(rat) == TRUE] = total_homo / total_hetro *
    5.5
  pval = NULL
  for (i in 1:((mx - 2) * 2)) {
    np = 1
    if (is.na(rat[i]) == FALSE) {
      np = t.test(rat, mu = rat[i])$p.value
    }
    pval = cbind(pval, np)
  }
  pval = p.adjust(pval, method = "fdr")
  pval = log(pval, 10) * (-1)
  rat[is.infinite(rat) == TRUE] = total_homo / total_hetro
  rat[is.na(rat) == TRUE] = total_homo / total_hetro
  y = pval > 3 & rat > (total_homo / total_hetro * 5)
  loc = 1
  for (i in 1:(mx - 2)) {
    for (j in 1:2) {
      if (y[loc] == FALSE) {
        col = "white"
      }
      if (y[loc] == TRUE) {
        col = rgb(249, 244, 93, maxColorValue = 255)
      }
      if (j == 1) {
        rect(xleft = i - 0.45, ybottom = +centromere_pos[i],
             xright = i + 0.45, 0, col = col, border = NA)
      }
      if (j == 2) {
        rect(xleft = i - 0.45, ybottom = +centromere_pos[i] -
          chr_size[i], xright = i + 0.45, 0, col = col,
             border = NA)
      }
      loc = loc + 1
    }
  }
  for (i in 1:mx) {
    chr = i
    tbl2 = tbl[tbl$chr.x == chr,]
    tr = tbl2[tbl2$snp == T,]
    fl = tbl2[tbl2$snp == F,]
    win = ceiling(chr_size[i] / window)
    pos = centromere_pos[i]
    library(gplots)
    cl = colorpanel(max, "grey", "red")
    cl2 = colorpanel(max2, "grey", "blue")
    for (j in 1:win) {
      top = j * window
      bottom = (j - 1) * window
      num_fl = sum(fl$start >= bottom & fl$start <= top)
      num_tr = sum(tr$start >= bottom & tr$start <= top)
      if (num_tr > max) {
        num_tr = max
      }
      if (num_fl > max2) {
        num_fl = max2
      }
      rect(i, pos, i + 0.4, pos - window, col = cl[num_tr],
           border = NA)
      rect(i, pos, i - 0.4, pos - window, col = cl2[num_fl],
           border = NA)
      pos = pos - window
    }
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 10, col = "black")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 8, col = "gray50")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 7, col = "gray53")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 6, col = "gray59")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 4, col = "gray75")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 2, col = "gray85")
    lines(c(i, i), c(centromere_pos[i], centromere_pos[i] -
      chr_size[i]), lwd = 1, col = "gray90")
    if (Organism == "Human") {
      points(i, 0, pch = 16, col = "grey13")
      if (i < 23) {
        text(i, centromere_pos[i] + 1.8e+07, i, cex = 0.8)
      }
      if (i == 23) {
        text(i, centromere_pos[i] + 1.8e+07, "X", cex = 0.8)
      }
      if (i == 24) {
        text(i, centromere_pos[i] + 1.8e+07, "Y", cex = 0.8)
      }
    }
    if (Organism == "Mouse") {
      points(i, 2e+08, pch = 16, col = "grey13")
      if (i < 20) {
        text(i, 2.12e+08, i, cex = 0.8)
      }
      if (i == 20) {
        text(i, 2.12e+08, "X", cex = 0.8)
      }
      if (i == 21) {
        text(i, 2.12e+08, "Y", cex = 0.8)
      }
    }
  }
}

MajorMinorCalc<-function(Table,minDP,maxDP,minAF){

  Table[is.na(Table)] = 0
  newTable = Table[Table$DP >= minDP,]
  newTable = newTable[newTable$DP <= maxDP,]
  AF1 = newTable$AD1/newTable$DP
  AF2 = newTable$AD2/newTable$DP
  newTable = data.frame("chr" = newTable$chr, "position" = newTable$position, "AF1" = AF1, "AF2" = AF2)
  frequncyTable = newTable[newTable$AF1 >= minAF,]
  frequncyTable = frequncyTable[frequncyTable$AF2 >= minAF,]
  orderedTable = Sort_major_minor(data=frequncyTable, col1=3, col2=4)
  MajorMinor = orderedTable$AF1/orderedTable$AF2
  orderedTable["MajorMinor"] = MajorMinor
  return(orderedTable)
}

Sort_major_minor<-function(data, col1, col2){
  for (i in 1:dim(data)[1]) {
    if (data[i,col1] < data[i,col2]) {
      save = data[i,col2]
      data[i,col2] = data[i,col1]
      data[i,col1] = save
    }
  }
  return(data)
}


print("Begin eSNPKaryotyping analysis")
if (organism == "homo_sapiens") {
  organism = "Human"
}
if (organism == "mus_musculus") {
  organism = "Mouse"
}

Edit_dbSNP_Files2(Directory = dbSNP_dir, File_Name = "/Common_SNPs_142_chr", Organism = organism)
table = EditVCF2(Directory = vcf_edited_dir, Organism = organism)
table$chr = as.numeric(table$chr)
table = table[order(table$chr, table$position),]
table = table[table$chr > 0,]
table2 = MajorMinorCalc(Table = table,
                        minDP = 20,
                        maxDP = 1000000,
                        minAF = 0.2)
pdf(paste(out_file_prefix, "Genome_Zygosity.pdf", sep = "_"), width = 5, height = 5)

#Plot Allelic ratio along the genome for duplication detection
PlotGenome2(orderedTable = table2,
            Window = 151,
            Ylim = 3,
            PValue = TRUE,
            Organism = organism)
dev.off()
#Intersect the observed SNPs with the common SNPs table from dbSNPs, Creates file with the LOH data called
tbl = DeletionTable2(Directory = vcf_edited_dir,
                     Table = table2,
                     dbSNP_Data_Directory = dbSNP_dir,
                     Genome_Fa_dict = ref_picard_dic,
                     Organism = organism,
                     bam_Directory = bams_path)
deletions = read.delim(paste(vcf_edited_dir, "Deletions.txt", sep = "/"), sep = "\t")
pdf(paste(out_file_prefix, "Zygosity_Sinle.pdf", sep = "_"), width = 5, height = 5)
# Plot each SNP, without any summarization
PlotZygositySinle2(Table = deletions, Organism = organism)
dev.off()
pdf(paste(out_file_prefix, "Zygosity_Blocks.pdf", sep = "_"), width = 5, height = 5)
# Plot blocks of heterozygous and homozygous SNPs
PlotZygosityBlocks2(Table = deletions,
                    Window = 1500000,
                    Max = 6,
                    Max2 = 60,
                    Organism = organism)
dev.off()

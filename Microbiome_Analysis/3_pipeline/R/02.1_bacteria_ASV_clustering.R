#!/usr/bin/env Rscript

#clear the object from memory
rm(list=ls())

#### CHANGE THIS PATHS ####
source <- "your_path_to_file_source"
taxa_database <- "your_path_to_database" #we used silva_nr_v132_train_set.fa.gz
##########################


# install required libraries before running this script


# load libraries
library(dada2)
library(gtools)
library(grid)
library(gridExtra)
library(xlsx)
library(tidyr)
library(vegan)

#set working directory to source file location
setwd(source)

#generate output files
output <- T
#assign taxonomy
taxonomy <- T

#list of bacteria runs
runs <- list.files("../cmd/demultiplexed/bacteria")

#folder to store interim results
temp <- "bacteria_interim/"
dir.create(temp)


# ------------------------------------------------------------------------
# Loop over all bacteria runs
# ------------------------------------------------------------------------
for (run in runs) {
  print(paste(run, "started"))
  
  # ------------------------------------------------------------------------
  # get the samples
  # ------------------------------------------------------------------------
  
  #Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
  path <- paste("../cmd/demultiplexed/bacteria/", run, sep="")
  fnFs <- sort(list.files(path, pattern="r1.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="r2.fastq", full.names = TRUE))
  
  #Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
  sample.names <- sapply(strsplit(basename(fnFs), ".", fixed=T), `[`, 1)
  
  #Place filtered files in filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  #give them a unique name
  assign(paste(run, "_fnFs", sep=""), fnFs)
  assign(paste(run, "_fnRs", sep=""), fnRs)
  assign(paste(run, "_filtFs", sep=""), filtFs)
  assign(paste(run, "_filtRs", sep=""), filtRs)
  
  # ------------------------------------------------------------------------
  #Quality filtering and trimming
  # ------------------------------------------------------------------------
  
  #filter
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,170),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)
  
  #unique name and save
  assign(paste(run, "_out", sep=""), out)
  saveRDS(out, paste(temp, run, "_out.RDS", sep=""))
  
  #Avoids to raise an error message if no reads from a sample pass the filter
  missing_fltFs <- c()
  missing_fltRs <- c()
  for (file in filtFs) {if(!file.exists(file)) {missing_fltFs <- c(missing_fltFs, file)}} #collect names of the missing files
  for (file in filtRs) {if(!file.exists(file)) {missing_fltRs <- c(missing_fltRs, file)}}
  filtFs <- setdiff(filtFs, missing_fltFs)
  filtRs <- setdiff(filtRs, missing_fltRs)

  # ------------------------------------------------------------------------
  # Learn the error rate
  # ------------------------------------------------------------------------
  
  #error rate
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  
  #unique name and save
  assign(paste(run, "_errF", sep=""), errF)
  assign(paste(run, "_errR", sep=""), errR)
  saveRDS(errF, paste(temp, run, "_errF.RDS", sep=""))
  saveRDS(errR, paste(temp, run, "_errR.RDS", sep=""))
  
  # ------------------------------------------------------------------------
  # Dereplicate the files
  # ------------------------------------------------------------------------
  
  #dereplicate
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  
  #unique name and save
  assign(paste(run, "_derepFs", sep=""), derepFs)
  assign(paste(run, "_derepRs", sep=""), derepRs)
  saveRDS(derepFs, paste(temp, run, "_derepFs.RDS", sep=""))
  saveRDS(derepRs, paste(temp, run, "_derepRs.RDS", sep=""))
  
  # ------------------------------------------------------------------------
  #Sample inference algorithm
  # ------------------------------------------------------------------------
  
  # inference algortihm
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
  #unique name and save
  assign(paste(run, "_dadaFs", sep=""), dadaFs)
  assign(paste(run, "_dadaRs", sep=""), dadaRs)
  saveRDS(dadaFs, paste(temp, run, "_dadaFs.RDS", sep=""))
  saveRDS(dadaRs, paste(temp, run, "_dadaRs.RDS", sep=""))
  
}

# ------------------------------------------------------------------------
#combine all bacteria runs
# ------------------------------------------------------------------------

#combine and save
fnFs <- vector()
for(i in runs){fnFs <- c(fnFs, get(paste(i, "_fnFs", sep="")))}
fnRs <- vector()
for(i in runs){fnRs <- c(fnRs, get(paste(i, "_fnRs", sep="")))}
saveRDS(fnFs, paste0(temp, "fnFs.RDS"))
saveRDS(fnRs, paste0(temp, "fnRs.RDS"))

filtFs <- vector()
for(i in runs){filtFs <- c(filtFs, get(paste(i, "_filtFs", sep="")))}
filtRs <- vector()
for(i in runs){filtRs <- c(filtRs, get(paste(i, "_filtRs", sep="")))}
saveRDS(filtFs, paste0(temp, "filtFs.RDS"))
saveRDS(filtRs, paste0(temp, "filtRs.RDS"))

out <- get(paste(runs[1], "_out", sep=""))
if(length(runs)>1){for(i in runs[2:length(runs)]){out <- rbind(out, get(paste(i, "_out", sep="")))}}
saveRDS(out, paste0(temp, "out.RDS"))

derepFs <- vector()
for(i in runs){derepFs <- c(derepFs, get(paste(i, "_derepFs", sep="")))}
derepRs <- vector()
for(i in runs){derepRs <- c(derepRs, get(paste(i, "_derepRs", sep="")))}
saveRDS(derepFs, paste0(temp, "derepFs.RDS"))
saveRDS(derepRs, paste0(temp, "derepRs.RDS"))

dadaFs <- vector()
for(i in runs){dadaFs <- c(dadaFs, get(paste(i, "_dadaFs", sep="")))}
dadaRs <- vector()
for(i in runs){dadaRs <- c(dadaRs, get(paste(i, "_dadaRs", sep="")))}
saveRDS(dadaFs, paste0(temp, "dadaFs.RDS"))
saveRDS(dadaRs, paste0(temp, "dadaRs.RDS"))

#names from all bacterial runs
sample.names <- sapply(strsplit(basename(fnFs), ".", fixed=T), `[`, 1)

# ------------------------------------------------------------------------
# Merge the forward and reverse reads
# ------------------------------------------------------------------------

#merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers, paste0(temp, "mergers.RDS"))

# ------------------------------------------------------------------------
# ASV table
# ------------------------------------------------------------------------

#create table and remove bimera
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab, paste0(temp, "seqtab.RDS"))
saveRDS(seqtab.nochim, paste0(temp, "seqtab_nochim.RDS"))

# ------------------------------------------------------------------------
# track table
# ------------------------------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
saveRDS(track, paste0(temp, "track.RDS"))

# ------------------------------------------------------------------------
# Assign taxanomy
# ------------------------------------------------------------------------
if(taxonomy) {
  taxa <- assignTaxonomy(seqtab.nochim, taxa_database, multithread=TRUE)
  saveRDS(taxa, paste0(temp, "taxa.RDS"))
}

# ------------------------------------------------------------------------
# Header names
# ------------------------------------------------------------------------

#change headers from the long sequence to simple ASV numbers
seq <- colnames(seqtab.nochim)
ASV <- c()
for (i in 1:length(seq)) {ASV[i] <- paste("ASV",i , sep="") }
ASV_seq <- cbind(ASV, seq)
ASV_seq[,1] <- paste(">", ASV_seq[,1], sep="")

#replace sequence headers by ASV numbers in seqtab.nochim
colnames(seqtab.nochim) <- ASV

if(taxonomy) {
  #replace the sequence headers by the correct ASV number in taxa
  for (i in 1:nrow(taxa)){
    for (j in 1:nrow(ASV_seq)){
      if(rownames(taxa)[i] == ASV_seq[j,2]){
        rownames(taxa)[i] <- ASV_seq[j,1]
      }
    }
  }
}

# ------------------------------------------------------------------------
# Generate output files
# ------------------------------------------------------------------------
if(output){
  
  setwd(source)
  dir.create("../../4_output/bacteria_ASV")
  path_out <- c("../../4_output/bacteria_ASV/")
  
  #quality profiles unfiltered and untrimmed
  pdf(paste(path_out,"bacteria_reads_quality_unfilt_untrim.pdf",sep=""))
    for (i in 1:length(sample.names)) {
      try(figure <- plotQualityProfile(c(fnFs[i],fnRs[i])))
      try(print(figure))
      try(rm(figure))
    }
  dev.off()
  
  #quality profiles filtered and trimmed
  pdf(paste(path_out,"bacteria_reads_quality_filt_trim.pdf",sep=""))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(c(filtFs[i],filtRs[i])))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()
  
  #error rate learning
  runs_errF <- c (paste(runs, "_errF", sep=""))
  runs_errR <- c (paste(runs, "_errR", sep=""))
  
  pdf(paste(path_out,"bacteria_error_rate.pdf",sep=""))
    for (i in 1:length(runs)) {
      plot(plotErrors(get(runs_errF[i]), nominalQ=TRUE))
      grid.text(paste(runs[i], "forward"),hjust=-1.2, vjust = -27.5, rot = 90)
      plot(plotErrors(get(runs_errR[i]), nominalQ=TRUE))
      grid.text(paste(runs[i], "reverse"),hjust=-1.2, vjust = -27.5, rot = 90)
      }
  dev.off()
  
  #rarefaction plots
  pdf(paste(path_out,"bacteria_rarefaction.pdf",sep=""))
    for (run in runs) {
      DAT <- seqtab.nochim[grep(run, rownames(seqtab.nochim)),]
      rarecurve(DAT, step = 20, label=F, main = paste("Rarefaction", run), ylab="ASVs", xlab="Sequencing depth")
    }
  dev.off()
  
  #tables
  write.table(seqtab.nochim, paste(path_out, "bacteria_ASV.tab", sep=""), sep="\t")
  if (taxonomy) {write.table(taxa, paste(path_out, "bacteria_taxa.tab", sep=""), sep="\t")}
  write.table(ASV_seq, paste(path_out, "bacteria_sequences.fasta", sep=""), quote=F, row.names = FALSE, col.names = FALSE)
  write.xlsx(track, paste(path_out, "bacteria_track.xlsx", sep=""))
}



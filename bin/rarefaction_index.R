#!/bin/Rscript
# Generates an index of jobs to distribute to bootstrap rarefaction curves.
# argument 1: character vector of populations <- DOES NOTHING FOR NOW
# argument 2: integer, maximum number of samples per group size
# argument 3: integer, number of processors to use
# argument 4: path to sample_groups
# argument 5: path to *.traw file

library(readr)
library(tidyr)
library(dplyr)
library(gplots)
library(foreach)
library(doMC)
args <- commandArgs(trailingOnly = TRUE)

## For testing
#args <- 1
#args[1] <- "jacksoni,sumatrae,altaica,tigris"
#args[2] <- 100
#args[3] <- 8
#args[4] <- "../../sample_groups.tsv"
#args[5] <- "fb-170208-1-snp-qq30-gq30-hwe-ss0.001-ldp.traw"
#setwd("/Volumes/bio-dap15.stanford.edu/zstor/2016-tiger-wgs/basic-stats/bs-test/outputs/rarefaction")
###

#pops <- strsplit(args[1], ',')[[1]]
registerDoMC(cores = args[3])

genotypes <- read_tsv(args[5], col_names=TRUE)
sample_groups <- read_tsv(args[4], col_names=c("ID", "pop"))
#sample_groups <- sample_groups %>% filter(pop %in%  pops)
colnames(genotypes) <- gsub("_.*", "", colnames(genotypes))
genotypes_cols <- colnames(genotypes)
sample_groups <- sample_groups[sample_groups$ID %in% genotypes_cols, ]
sample_groups$idcol <- match(sample_groups$ID, genotypes_cols)

rarefaction_index <- foreach(group.i = unique(sample_groups$pop), .combine = c) %do% {
  sg.i <- sample_groups %>% filter(pop == group.i)
  foreach(sample_size.j = 1:nrow(sg.i), .combine = c) %dopar% {
    if(choose(nrow(sg.i), sample_size.j) > args[2]){
      sample_combos <- matrix( , sample_size.j, as.integer(args[2]))
      for(i in 1:args[2]) {sample_combos[ ,i] <- sample(sg.i$idcol, sample_size.j)}
    } else {
      sample_combos <- combn(sg.i$idcol, sample_size.j)
    }
    n_combos <- dim(sample_combos)[2]
    foreach(combo.k = 1:n_combos) %do% {
      #genotypes.k <- genotypes[ ,sample_combos[ ,combo.k]]
      #n_poly <- sum(apply(genotypes.k, 1, function(x) any(x == 1, na.rm=TRUE) || (any(x == 0, na.rm=TRUE) & any(x == 2, na.rm=TRUE))))
      list(pop = group.i, n_ind = sample_size.j, itt = combo.k, samples = c(sample_combos[ ,combo.k]))
    }
  }
}

save(rarefaction_index, file = paste(args[6], ".rfi.Rdata", sep = ""))
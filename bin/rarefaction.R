#!/bin/Rscript
# Saves an R data object with a table 'frq_data'. 'frq_data' has a column for locus, and the minor allele frequency for each population
# argument 1: integer, node number
# argument 2: integer, total number of nodes
# argument 3: integer, number of processors to use
# argument 4: path to *.traw file
# argument 5: path to rfi-*.Rdata file
# argument 6: prefix for outfile

library(readr)
library(tidyr)
library(dplyr)
library(gplots)
library(foreach)
library(doMC)
args <- commandArgs(trailingOnly = TRUE)

## For testing
#args <- 1
#args[1] <- 1
#args[2] <- 8
#args[3] <- 8
#args[4] <- "/zstor/2016-tiger-wgs/basic-stats/bs-test/outputs/plink/fb-170208-1-snp-qq30-gq30-hwe-ss0.001.traw"
#args[5] <- "/zstor/2016-tiger-wgs/basic-stats/bs-test/outputs/rarefaction/rfi.fb-170208-1-snp-qq30-gq30-hwe-ss0.001.Rdata"
#args
###

registerDoMC(cores = args[3])

load(args[5])
genotypes <- read_tsv(args[4], col_names=TRUE)
colnames(genotypes) <- gsub("_.*", "", colnames(genotypes))

# Get jobs for this node
n_jobs <- length(rarefaction_index)
jobs <- data_frame(index=1:n_jobs, job = cut(1:n_jobs, as.integer(args[2]), labels=FALSE))
node.jobs <- jobs[jobs$job == as.integer(args[1]),]

rarefaction <- foreach(nj.i = node.jobs$index, .combine = rbind) %dopar% {
  task.i <- rarefaction_index[[nj.i]]
  genotypes.i <- genotypes[ ,task.i$samples]
  n_poly <- sum(apply(genotypes.i, 1, function(x) any(x == 1, na.rm=TRUE) || (any(x == 0, na.rm=TRUE) & any(x == 2, na.rm=TRUE))))
  data_frame(pop = task.i$pop, n_ind = task.i$n_ind, itt = task.i$itt, n_poly = n_poly, inds = paste(colnames(genotypes)[task.i$samples], collapse = ','))
}

save(rarefaction, file = paste(args[6],".rf-", args[1], ".Rdata", sep = ''))
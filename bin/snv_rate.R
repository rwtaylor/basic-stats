#!/bin/Rscript
# Calculates snv rates for each individual in a plink traw file
# argument 1: integer, number of processors to use
# argument 2: path to *.traw file

library(readr)
library(tidyr)
library(dplyr)
library(foreach)
library(doMC)
args <- commandArgs(trailingOnly = TRUE)

## For testing
#args <- 1
#args[1] <- 8
#args[2] <- "fb-170208-1-snp-qq30-gq30-hwe-ss0.001-ldp.traw"
###
registerDoMC(cores = args[1])

genotypes <- read_tsv(args[2], col_names=TRUE)
colnames(genotypes) <- gsub("_.*", "", colnames(genotypes))

het_loci_counts <- apply(genotypes[ ,7:ncol(genotypes)], 2, function(x) sum(x == 1, na.rm = TRUE))



#!/bin/Rscript
# Saves an R data object with a table 'frq_data'. 'frq_data' has a column for locus, and the minor allele frequency for each population
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
library(gplots)
args <- commandArgs(trailingOnly = TRUE)

## For testing
#args <- 1
#args[1] <- "jacksoni,sumatrae,altaica,tigris"
#args[2] <- 100
#args[3] <- 8
#args[4] <- "../../sample_groups.tsv"
#args[5] <- "fb-170208-1-snp-qq30-gq30-hwe-ss0.001-ldp.traw"
#args[6] <- "fb-170208-1-snp-qq30-gq30-hwe-ss0.001-ldp"
###

pops <- strsplit(args[1], ',')[[1]]
registerDoMC(cores = args[3])

genotypes <- read_tsv(args[5], col_names=TRUE)
sample_groups <- read_tsv(args[4], col_names=c("ID", "pop"))
#sample_groups <- sample_groups %>% filter(pop %in%  pops)
colnames(genotypes) <- gsub("_.*", "", colnames(genotypes))
genotypes_cols <- colnames(genotypes)
sample_groups <- sample_groups[sample_groups$ID %in% genotypes_cols, ]
sample_groups$idcol <- match(sample_groups$ID, genotypes_cols)

pop_polymorphic_sites <- foreach(group.i = unique(sample_groups$pop), .combine = cbind) %do% {
  sg.i <- sample_groups %>% filter(pop == group.i)
  genotypes.i <- genotypes[ ,sg.i$idcol]
  out <- data_frame(pop = apply(genotypes.i, 1, function(x) any(x == 1, na.rm=TRUE) || (any(x == 0, na.rm=TRUE) & any(x == 2, na.rm=TRUE))))
  names(out) <- group.i
  out
}

venn_diagram <- venn(pop_polymorphic_sites[ , pops], show.plot = FALSE)

pdf(file = paste(args[6], ".venn.pdf", sep = ''))
plot(venn_diagram)
dev.off()

png(file = paste(args[6], ".venn.png", sep = ''))
plot(venn_diagram)
dev.off()

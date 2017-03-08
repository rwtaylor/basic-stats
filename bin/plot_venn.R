#!/bin/Rscript
# Plots venn diagram
# argument 1: cpus
# argument 2: comma separated list of populations
# argument 3: path to sample_groups
# argument 4: path to *.traw file
# argument 5: output prefix

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
#args[1] <- 8
#args[2] <- "jacksoni,sumatrae,altaica,tigris"
#args[3] <- "../../sample_groups.tsv"
#args[4] <- "fb-170208-1-snp-qq30-gq30-hwe-ss0.001-ldp.traw"
#args[5] <- "fb-170208-1-snp-qq30-gq30-hwe-ss0.001-ldp"
###

registerDoMC(cores = args[1])
pops <- strsplit(args[2], ',')[[1]]
genotypes <- read_tsv(args[4], col_names=TRUE)
sample_groups <- read_tsv(args[3], col_names=c("ID", "pop"))

sample_groups <- sample_groups %>% filter(pop %in% pops)

#sample_groups <- sample_groups %>% filter(pop %in%  pops)
colnames(genotypes) <- gsub("_.*", "", colnames(genotypes))
genotypes_cols <- colnames(genotypes)
sample_groups <- sample_groups[sample_groups$ID %in% genotypes_cols, ]
sample_groups$idcol <- match(sample_groups$ID, genotypes_cols)

pop_polymorphic_sites <- foreach(group.i = unique(sample_groups$pop), .combine = cbind) %dopar% {
  sg.i <- sample_groups %>% filter(pop == group.i)
  genotypes.i <- genotypes[ ,sg.i$idcol]
  out <- data_frame(pop = apply(genotypes.i, 1, function(x) any(x == 1, na.rm=TRUE) || (any(x == 0, na.rm=TRUE) & any(x == 2, na.rm=TRUE))))
  names(out) <- group.i
  out
}

venn_diagram <- venn(pop_polymorphic_sites[ , pops], show.plot = FALSE)

pdf(file = paste(args[5], ".venn.pdf", sep = ''))
plot(venn_diagram)
dev.off()

png(file = paste(args[5], ".venn.png", sep = ''))
plot(venn_diagram)
dev.off()

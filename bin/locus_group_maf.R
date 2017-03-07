#!/bin/Rscript
# Saves an R data object with a table 'frq_data'. 'frq_data' has a column for locus, and the minor allele frequency for each population
# argument 1: character vector of populations
# argument 2: character vector of frq files for each population
# argument 3: character string of output file name

library(readr)
library(tidyr)
library(dplyr)
library(gplots)
args <- commandArgs(trailingOnly = TRUE)

## For testing
#args <- 1
#args[1] <- "jacksoni,sumatrae,altaica,tigris"
#args[2] <- "fb-test-snp-q20-gq20-jacksoni.frq,fb-test-snp-q20-gq20-sumatrae.frq,fb-test-snp-q20-gq20-altaica.frq,fb-test-snp-q20-gq20-tigris.frq"
#args[3] <- "fb-test-snp-q20-gq20.lgmaf"

pops <- strsplit(args[1], ',')[[1]]
frq_files <- strsplit(args[2], ',')[[1]]

print(pops)
print(frq_files)

import <- read_tsv(frq_files[1], skip = 1, col_names=c("chrom", "pos", "n_alleles", "n_chr", "af1", "af2", "af3", "af4"))
frq_data <- import %>% separate(af1, c("a1", "f1"), ":") %>%
    separate(af2, c("a2", "f2"), ":") %>%
    separate(af3, c("a3", "f3"), ":") %>%
    separate(af4, c("a4", "f4"), ":") %>%
    mutate(f1 = as.numeric(f1), f2 = as.numeric(f2), f3 = as.numeric(f3), f4 = as.numeric(f4)) %>%
    rowwise() %>%
    mutate(maf = min(c(f1,f2,f3,f4), na.rm = TRUE)) %>%
    mutate(maf_bi = ifelse((is.na(f4) & is.na(f3)), min(c(f1,f2)), NA)) %>%
    select(chrom, pos, maf)
names(frq_data) <- c("chrom", "pos", pops[1])


for(i in 2:length(pops)){
    import <- read_tsv(frq_files[i], skip = 1, col_names=c("chrom", "pos", "n_alleles", "n_chr", "af1", "af2", "af3", "af4"))
    frq_data.i <- import %>% separate(af1, c("a1", "f1"), ":") %>%
        separate(af2, c("a2", "f2"), ":") %>%
        separate(af3, c("a3", "f3"), ":") %>%
        separate(af4, c("a4", "f4"), ":") %>%
        mutate(f1 = as.numeric(f1), f2 = as.numeric(f2), f3 = as.numeric(f3), f4 = as.numeric(f4)) %>%
        rowwise() %>%
        mutate(maf = min(c(f1,f2,f3,f4), na.rm = TRUE)) %>%
        mutate(maf_bi = ifelse((is.na(f4) & is.na(f3)), min(c(f1,f2)), NA)) %>%
        select(chrom, pos, maf)
    names(frq_data.i) <- c("chrom", "pos", pops[i])
    frq_data <- left_join(frq_data, frq_data.i)
}

print("Finished joining")

print("Writing tsv")
write.table(frq_data, sep='\t', col.names = NA, file = paste(args[3], ".lgmaf", sep = ''), quote = FALSE)
print("Finished writing tsv, now writing Rdata")
save(frq_data, file = paste(args[3], ".Rdata", sep = ''))
print("Finished writing Rdata")


#!/usr/bin/Rscript

library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)


import <- read_tsv(args[2], skip=1, col_names = c("chrom", "pos", "n_alleles", "n_chr", "af1", "af2", "af3", "af4"))

frq <- import %>% separate(af1, c("a1", "f1"), ":") %>%
    separate(af2, c("a2", "f2"), ":") %>%
    separate(af3, c("a3", "f3"), ":") %>%
    separate(af4, c("a4", "f4"), ":") %>%
    mutate(f1 = as.numeric(f1), f2 = as.numeric(f2), f3 = as.numeric(f3), f4 = as.numeric(f4)) %>%
    rowwise() %>%
    mutate(maf = min(c(f1,f2,f3,f4), na.rm = TRUE)) %>%
    mutate(maf_bi = ifelse((is.na(f4) & is.na(f3)), min(c(f1,f2)), NA))

p1 <- ggplot(frq, aes(x = maf)) + geom_histogram(binwidth=0.025) + ggtitle("Minor Allele Frequency Histogram (all loci)")
ggsave(p1, file = paste(args[1], "-maf-hist.pdf", sep = ''))
ggsave(p1, file = paste(args[1], "-maf-hist.png", sep = ''))

p2 <- ggplot(frq, aes(x = maf_bi)) + geom_histogram(binwidth=0.025) + ggtitle("Minor Allele Frequency Histogram (bi-allelic loci)")
ggsave(p2, file = paste(args[1], "-bimaf-hist.pdf", sep = ''))
ggsave(p2, file = paste(args[1], "-bimaf-hist.png", sep = ''))

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
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

## For testing
#args <- 1
#args[1] <- 2
#args[2] <- "fb-170208-1-snp-qq30-gq30-hwe-ss0.001"
#setwd("/Volumes/bio-dap15.stanford.edu/zstor/2016-tiger-wgs/basic-stats/bs-test/outputs/rarefaction")
###

load(paste(args[2], ".rf-1.Rdata", sep = ''))

rf_concat <- rarefaction

for(i in 2:args[1]) {
  load(paste(args[2], ".rf-",i,".Rdata", sep = ''))
  rf_concat <- rbind(rf_concat, rarefaction)
}

rf_means = rf_concat %>% group_by(pop, n_ind) %>% summarize(mean_snp = mean(n_poly))
save(rf_means, rf_concat, file = paste(args[2], ".rf_concat.Rdata", sep = ''))

p <- ggplot(rf_concat, aes(n_ind, n_poly, color = pop)) + geom_point() + geom_line(data = rf_means, aes(n_ind, mean_snp), color = "black") + facet_wrap( ~ pop) + xlab("N Individuals") + ylab("N SNPs") + scale_color_discrete(name = "Population")

ggsave(p, file = paste(args[2],".png", sep = ''))
ggsave(p, file = paste(args[2],".pdf", sep = ''))

p2 <- ggplot(rf_means, aes(n_ind, mean_snp*1000, color = pop)) + geom_line(size = 1) + xlab("N Individuals") + ylab("N SNPs") + scale_color_discrete(name = "Population")

ggsave(p2, file = paste(args[2],".lm.png", sep = ''))
ggsave(p2, file = paste(args[2],".lm.pdf", sep = ''))

p3 <- ggplot(rf_means %>% filter(pop != "zoo_tigris"), aes(n_ind, log(mean_snp*1000), color = pop)) + geom_line(size = 1) + xlab("N Individuals") + ylab("N SNPs") + scale_color_discrete(name = "Population")

ggsave(p3, file = paste(args[2],".nozoo.mean.png", sep = ''))
ggsave(p3, file = paste(args[2],".nozoo.mean.pdf", sep = ''))

p4 <- ggplot(rf_concat %>% filter(pop %in% c("tigris", "sumatrae", "jacksoni", "altaica", "zoo_tigris")), aes(n_ind, n_poly*1000, group = pop, color = pop)) + geom_smooth(method = "lm", formula = ' y ~ log(x)') + xlab("N Individuals") + ylab("N SNPs") + scale_color_discrete(name = "Population")

ggsave(p4, file = paste(args[2],".lm.png", sep = ''))
ggsave(p4, file = paste(args[2],".lm.pdf", sep = ''))

p5 <- ggplot(rf_concat %>% filter(pop %in% c("tigris", "sumatrae", "jacksoni", "altaica")), aes(n_ind, n_poly*1000, group = pop, color = pop)) + geom_smooth(method = "lm", formula = ' y ~ log(x)') + xlab("N Individuals") + ylab("N SNPs") + scale_color_discrete(name = "Population")

ggsave(p5, file = paste(args[2],".nozoo.lm.png", sep = ''))
ggsave(p5, file = paste(args[2],".nozoo.lm.pdf", sep = ''))

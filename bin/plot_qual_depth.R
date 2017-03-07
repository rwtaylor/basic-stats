#!/usr/bin/Rscript

library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("q20-gq20", "fb-test-snp-q20-gq20-all.lqual", "fb-test-snp-q20-gq20-all.ldepth.mean")

colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")

import <- read_tsv(args[2], skip=1, col_names = c("chrom", "pos", "qual"))

qual <- import %>% arrange(-qual) %>% mutate(num = 1:n())

p1 <- ggplot(qual, aes(x = qual, y = num)) + geom_line() + xlab("QUAL score") + ylab("Cumulative loci")
ggsave(p1, file = paste(args[1], "qcc.pdf", sep = '-'))
ggsave(p1, file = paste(args[1], "qcc.png", sep = '-'))

p2 <- ggplot(qual, aes(x = qual, y = num / max(num))) + geom_line() + xlab("QUAL score") + ylab("Proportion of loci")
ggsave(p2, file = paste(args[1], "qcp.pdf", sep = '-'))
ggsave(p2, file = paste(args[1], "qcp.png", sep = '-'))

p3 <- ggplot(qual, aes(x = qual)) + geom_histogram() + xlab("QUAL score") + ylab("loci") + ggtitle("QUAL histogram")
ggsave(p3, file = paste(args[1], "qual-hist.pdf", sep = '-'))
ggsave(p3, file = paste(args[1], "qual-hist.png", sep = '-'))

depth <- read_tsv(args[3], skip = 1, col_names = c("chrom", "pos", "mean_depth", "var_depth"))

qual <- left_join(qual, depth)

qual_bins <- seq(0,max(qual$qual), by = 10)
depth_bins <- c(5,10,15,20,30)

qual_depth <- qual %>% filter(mean_depth >= 0) %>% mutate(qbin = cut(qual, breaks = qual_bins, labels = qual_bins[-1])) %>% group_by(qbin) %>% summarize(n = n()) %>% mutate(count = cumsum(n), depth = 0, qbin = as.numeric(qbin))

for(i in depth_bins){
    qual_depth.i <- qual %>% filter(mean_depth >= i) %>% mutate(qbin = cut(qual, breaks = qual_bins, labels = qual_bins[-1])) %>% group_by(qbin) %>% summarize(n = n()) %>% mutate(count = cumsum(n), depth = i, qbin = as.numeric(qbin))
    qual_depth <- rbind(qual_depth, qual_depth.i)
}

p4 <- ggplot(data = qual_depth, aes(x = qbin, y = count, color = as.factor(depth))) + geom_line() + xlab("QUAL") + ylab("Cumulative loci") + scale_color_discrete(name = "Min mean depth")

ggsave(p4, file = paste(args[1], "qdcc.pdf", sep = '-'))
ggsave(p4, file = paste(args[1], "qdcc.png", sep = '-'))

p5 <- ggplot(data = qual_depth, aes(x = qbin, y = count / max(count), color = as.factor(depth))) + geom_line() + xlab("QUAL") + ylab("Proportion of loci") + scale_color_discrete(name = "Min mean depth")

ggsave(p5, file = paste(args[1], "qdcp.pdf", sep = '-'))
ggsave(p5, file = paste(args[1], "qdcp.png", sep = '-'))


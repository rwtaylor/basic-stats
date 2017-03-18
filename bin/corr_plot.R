#!/usr/bin/Rscript

library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

ids <- read_tsv(file = "../outputs/plink_stats/fb-170208-snp-q30-gq30-hwe-ss0.001-ldp.rel.id", col_names=c("ID", "ID2"))

dat <- read_tsv(file = "../outputs/plink_stats/fb-170208-snp-q30-gq30-hwe-ss0.001-ldp.rel", col_names=ids$ID)

col3 <- colorRampPalette(c("red", "white", "blue"))

m <- as.matrix(dat)
m[upper.tri(m)] <- t(m)[upper.tri(m)]

m2 <- m
diag(m2) <- 0
corrplot(m2[-34,-34], is.corr=FALSE, diag = FALSE, method = "color", order = 'hclust', addrect = 4, rect.col = "black", tl.col="black")
#corrplot(m, is.corr=FALSE, diag = TRUE, method = "square", order = 'hclust')

#corrplot.mixed(m[-34,-34], is.corr=FALSE, lower = "number" ,upper = "color")

rel <- read_tsv(file = "../outputs/plink_stats/fb-170208-snp-q30-gq30-hwe-ss0.001-ldp.grm", col_names = c("ID1", "ID2", "nsnps", "rel"))
inds <- read_tsv(file = "../outputs/plink_stats/fb-170208-snp-q30-gq30-hwe-ss0.001-ldp.grm.id", col_names = c("Ind", "ID"))
inds$ID <- 1:nrow(inds)


rel <- left_join(rel, inds %>% select(Ind1 = Ind, ID1 = ID))
rel <- left_join(rel, inds %>% select(Ind2 = Ind, ID2 = ID))



inbreeding <- data_frame(f = diag(m), ID = ids$ID)

ggplot(inbreeding, aes(x = ID, y = f)) + geom_point()+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5)) + xlab("Individual") + ylab("Inbreeding coefficient (F)")

ggplot(inbreeding %>% filter(ID != "ZOO7"), aes(x = ID, y = f)) + geom_point()+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5)) + xlab("Individual") + ylab("Inbreeding coefficient (F)")




ids <- read_tsv(file = "../outputs/plink_stats/fb-170208-snp-q30-gq30-hwe-ss0.001.rel.id", col_names=c("ID", "ID2"))

dat <- read_tsv(file = "../outputs/plink_stats/fb-170208-snp-q30-gq30-hwe-ss0.001.rel", col_names=ids$ID)

col3 <- colorRampPalette(c("red", "white", "blue"))

m <- as.matrix(dat)
dimnames(m)[[1]] <- dimnames(m)[[2]]
m[upper.tri(m)] <- t(m)[upper.tri(m)]

m2 <- m
diag(m2) <- 0
corrplot(m2[-34,-34], is.corr=FALSE, diag = FALSE, method = "color", order = 'hclust', addrect = 4, rect.col = "black", tl.col="black")
corrplot(m, is.corr=FALSE, diag = TRUE, method = "square", order = 'hclust')

rel <- read_tsv(file = "../outputs/plink_stats/fb-170208-snp-q30-gq30-hwe-ss0.001-ldp.grm", col_names = c("ID1", "ID2", "nsnps", "rel"))
inds <- read_tsv(file = "../outputs/plink_stats/fb-170208-snp-q30-gq30-hwe-ss0.001-ldp.grm.id", col_names = c("Ind", "ID"))
inds$ID <- 1:nrow(inds)


rel <- left_join(rel, inds %>% select(Ind1 = Ind, ID1 = ID))
rel <- left_join(rel, inds %>% select(Ind2 = Ind, ID2 = ID))



inbreeding <- data_frame(f = diag(m), ID = ids$ID)

ggplot(inbreeding %>% filter(ID != "ZOO7"), aes(x = ID, y = f)) + geom_point()+ theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5)) + xlab("Individual") + ylab("Inbreeding coefficient (F)")









ids <- read_delim(file = "../outputs/ldak/fb-170208-snp-q30-gq30-hwe-ss0.001.grm.id", col_names=c("ID", "ID2"), delim = " ")
dat <- read_delim(file = "../outputs/ldak/fb-170208-snp-q30-gq30-hwe-ss0.001.grm.raw", delim = " ", col_names=FALSE)
m <- as.matrix(dat[-ncol(dat)])
dimnames(m)[[1]] <- ids$ID
dimnames(m)[[2]] <- ids$ID

m2 <- m[-which(dimnames(m)[[1]] == "ZOO7"), -which(dimnames(m)[[1]] == "ZOO7")]

corrplot(m2, is.corr=FALSE, diag = FALSE, method = "color", order = 'hclust', addrect = 4, rect.col = "black", tl.col="black")







col3 <- colorRampPalette(c("red", "white", "blue"))



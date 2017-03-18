
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)



pc_cols <- paste("PC", 1:20, sep = '')
dat <- read_delim(file="../outputs/ldak/fb-170208-snp-q30-gq30-hwe-ss0.001.vect", delim=" ", col_names=c("ID", "ID2", pc_cols))
dat <- dat %>% select(-ID2)

sgs <- read_tsv("../sample_groups.tsv", col_names=c("ID", "Pop"))

dat <- left_join(dat, sgs)

ggplot(dat, aes(x = PC1, y = PC2, color = Pop)) + geom_point(size = 3)

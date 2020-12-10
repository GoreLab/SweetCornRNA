library(tidyverse)


carotDat <- read.table("~/Dropbox/Sweet corn - papers and data/Data/cleaned_HPLC_carot_14.15_20170609.txt", header = T, sep = "\t", na.strings = "NA")
tocoDat <- read.table("~/Dropbox/Sweet corn - papers and data/Data/cleaned_HPLC_toco_14.15_20170406.txt", header = T, sep = "\t", na.strings = "NA")

length(unique(carotDat$Geno))
length(unique(tocoDat$Geno))

toco_final <- read.csv("~/Non_icloud_files/toco_analysis_2017/toco_back_transformed_blups.csv")
toco_final$Sample.ID

carot_final <- read.csv("~/Documents/GitHub/SweetCornRNA/data/carotenoid_back_transformed_BLUPs-converted.csv",
                        skip = 1)
carot_final <- carot_final[1:(nrow(carot_final)-1),]
carot_final$Sample.ID
# toco started with 419 lines
# 384 with toco data in final paper

# carotenoids started with 422 lines
# carotenoids ended with 308 lines in paper

lines_with_carot_and_or_toco <- unique(c(as.character(carot_final$Sample.ID), as.character(toco_final$Sample.ID)))

# genos
genos <- read.table("~/Documents/GitHub/SweetCornRNA/data/Ion_163K_401_v4.hmp.txt", nrows = 1 )
genos <- genos[1,12:ncol(genos)]
genos <- unname(unlist(genos))
genos.split <- strsplit(as.character(genos), "[.]")
genos.simple <- unlist(map(genos.split, 1))

genos.simple <- gsub('^X','',genos.simple)


genos.simple
lines_with_carot_and_or_toco
genos_no_phenos <- setdiff(genos.simple, lines_with_carot_and_or_toco) # in x but not y
phenos_no_genos <- setdiff(lines_with_carot_and_or_toco, genos.simple)

write.csv(as.data.frame(genos.simple ), "~/Desktop/geno_pheno_matches.csv", row.names = F)


genos2 <- read.table("~/Documents/GitHub/SweetCornRNA/data/GBSvmatch_v2_v4_pos_20180606.txt", nrows = 3 )


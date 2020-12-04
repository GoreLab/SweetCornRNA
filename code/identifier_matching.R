# Match RNA rlog counts, metadata, and toco/carotenoid BLUP identifiers for analysis
# Jenna Hershberger
# jmh579@cornell.edu
# 12/02/2020

library(tidyverse)
library(readxl)

rna_metadata <- read_excel("./data/Sweet_Corn_TagSeq_Sample_Tracking_All_Sequencing.xlsx", sheet = "Sample_Tracking")
field_coordinates <- read.csv("./data/rnaseq_trial_2019_upload.csv")
rlog_ids <- read.table("~/Google Drive/Work/Sweetcorn_RNA/Sweetcorn_TagSeq/htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.txt", nrows = 1)
hplc_blup_ids <- read.csv("~/Google Drive/Work/Sweetcorn_RNA/Tocochromanols/toco_transformed_blups.csv")


hplc_blup_ids <- hplc_blup_ids %>% rename(accession_name = Sample.ID)
rlog_ids <- unlist(as.character(t(rlog_ids)))
rlog_ids <- rlog_ids[6:length(rlog_ids)]
rlog.df <- as.data.frame(rlog_ids) %>%
  separate(rlog_ids, into = c("SC", "RNA", "plate_number", "well", "plot_name", "accession", "extra"), remove = F) %>%
  mutate(accession_name = ifelse(is.na(extra), accession, paste(accession, extra, sep = "_"))) %>%
  mutate(plate_name = paste(SC, RNA, plate_number, sep = "_")) %>%
  mutate(plate_number = as.numeric(plate_number)) %>%
  dplyr::select(rlog_ids, plate_name, plate_number, well, plot_name, accession_name) %>%
  full_join(field_coordinates, by = c("plot_name", "accession_name"))
# fix accession_name s: "p39(new)
rlog.df$accession_name


# make sure accession names match hplc blups
# full_join all pheno datasets with rlog.df and field layout info

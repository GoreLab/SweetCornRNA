# Match RNA rlog counts, metadata, and toco/carotenoid BLUP identifiers for analysis
# Jenna Hershberger
# jmh579@cornell.edu
# 12/02/2020

library(tidyverse)
library(readxl)

rna_metadata <- read_excel("./data/Sweet_Corn_TagSeq_Sample_Tracking_All_Sequencing.xlsx", sheet = "Sample_Tracking")
field_coordinates <- read.csv("./data/rnaseq_trial_2019_upload.csv")
rlog_ids <- read.table("~/Google Drive/Work/Sweetcorn_RNA/Sweetcorn_TagSeq/htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.txt", nrows = 1)
toco_hplc_blup_ids <- read.csv("~/Google Drive/Work/Sweetcorn_RNA/Tocochromanols/toco_transformed_blups.csv")
carot_hplc_blup_ids <- read.csv("~/Google Drive/Work/Sweetcorn_RNA/Carotenoids/carotenoid_transformed_BLUPs-converted.csv", skip = 1)
ion_blup_ids <- read_excel("~/Google Drive/Work/Sweetcorn_RNA/Ionomics/SuppTableS3.BLUPs_Oct20.xlsx", sheet = "BLUPs")

rna_metadata$plot_name <- rna_metadata$sample_name

rlog_ids <- unlist(as.character(t(rlog_ids)))
rlog_ids <- rlog_ids[6:length(rlog_ids)]

toco_hplc_blup_ids <- toco_hplc_blup_ids %>% rename(accession_name = Sample.ID)
carot_hplc_blup_ids <- carot_hplc_blup_ids %>% rename(accession_name = Sample.ID)
ion_blup_ids <- ion_blup_ids %>% rename(accession_name = `Sample ID`)

# fix accession_name s in field file:
field_coordinates <- field_coordinates %>%
  mutate(accession_name_field = accession_name) %>%
  mutate(accession_name = str_replace(accession_name, "-", "_")) %>%
  mutate(accession_name = recode(accession_name,
                                 "LUTHER HILL" = "LUTHER_HILL",
                                 "P39(New)" = "P39_New"))

# fix accession_name s in toco hplc blups:
toco_hplc_blup_ids <- toco_hplc_blup_ids %>%
  mutate(accession_name_toco = accession_name) %>%
  mutate(accession_name = str_replace(accession_name, "Strain_", "")) %>%
  mutate(accession_name = recode(accession_name, # old = new
                                 "C23_su" = "C23",
                                 "C8_pseudo" = "C8p",
                                 "Florida_56" = "Fa56a",
                                 "IaEv191" = "IaEV191",
                                 "IL101T" = "Il101t",
                                 "IL677a" = "Il677a",
                                 "Il731a" = "IL731a",
                                 "Il767b" = "IL767b",
                                 "Il779a" = "IL779a",
                                 "IL802b_Ht2" = "IL802b",
                                 "Luther_Hill" = "LUTHER_HILL",
                                 "Me244_Wb" = "Me244_wb",
                                 "Olcott" = "OLCOTT",
                                 "P39A" = "P39a",
                                 "P39_Goodman_Buckler" = "P39Goodman_Buckler",
                                 "P39Le_253_60" = "P39Le_253",
                                 "T62S" = "T62s",
                                 "304A_408_68B" = "304A",
                                 "675A_415_68B" = "675A",
                                 "T20_2_68B" = "T20",
                                 "T24_395_68B" = "T24",
                                 "T32_397_68r" = "T32",
                                 "T33_399_68B" = "T33",
                                 "T35_388_68A" = "T35"))

# use toco changes to update carotenoid names
carot_hplc_blup_ids <- carot_hplc_blup_ids %>%
  mutate(accession_name_carot = accession_name) %>%
  mutate(accession_name = str_replace(accession_name, "Strain_", "")) %>%
  mutate(accession_name = recode(accession_name, # old = new
                                 "C23_su" = "C23",
                                 "C8_pseudo" = "C8p",
                                 "Florida_56" = "Fa56a",
                                 "IL101T" = "Il101t",
                                 "IL677a" = "Il677a",
                                 "IL802b_Ht2" = "IL802b",
                                 "Me244_Wb" = "Me244_wb",
                                 "P39_Goodman_Buckler" = "P39Goodman_Buckler",
                                 "304A_408_68B" = "304A",
                                 "675A_415_68B" = "675A",
                                 "T20_2_68B" = "T20",
                                 "T24_395_68B" = "T24",
                                 "T32_397_68r" = "T32",
                                 "T33_399_68B" = "T33",
                                 "T35_388_68A" = "T35"))

ion_blup_ids <- ion_blup_ids %>%
  mutate(accession_name_ion = accession_name) %>%
  mutate(accession_name = str_replace(accession_name, "Strain_", "")) %>%
  mutate(accession_name = recode(accession_name, # old = new
                                 "IL101T" = "Il101t",
                                 "IL14H" = "Il14H",
                                 "IL677a" = "Il677a",
                                 "IL802b_Ht2" = "IL802b",
                                 "Luther_Hill" = "LUTHER_HILL",
                                 "Me244_Wb" = "Me244_wb",
                                 "Olcott" = "OLCOTT",
                                 "P39_Goodman_Buckler" = "P39Goodman_Buckler",
                                 "P39Le_253_60" = "P39Le_253",
                                 "T62S" = "T62s",
                                 "304A_408_68B" = "304A",
                                 "675A_415_68B" = "675A",
                                 "T20_2_68B" = "T20",
                                 "T24_395_68B" = "T24",
                                 "T32_397_68r" = "T32",
                                 "T33_399_68B" = "T33",
                                 "T35_388_68A" = "T35"))

# full_join all pheno datasets with rlog ids and field layout info
master.key <- as.data.frame(rlog_ids) %>%
  separate(rlog_ids, into = c("SC", "RNA", "plate_number", "well", "plot_name", "accession", "extra"), remove = F) %>%
  mutate(accession_name = ifelse(is.na(extra), accession, paste(accession, extra, sep = "_"))) %>%
  mutate(plate_name = paste(SC, RNA, plate_number, sep = "_")) %>%
  mutate(plate_number = as.numeric(plate_number)) %>%
  mutate(sample_name_metadata = plot_name) %>%
  dplyr::select(rlog_ids, sample_name_metadata, plate_name, plate_number, well, plot_name, accession_name) %>%
  full_join(field_coordinates, by = c("plot_name", "accession_name")) %>%
  full_join(toco_hplc_blup_ids[,c("accession_name", "accession_name_toco")], by = "accession_name") %>%
  full_join(carot_hplc_blup_ids[,c("accession_name", "accession_name_carot")], by = "accession_name") %>%
  full_join(ion_blup_ids[,c("accession_name", "accession_name_ion")], by = "accession_name") %>%
  full_join(rna_metadata, by = "plot_name") %>%
  filter(!is.na(well)) %>%
  arrange(plate_name, well) %>%
  dplyr::select(rlog_ids, sample_name_metadata, plot_name, accession_name, accession_name_toco,
                accession_name_carot, accession_name_ion, plate_name, plate_number, well,
                tissue_harvest_date, nucleic_acid_extraction_date)

write.csv(master.key, "./output/master_key.csv", row.names = F)


# make sure accession names match hplc blups
# rna_not_hplc <- master.key %>%
#   full_join(ion_blup_ids, by = "accession_name") %>%
#   filter(is.na(`Endosperm mutation`)) %>% dplyr::select(accession_name) %>% distinct() %>%
#   arrange(accession_name)
#
# hplc_not_rna <- master.key %>%
#   full_join(ion_blup_ids, by = "accession_name") %>%
#   filter(is.na(well)) %>% dplyr::select(accession_name) %>% distinct() %>% arrange(accession_name)

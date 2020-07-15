#!/usr/bin/env Rscript

# RNA extraction QC
# Jenna Hershberger
# jmh579@cornell.edu
# 07/14/2020

# required input:
# args[1] = working directory ex// "/Users/Jenna/Documents/GitHub/SweetCornRNA"


# example:
# Rscript rna_qc_app.R "/Users/Jenna/Documents/GitHub/SweetCornRNA"
library(tidyverse)
library(googledrive)
library(readxl)
library(lubridate)

# for debugging
#args = "/Users/Jenna/Documents/GitHub/SweetCornRNA"

args = commandArgs(trailingOnly=TRUE)

# test if there are three arguments: if not, return an error
if (length(args)<1) {
  stop("An input argument is required: Working directory", call.=FALSE)
}

#working.directory <- args[1]
working.directory <- getwd()
setwd(working.directory) #"/Users/Jenna/Documents/GitHub/SweetCornRNA"


#### GRINDING LOG ####
drive_download(file = "SC_grinding_2020", overwrite = T) # saves to working directory
grinding_log <- read_excel("SC_grinding_2020.xlsx", na = "NA") %>%
  separate(Polycon_barcode, into = c("Sample_barcode", "rep"), sep = "_", remove = FALSE) %>%
  mutate(Sample_barcode = case_when(Sample_barcode == "19A0012" & Notes == "from sample harvested 8/22" ~ Polycon_barcode,
                                    Sample_barcode == "19A0012" & Notes != "from sample harvested 8/22" ~ "19A0012_2",
                                    Sample_barcode == "19A0096" & Notes == "from sample harvested 9/6/2019" ~ Polycon_barcode,
                                    Sample_barcode == "19A0096" & Notes == "from sample harvested 9/11/2019" ~ "19A0096_2",
                                    TRUE ~ Sample_barcode)) %>%
  rename(grind_notes = Notes) %>%
  filter(!is.na(`5ml_Rack_Number`)) %>%
  dplyr::select(Sample_barcode, Grind_Date, `5ml_Rack_Number`, grind_notes)

#### EXTRACTIONS LOG ####
drive_download(file = "rna_extractions_log.xlsx", overwrite = T) # saves to working directory
raw_firstpass <- read_excel("rna_extractions_log.xlsx", sheet = "FirstPass", col_names = TRUE,
                            range = cell_cols("A:H"), na = "NA") %>%
  rename(extraction_notes = Notes)
raw_QC <- read_excel("rna_extractions_log.xlsx", sheet = "QC", col_names = TRUE,
                     range = cell_cols("A:P"), na = "NA") %>% rename(Sample_barcode = Name)
raw_gel <- read_excel("rna_extractions_log.xlsx", sheet = "Gels", col_names = TRUE,
                     range = cell_cols("A:G"), na = "NA") %>%
  filter(Gel_num != 1 & Gel_num != 3)


#### MATCH AND SUMMARIZE ####
matched <- raw_firstpass %>%
  full_join(raw_QC, by = "Sample_barcode") %>%
  full_join(raw_gel, by = "Sample_barcode") %>%
  full_join(grinding_log, by = "Sample_barcode")

summary_df <- matched %>%
  dplyr::select(Sample_barcode, Extraction_day_1_date, done_box_number, extraction_notes,
                QC_date, `260/280`:Gel_date, gel_pass, gel_notes, Grind_Date, `5ml_Rack_Number`,
                grind_notes)
redos <- summary_df %>%
  filter(gel_pass == FALSE | QC_pass == FALSE) %>%
  mutate(qc_matches = ifelse(gel_pass == QC_pass, TRUE, FALSE))

#write.csv(redos, "redos.csv", row.names = F)
write.csv(redos, paste0(working.directory, "/rna_qc_redos", today(), ".csv"), row.names = F)
write.csv(summary_df, paste0(working.directory, "/rna_qc_summary_", today(), ".csv"), row.names = F)



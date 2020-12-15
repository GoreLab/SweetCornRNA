# Formatting metadata for Buell Lab

library(tidyverse)
library(readxl)

plate_layout <- read.csv("~/Documents/Github/SweetCornRNA/output/SC_RNA_plate_master_for_submission.csv", skip = 1) %>%
  rename(sample_name = Sample.Name)
field_layout <- read.csv("~/Documents/GitHub/GoreLabBase/2019/formatted_for_upload/trials/rnaseq_trial_2019_upload.csv") %>%
  rename(sample_name = plot_name)

# 19A0012 == 19A0012_2
# 19A0096 == 19A0096_2
# 19A0059 == 19A0059_1
qc_summary <- read.csv("~/Documents/Github/SweetCornRNA/output/rna_qc_summary_2020-07-23.csv") %>%
  rename(sample_name = Sample_barcode) %>%
  mutate(sample_name = recode(sample_name,
                              "19A0012_2" = "19A0012",
                              "19A0096_2" = "19A0096",
                              "19A0059_1" = "19A0059"
  ))
harvest_summary <- read_excel("~/Documents/GitHub/SweetCornRNA/data/crosses_2019.xlsx",
                              col_types = c("date", "text", "numeric", "date", "text")) %>%
  rename(sample_name = plot_id) %>% drop_na(harvested)

layout_temp <- full_join(plate_layout, field_layout, by = "sample_name") %>%
  janitor::clean_names()

layout_temp.2 <- full_join(layout_temp, qc_summary, by = "sample_name")

layout_full <- full_join(layout_temp.2, harvest_summary, by = "sample_name") %>%
  janitor::clean_names()

#write.csv(layout_full, "~/Documents/Github/SweetCornRNA/output/SC_metadata.csv", row.names = F)

meta_template <- read_excel("~/Documents/Github/SweetCornRNA/data/Sweet_Corn_TagSeq_Sample_Tracking_All_Sequencing.xlsx")

meta <- layout_full %>%
  # filter(sample_name == "19A0012" & notes != "Also harvested 8/1/2019",
  #        )
  rename(ncbi_taxon_id = ncbi_tax_id) %>%
  mutate(library_name = NA,
         organism = "Zea mays",
         cultivar_or_isolate_or_ecotype = "Cultivar",
         tissue = "seed",
         tissue_developmental_stage = "400 GDD after pollination",
         tissue_harvest_date = lubridate::ymd(harvested),
         plant_age = as.numeric(tissue_harvest_date - lubridate::ymd("2019-06-13")),
         treatment = "None",
         specimen_source = "Gore lab",
         growing_location = "Musgrave Farm, Aurora, NY",
         growing_conditions = "Non-irrigated",
         geographic_location = "Musgrave Farm, Aurora, NY",
         biological_replicate = 1,
         tissue_harvested_by = "Jenna Hershberger",
         sample_tissue_comments = notes,
         nucleic_acid_extraction_method = "Hot Borate",
         dnase_rnase_treatment = "yes",
         nucleic_acid_clean_method = NA,
         `nucleic_acid_concentration_via_qubit_ng/ul` = ng_m_l,
         `nanodrop_260/280` = x260_280,
         `nanodrop_260/230` = x260_230,
         nucleic_acid_volume_ul = volume_ul,
         nucleic_acid_buffer = "Water",
         nucleic_acid_extracted_by = "Jenna Hershberger",
         nucleic_acid_extraction_date = extraction_day_1_date,
         nucleic_acid_image = ifelse(is.na(gel_num), FALSE, TRUE),
         Plate = plate_name,
         Well = well,
         nucleic_acid_comments = extraction_notes,
         well_letter = substr(Well, 1,1),
         well_number = parse_number(as.character(Well)),
         sample_title_DO_NOT_EDIT = NA,
         library_method = NA,
         library_pool = NA,
         `library_concentration_via_qubit_ng/ul` = NA,
         library_volume_ul = NA,
         library_buffer = NA,
         library_made_by = NA,
         library_preparation_date = NA,
         library_image = NA,
         library_insert_size_bp = NA,
         library_strategy = NA,
         library_source = NA,
         library_selection = NA,
         library_index = NA,
         library_comments = NA,
         sequencing_submission_form = NA,
         date_submitted_for_sequencing = lubridate::ymd("2020-07-31"),
         sequencing_center = "Cornell University Institute of Biotechnology (COR)",
         sequencing_platform = "ILLUMINA",
         `instrument_model *See "Instrument_Model_Options" tab` = "NextSeq 500",
         sequence_length_nt = 86,
         layout = "single",
         total_reads_requested_million = NA,
         sequencing_date = NA,
         purity_filtered_reads = NA,
         flow_cell = NA) %>%
  dplyr::arrange(Plate, well_letter, well_number) %>%
  dplyr::select(colnames(meta_template)[1:32]) %>%
  drop_na(Well, Plate)

write.csv(meta, "~/Documents/GitHub/SweetCornRNA/output/metadata_for_Josh.csv", row.names = F, na = "")









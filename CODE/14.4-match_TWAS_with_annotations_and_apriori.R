# Compile TWAS results and match with annotations, a priori candidate lists, and NAM JL QTL intervals
# Run "match_apriori_with_annotations.R" first
# 5/11/2021

library(tidyverse)
library(magrittr)
library(fs)

source("14.3-TWAS_top_hits_apriori_FUNCTIONS.R")

#### Read in files ####
trait_key <- read.csv("RAWDATA/Metadata/trait_key.csv")
element.traits <- c("Boron", "Cadmium", "Calcium", "Copper", "Magnesium",
                    "Manganese", "Molybdenum", "Nickel", "Phosphorus", "Potassium",
                    "Rubidium", "Sulfur", "Zinc",
                    "Iron", "Strontium")

# Annotations and a priori candidates
apriori.df <- read.csv("./RESULT/14.2-Match_apriori_with_annotations/Apriori_annotated_20210518.csv", na = "")
annotations <- read.csv("./RESULT/14.2-Match_apriori_with_annotations/annotations_with_JL.csv") %>% dplyr::select(-Chr)

# TWAS results
twas_results_dir <- "./RESULT/7.1-TWAS"
twas_files <- list.files(twas_results_dir, pattern = "B73", full.names = T)
twas.df <- map_dfr(twas_files, read_csv_adapted_B73) %>%
  process_infile_B73(input.df = ., exclude.vec = element.traits, apriori.df = apriori.df,
                     annotations.df = annotations) %>%
  left_join(trait_key)

twas.JL <- twas.df %>%
  filter(str_detect(as.character(NAM.JL.QTL.traits), as.character(trait))) %>%
  dplyr::select(RefGen_v4.Gene.ID:Wu_21, pathway, QTL_ID:NAM.JL.QTL.traits,
                starts_with("Common.")) %>%
  distinct()
twas.noJL <- twas.df %>%
  dplyr::select(RefGen_v4.Gene.ID:Wu_21, pathway) %>% distinct()

twas.2 <- left_join(twas.noJL, twas.JL) %>% arrange(trait, rank)

write.csv(twas.2, "./RESULT/14.4-match_TWAS_with_annotations_and_apriori/TWAS_all_genes.csv", row.names = F)

#### Pathway FDR calculation ####
twas_carot <- twas.2 %>%
  filter(!is.na(Diepenbrock_21),
         pathway == "carot") %>%
  mutate(pathway.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr")) %>%
  dplyr::select(RefGen_v4.Gene.ID:overall.fdr, pathway.fdr, everything())

twas_toco <- twas.2 %>%
  filter(!is.na(Wu_21),
         pathway == "toco") %>%
  mutate(pathway.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr")) %>%
  dplyr::select(RefGen_v4.Gene.ID:overall.fdr, pathway.fdr, everything())

carot_ap_passing <- twas_carot %>% filter(pathway.fdr < 0.05)
toco_ap_passing <- twas_toco %>% filter(pathway.fdr < 0.05)

rbind(carot_ap_passing, toco_ap_passing) %>%
  write.csv(., "./RESULT/14.4-match_TWAS_with_annotations_and_apriori/noharv_B73_pathway_ap_FDR0.05.csv", row.names = F)

## by overall FDR ##
fdr.overall.0.05 <- twas.2 %>%
  filter(overall.fdr < 0.05) %>%
  ungroup() %>%
  distinct()
length(unique(fdr.overall.0.05$RefGen_v4.Gene.ID))
fdr.overall.0.05 %>%
  write.csv(., "./RESULT/14.4-match_TWAS_with_annotations_and_apriori/noharv_B73_overallFDR_0.05.csv", row.names = F, na = "")

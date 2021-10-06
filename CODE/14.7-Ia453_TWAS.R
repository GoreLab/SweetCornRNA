# Compile TWAS results for Ia453
# Run "match_apriori_with_annotations.R" first
# 8/20/2021

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
annotations <- read.csv("./RESULT/14.2-Match_apriori_with_annotations/Ia453_annotations_all.csv")

# TWAS results
twas_results_dir <- "./RESULT/7.1-TWAS"
twas_files <- list.files(twas_results_dir, pattern = "Ia453", full.names = T)
twas.df <- map_dfr(twas_files, read_csv_adapted_Ia453) %>%
  process_infile_Ia453(input.df = ., exclude.vec = element.traits,
                     annotations.df = annotations) %>%
  left_join(trait_key)

write.csv(twas.df, "./RESULT/14.7-Ia453_TWAS/Ia453_TWAS_all_genes.csv", row.names = F)


#### Subsets ####

## by rank ##
0.01 * 18765 # top 1% = 188 genes for B73
0.005 * 18765 # top 0.5% = 94 genes for B73

0.01 * 18477 # top 1% = 185 genes for Ia453
0.005 * 18477 # top 0.5% = 92 genes for Ia453

## by overall FDR ##
fdr.overall.0.05 <- twas.df %>%
  filter(overall.fdr < 0.05) %>%
  ungroup() %>%
  distinct() %>%
  arrange(overall.fdr) %>%
  dplyr::select(Ia453_id, neg.log.P:rank_percent, Ia453_start:Ia453_arabidopsis_hit_annotation)
length(unique(fdr.overall.0.05$Ia453_id))

# BLAST results (Ia453 transcripts BLASTed to B73v4)
blast <- read.csv("RAWDATA/Annotation/Ia453_B73v4_BLAST_key.csv", na.strings = c("", "NA")) %>%
  rename(B73v4_id = B73v4_blast)

# all B73 annotations
all_annotations_B73v4 <- read.csv("RAWDATA/Annotation/Zea_mays.B73_RefGen_v4.59_anno.csv")

blast_annotated <- blast %>%
  left_join(all_annotations_B73v4, by = c("B73v4_id" = "ID")) %>%
  drop_na(Ia453_id)

# B73 annotated
B73.apriori.df <- read.csv("./RESULT/14.2-Match_apriori_with_annotations/Apriori_annotated_20210518.csv",
                           check.names = F) %>%
  rename(B73v4_id = RefGen_v4.Gene.ID) %>%
  dplyr::select(B73v4_id, Apriori, Diepenbrock_21, Wu_21) %>%
  distinct()

fdr.overall.0.05_annotated <- fdr.overall.0.05 %>%
  left_join(blast_annotated, by = "Ia453_id") %>%
  left_join(B73.apriori.df, by = "B73v4_id") %>%
  dplyr::select(Ia453_id, B73v4_id, neg.log.P:rank_percent, Ia453_chr,
                Ia453_start:Ia453_arabidopsis_hit_annotation, e.value:percent_identity,
                gene_name:logic_name, Apriori:Wu_21) %>%
  mutate(Apriori = ifelse(!is.na(Apriori), TRUE, FALSE),
         across(.cols = everything(), ~na_if(., "")))

fdr.overall.0.05_annotated %>%
  write.csv(., "./RESULT/14.7-Ia453_TWAS/noharv_Ia453_overallFDR_0.05.csv", row.names = F, na = "")

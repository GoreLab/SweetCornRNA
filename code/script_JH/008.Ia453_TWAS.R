# 008
# Compile TWAS results for Ia453
# Run "match_apriori_with_annotations.R" first
# 8/20/2021

library(tidyverse)
library(magrittr)
library(fs)

source("code/script_JH/002.TWAS_top_hits_apriori_FUNCTIONS.R")

#### Read in files ####
trait_key <- read.csv("data/trait_key.csv")
element.traits <- c("Boron", "Cadmium", "Calcium", "Copper", "Magnesium",
                    "Manganese", "Molybdenum", "Nickel", "Phosphorus", "Potassium",
                    "Rubidium", "Sulfur", "Zinc",
                    "Iron", "Strontium")

# Annotations and a priori candidates
# apriori.df <- read.csv("./output/Ia453_annotations_apriori.csv", na = "")
annotations <- read.csv("./output/Ia453_annotations_all.csv")

# TWAS results
twas_results_dir <- "./output/NoHarvDate/7.1-TWAS"
twas_files <- list.files(twas_results_dir, pattern = "Ia453", full.names = T)
twas.df <- map_dfr(twas_files, read_csv_adapted_Ia453) %>%
  process_infile_Ia453(input.df = ., exclude.vec = element.traits,
                     annotations.df = annotations) %>%
  left_join(trait_key)

# twas.JL <- twas.df %>%
#   filter(str_detect(as.character(NAM.JL.QTL.traits), as.character(trait))) %>%
#   dplyr::select(Gene.ID:Wu_21, pathway, QTL_ID:NAM.JL.QTL.traits,
#                 starts_with("Common.")) %>%
#   distinct()
# twas.noJL <- twas.df %>%
#   dplyr::select(Gene.ID:Wu_21, pathway) %>% distinct()
#
# twas.2 <- left_join(twas.noJL, twas.JL) %>% arrange(trait, rank)

# #### Match TWAS results to annotations and NAM JL QTL intervals ####
# twas.JL <- twas.df %>%
#   filter(trait %in% trait_key$trait) %>%
#   match_JL_to_TWAS(., NAM_JL) %>%
#   distinct()
# twas.JL <- twas.JL %>% rename(overall.fdr = p.fdr) # RUN NEXT 5/11

write.csv(twas.df, "./output/NoHarvDate/Ia453_TWAS_all_genes.csv", row.names = F)
#
# apriori_not_in_twas <- twas.df %>% filter(!trait %in% trait_key$trait)
# write.csv(apriori_not_in_twas, "output/Ia453_apriori_not_in_twas.csv", row.names = F)

#
# #### Pathway FDR calculation ####
# twas_carot <- twas.2 %>%
#   filter(!is.na(Diepenbrock_21),
#          pathway == "carot") %>%
#   mutate(pathway.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr")) %>%
#   dplyr::select(Gene.ID:overall.fdr, pathway.fdr, everything())
#
# twas_toco <- twas.2 %>%
#   filter(!is.na(Wu_21),
#          pathway == "toco") %>%
#   mutate(pathway.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr")) %>%
#   dplyr::select(Gene.ID:overall.fdr, pathway.fdr, everything())
#
# carot_ap_passing <- twas_carot %>% filter(pathway.fdr < 0.05)
# toco_ap_passing <- twas_toco %>% filter(pathway.fdr < 0.05)
#
# rbind(carot_ap_passing, toco_ap_passing) %>%
#   write.csv(., "./output/NoHarvDate/noharv_Ia453_pathway_ap_FDR0.05.csv", row.names = F)


#### Subsets ####

## by rank ##
0.01 * 18765 # top 1% = 188 genes for B73
0.005 * 18765 # top 0.5% = 94 genes for B73

0.01 * 18477 # top 1% = 185 genes for Ia453
0.005 * 18477 # top 0.5% = 92 genes for Ia453
#
# twas_top_rank_0.01 <- twas.2 %>%
#   slice_max(order_by = neg.log.P, n = 185) %>% ungroup()
# twas_top_rank_0.005 <- twas.2 %>%
#   slice_max(order_by = neg.log.P, n = 92) %>% ungroup()
#
# write.csv(twas_top_rank_0.01, "./output/NoHarvDate/noharv_Ia453_top_rank_0.01.csv", row.names = F, na = "")
# write.csv(twas_top_rank_0.005, "./output/NoHarvDate/noharv_Ia453_top_rank_0.005.csv", row.names = F, na = "")
#
# twas_top_rank_0.01 %>% filter(Apriori) %>%
#   write.csv(., "./output/NoHarvDate/noharv_Ia453_top_rank_0.01_apriori.csv", row.names = F, na = "")
# twas_top_rank_0.005 %>% filter(Apriori) %>%
#   write.csv(., "./output/NoHarvDate/noharv_Ia453_top_rank_0.005_apriori.csv", row.names = F, na = "")


## by overall FDR ##
fdr.overall.0.05 <- twas.df %>%
  filter(overall.fdr < 0.05) %>%
  ungroup() %>%
  distinct() %>%
  arrange(overall.fdr) %>%
  dplyr::select(Ia453_id, neg.log.P:rank_percent, Ia453_start:Ia453_arabidopsis_hit_annotation)
length(unique(fdr.overall.0.05$Ia453_id))

# BLAST results (Ia453 transcripts BLASTed to B73v4)
blast <- read.csv("data/Ia453_B73v4_BLAST_key.csv", na.strings = c("", "NA")) %>%
  rename(B73v4_id = B73v4_blast)

# all B73 annotations
all_annotations_B73v4 <- read.csv("data/B73_v4/Zea_mays.B73_RefGen_v4.59_anno (1).csv")

blast_annotated <- blast %>%
  left_join(all_annotations_B73v4, by = c("B73v4_id" = "ID")) %>%
  drop_na(Ia453_id)

# B73 annotated
B73.apriori.df <- read.csv("./output/Apriori_annotated_20210518.csv", check.names = F) %>%
  rename(B73v4_id = RefGen_v4.Gene.ID) %>%
  dplyr::select(B73v4_id, Apriori, Diepenbrock_21, Wu_21) %>%
  distinct()

fdr.overall.0.05_annotated <- fdr.overall.0.05 %>% left_join(blast_annotated, by = "Ia453_id") %>%
  left_join(B73.apriori.df, by = "B73v4_id") %>%
  dplyr::select(Ia453_id, B73v4_id, neg.log.P:rank_percent, Ia453_chr,
                Ia453_start:Ia453_arabidopsis_hit_annotation, e.value:percent_identity,
                gene_name:logic_name, Apriori:Wu_21) %>%
  mutate(Apriori = ifelse(!is.na(Apriori), TRUE, FALSE),
         across(.cols = everything(), ~na_if(., "")))

fdr.overall.0.05_annotated %>%
  write.csv(., "./output/NoHarvDate/noharv_Ia453_overallFDR_0.05.csv", row.names = F, na = "")














twas.2 %>%
  filter(overall.fdr < 0.1) %>%
  ungroup() %>%
  distinct() %>%
  write.csv(., "./output/NoHarvDate/noharv_Ia453_overallFDR_0.1.csv", row.names = F, na = "")


length(unique(fdr.overall.0.05$Gene.ID))
length(unique(rbind(carot_ap_passing, toco_ap_passing)))

length(unique(c(as.character(fdr.overall.0.05$Gene.ID), as.character(carot_ap_passing$Gene.ID),
                as.character(toco_ap_passing$Gene.ID))))

fdr.overall.0.05 %>% filter(Apriori) #%>% dplyr::select(Gene.ID) %>% distinct()

# top1_JL <- read.csv("./output/NoHarvDate/noharv_Ia453_top1_JL.csv", na = "")
# apriori_top1 <- top1_JL %>% filter(Apriori)
# apriori_top1 %>% dplyr::select(gene_name, description, Common.Support.Interval.ID) %>% distinct()
# # Number of a priori genes per trait (not counting JL-QTL support intervals)
# apriori_top1 %>% dplyr::select(trait, "Diepenbrock_21", "Wu_21") %>%
#   mutate(apriori2 = case_when(
#     !is.na(Diepenbrock_21) ~ TRUE,
#     !is.na(Wu_21) ~ TRUE,
#     TRUE ~ FALSE
#   )) %>%
#   filter(apriori2) %>%
#   group_by(trait) %>%
#   summarize(n()) %>%
#   print(n="inf")
#
# # total number a priori genes
# apriori.df %>% filter(!is.na(Diepenbrock_21) | !is.na(Wu_21))%>%
#   dplyr::select(Gene.ID)  %>% unique()

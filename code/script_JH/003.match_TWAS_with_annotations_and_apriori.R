# Compile TWAS results and match with annotations, a apriori candidate lists, and NAM JL QTL intervals
# Run "match_apriori_with_annotations.R" first
# 5/11/2021

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
apriori.df <- read.csv("./output/Apriori_annotated_20210518.csv", na = "")
annotations <- read.csv("./output/annotations_with_JL.csv") %>% dplyr::select(-Chr)

# # NAM JL QTL intervals
# NAM_JL <- read.csv("./output/NAM_JL_toco_carot_joined.csv", check.names = F)

# TWAS results
twas_results_dir <- "./output/NoHarvDate/7.1-TWAS"
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

# #### Match TWAS results to annotations and NAM JL QTL intervals ####
# twas.JL <- twas.df %>%
#   filter(trait %in% trait_key$trait) %>%
#   match_JL_to_TWAS(., NAM_JL) %>%
#   distinct()
# twas.JL <- twas.JL %>% rename(overall.fdr = p.fdr) # RUN NEXT 5/11

write.csv(twas.2, "./output/NoHarvDate/TWAS_all_genes.csv", row.names = F)
#TODO 5.18 many rows per gene bc overlap in JL QTL. Why multiple "marker at the peak" per QTL?
# 5.19 use CSI not SI; keep peak SNP only for trait that matches in TWAS vs NAM JL
# If filter by trait in NAM QTL traits, lose ratio results (lyce)

apriori_not_in_twas <- twas.2 %>% filter(!trait %in% trait_key$trait)
write.csv(apriori_not_in_twas, "output/apriori_not_in_twas.csv", row.names = F)


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
  write.csv(., "./output/NoHarvDate/noharv_B73_pathway_ap_FDR0.05.csv", row.names = F)


#### Subsets ####

## by rank ##
0.01 * 18831 # top 1% = 188 genes for B73
0.005 * 18831 # top 0.5% = 194 genes for B73
twas_top_rank_0.01 <- twas.2 %>%
  slice_max(order_by = neg.log.P, n = 188) %>% ungroup()
twas_top_rank_0.005 <- twas.2 %>%
  slice_max(order_by = neg.log.P, n = 94) %>% ungroup()

write.csv(twas_top_rank_0.01, "./output/NoHarvDate/noharv_B73_top_rank_0.01.csv", row.names = F, na = "")
write.csv(twas_top_rank_0.005, "./output/NoHarvDate/noharv_B73_top_rank_0.005.csv", row.names = F, na = "")

twas_top_rank_0.01 %>% filter(Apriori) %>%
  write.csv(., "./output/NoHarvDate/noharv_B73_top_rank_0.01_apriori.csv", row.names = F, na = "")
twas_top_rank_0.005 %>% filter(Apriori) %>%
  write.csv(., "./output/NoHarvDate/noharv_B73_top_rank_0.005_apriori.csv", row.names = F, na = "")


## by overall FDR ##
fdr.overall.0.05 <- twas.2 %>%
  filter(overall.fdr < 0.05) %>%
  ungroup() %>%
  distinct()
length(unique(fdr.overall.0.05$RefGen_v4.Gene.ID))
fdr.overall.0.05 %>%
  write.csv(., "./output/NoHarvDate/noharv_B73_overallFDR_0.05.csv", row.names = F, na = "")

twas.2 %>%
  filter(overall.fdr < 0.1) %>%
  ungroup() %>%
  distinct() %>%
  write.csv(., "./output/NoHarvDate/noharv_B73_overallFDR_0.1.csv", row.names = F, na = "")


length(unique(fdr.overall.0.05$RefGen_v4.Gene.ID))
length(unique(rbind(carot_ap_passing, toco_ap_passing)))

length(unique(c(as.character(fdr.overall.0.05$RefGen_v4.Gene.ID), as.character(carot_ap_passing$RefGen_v4.Gene.ID),
                as.character(toco_ap_passing$RefGen_v4.Gene.ID))))

fdr.overall.0.05 %>% filter(Apriori) #%>% dplyr::select(RefGen_v4.Gene.ID) %>% distinct()

# top1_JL <- read.csv("./output/NoHarvDate/noharv_B73_top1_JL.csv", na = "")
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
#   dplyr::select(RefGen_v4.Gene.ID)  %>% unique()

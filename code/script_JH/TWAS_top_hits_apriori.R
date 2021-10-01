# Compile TWAS results and match with annotations, a apriori candidate lists, and NAM JL QTL intervals
# Run TWAS_NAM_JL_QTL.R first
# 4/2/2021

library(tidyverse)
library(magrittr)
library(fs)
trait_key <- read.csv("data/trait_key.csv")


read_csv_adapted_B73 <- function(path){
  file.df <- read.csv(file = path)
  pat <- "(?<=B73_)(.*)(?=.csv)"
  trait.name <- str_extract(path, pat)
  file.df$trait <- trait.name
  return(file.df)
}

# read_csv_adapted_PH207 <- function(path){
#   file.df <- read.csv(file = path)
#   pat <- "(?<=PH207_)(.*)(?=.csv)"
#   trait.name <- str_extract(path, pat)
#   file.df$trait <- trait.name
#   return(file.df)
# }
#test.df <- read_csv_adapted_PH207("./output/NoHarvDate/7.1-TWAS/TwasResult_carot_PH207_Lutein.csv")

process_infile_B73 <- function(input.df, exclude.vec, apriori.df, annotations.df){
  output.df <- input.df %>%
    rename(RefGen_v4.Gene.ID = Gene) %>%
    filter(!trait %in% exclude.vec) %>%
    arrange(trait, -neg.log.P) %>%
    group_by(trait) %>%
    mutate(p.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr"),
           rank = rank(-neg.log.P),
           rank_percent = rank/18831*100) %>%
    full_join(apriori.df) %>%
    filter(!trait %in% exclude.vec) %>%
    dplyr::select(RefGen_v4.Gene.ID:rank_percent, Apriori, Diepenbrock_21, #Diepenbrock_17,
                  Wu_21) %>%
    left_join(annotations.df) %>%
    dplyr::select(RefGen_v4.Gene.ID, neg.log.P:rank_percent, chr:logic_name, Apriori, Diepenbrock_21, #Diepenbrock_17,
                  Wu_21)
}

# process_infile_PH207 <- function(input.df, exclude.vec, apriori.df, annotations.df){
#   output.df <- input.df %>%
#     dplyr::select(Gene, neg.log.P, trait) %>%
#     filter(!trait %in% exclude.vec) %>%
#     arrange(trait, -neg.log.P) %>%
#     group_by(trait) %>%
#     distinct() %>%
#     mutate(p.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr"),
#            rank = rank(-neg.log.P),
#            rank_percent = rank/18831*100) %>%
#     left_join(input.df, by = c("Gene", "trait", "neg.log.P")) %>%
#     filter(!trait %in% exclude.vec) %>%
#     rename(PH207_gene = Gene) %>%
#     full_join(apriori.df, by = "PH207_gene") %>%
#     dplyr::select(PH207_gene, RefGen_v4.Gene.ID, neg.log.P, p.fdr, trait, rank, rank_percent, Apriori, Diepenbrock_21, Diepenbrock_17, Wu_21) %>%
#     full_join(annotations.df, by = c("PH207_gene", "RefGen_v4.Gene.ID")) %>%
#     # rename(B73_chr = chr,
#     #        B73_start = start,
#     #        B73_end = end,
#     #        B73_strand = strand,
#     #        B73_type = type) %>%
#     dplyr::select(PH207_gene, RefGen_v4.Gene.ID, neg.log.P:rank_percent,
#                   #B73_chr:logic_name, starts_with("PH207"), Apriori, Diepenbrock_21, Diepenbrock_17, Wu_21)
#                   chr:logic_name, starts_with("PH207"), Apriori, Diepenbrock_21, Diepenbrock_17, Wu_21) %>%
#     drop_na(PH207_gene) %>%
#     distinct()
# }

element.traits <- c("Boron", "Cadmium", "Calcium", "Copper", "Magnesium",
                    "Manganese", "Molybdenum", "Nickel", "Phosphorus", "Potassium",
                    "Rubidium", "Sulfur", "Zinc",
                    "Iron", "Strontium")


#### Annotations and a priori candidates ####

# PH207_coordinates <- read.csv("./data/PH207/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0.v03eksc_anno.csv") %>%
#   rename(PH207_gene = ID, PH207_chr = chr, PH207_start = start, PH207_end = end, PH207_strand = strand, PH207_type = type)
# PH207_B73_key <- read_csv("./data/PH207/B73_vs_PH207_collinear_genes.csv", na = "", col_types = "ccccc") %>%
#   pivot_longer(cols = starts_with("B73"), names_to = "label", values_to = "RefGen_v4.Gene.ID") %>%
#   dplyr::select(-label) %>%
#   #add_row(PH207_gene = "Zm00008a037992", RefGen_v4.Gene.ID = "Zm00001d003513") %>% drop_na() %>%
#   full_join(PH207_coordinates)

v4_annotations <- read.csv("./data/B73_v4/Zea_mays.B73_RefGen_v4.59_anno (1).csv") %>%
  rename(RefGen_v4.Gene.ID = ID)

# full_annotations <- full_join(v4_annotations, PH207_B73_key)

apriori.df <- read.csv("./output/apriori_geneIDs_v4_20210510.csv") %>%
  distinct() %>%
  left_join(v4_annotations) %>%
  # left_join(full_annotations, by = "RefGen_v4.Gene.ID") %>%
  drop_na(RefGen_v4.Gene.ID) %>%
  mutate(Apriori = TRUE,
         Annotation = as.character(Annotation),
         RefGen_v4.Gene.ID = as.character(RefGen_v4.Gene.ID)) %>%
  dplyr::select(RefGen_v4.Gene.ID, chr:logic_name, starts_with("PH207"),
                Apriori, Paper, Annotation, -Trait_group) %>%
  pivot_wider(id_cols = c(RefGen_v4.Gene.ID:Apriori),
              names_from = Paper,
              values_from = Annotation)

write.csv(apriori.df, "./output/Apriori_list_20210510.csv", row.names = F, na = "")
#apriori.df <- read.csv("./output/Apriori_list.csv", na = "")

# top 1% = 188 genes for B73
0.01 * 18831
# top 1% = 180 genes for PH207
0.01 * 17996

# #### Harvest date in model ####
# harv_twas_results_dir <- "./output/HarvDate/7.1-TWAS"
# # B73
# harv_twas_files_b73 <- list.files(harv_twas_results_dir, pattern = "B73", full.names = T)
# harv_twas_b73.df <- map_dfr(harv_twas_files_b73, read_csv_adapted_B73) %>%
#   process_infile_B73(input.df = ., exclude.vec = element.traits, apriori.df = apriori.df,
#                  annotations.df = v4_annotations)
# harv_twas_b73_sliced.df <- harv_twas_b73.df %>%
#   slice_max(order_by = neg.log.P, n = 188) %>% ungroup()
#
# # PH207
# harv_twas_files_ph207 <- list.files(harv_twas_results_dir, pattern = "PH207", full.names = T)
# harv_twas_ph207.df <- map_dfr(harv_twas_files_ph207, read_csv_adapted_PH207) %>%
#   process_infile_PH207(input.df = ., exclude.vec = element.traits, apriori.df = apriori.df,
#                  annotations.df = full_annotations)
# harv_twas_ph207_sliced.df <- harv_twas_ph207.df %>%
#   slice_max(order_by = neg.log.P, n = 180) %>% ungroup()

#### Harvest date NOT in model ####
noharv_twas_results_dir <- "./output/NoHarvDate/7.1-TWAS"
# B73
noharv_twas_files_b73 <- list.files(noharv_twas_results_dir, pattern = "B73", full.names = T)
noharv_twas_b73.df <- map_dfr(noharv_twas_files_b73, read_csv_adapted_B73) %>%
  process_infile_B73(input.df = ., exclude.vec = element.traits, apriori.df = apriori.df,
                 annotations.df = v4_annotations)
noharv_twas_b73_sliced.df <- noharv_twas_b73.df %>%
  slice_max(order_by = neg.log.P, n = 188) %>% ungroup()
noharv_twas_b73_sliced_0.05.df <- noharv_twas_b73.df %>%
  slice_max(order_by = neg.log.P, n = 94) %>% ungroup()

# # PH207
# noharv_twas_files_ph207 <- list.files(noharv_twas_results_dir, pattern = "PH207", full.names = T)
# noharv_twas_ph207.df <- map_dfr(noharv_twas_files_ph207, read_csv_adapted_PH207) %>%
#   process_infile_PH207(input.df = ., exclude.vec = element.traits, apriori.df = apriori.df,
#                  annotations.df = full_annotations) %>% distinct()
# noharv_twas_ph207_sliced.df <- noharv_twas_ph207.df %>%
#   slice_max(order_by = neg.log.P, n = 180) %>% ungroup()


#### Match to NAM JL QTL intervals ####

NAM_JL <- read.csv("./output/NAM_JL_toco_carot_joined.csv", check.names = F) %>%
  mutate(Phenotype = as.character(Phenotype))
#SC_TWAS <- read.csv("./output/TWAS_top1percent.csv", na.strings = "")


match_JL_to_TWAS <- function(SC_TWAS, NAM_JL){
  SC_TWAS <- as_tibble(SC_TWAS)
  for(i in 1:nrow(SC_TWAS)){
    keep.NAM_JL <- NAM_JL %>%
      filter(Phenotype == unlist(SC_TWAS[i, "trait"]),
             `Chr` == unlist(SC_TWAS[i,"chr"]),
             # account for partial overlaps:  start < start and end > start or start < end and end > end
             # a complete overlap should be true for both
             (`Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)` <= unlist(SC_TWAS[i,"start"])) &
               (`Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)` >= unlist(SC_TWAS[i,"start"])) |
               (`Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)` <= unlist(SC_TWAS[i,"end"])) &
               (`Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)` >= unlist(SC_TWAS[i,"end"]))) %>%
      mutate(Overlap = ifelse((`Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)` <= unlist(SC_TWAS[i,"start"])) &
                                (`Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)` >= unlist(SC_TWAS[i,"end"])),
                              "Complete", "Partial"))

    if(nrow(keep.NAM_JL) < 1){
      keep.NAM_JL <- t(matrix(data = rep(NA, ncol(NAM_JL)+1)))
      colnames(keep.NAM_JL) <- c(colnames(NAM_JL), "Overlap")
    }
    new.i <- cbind(SC_TWAS[i,], as_tibble(keep.NAM_JL))

    if(i == 1){
      new.SC_TWAS <- new.i
    } else{
      new.SC_TWAS <- rbind(new.SC_TWAS, new.i)
    }
  }
  new.SC_TWAS %<>% dplyr::select(-Phenotype, -Chr)
  return(new.SC_TWAS)
}

# harv_B73_top1_JL <- harv_twas_b73_sliced.df %>%
#   ungroup() %>%
#   match_JL_to_TWAS(., NAM_JL) %>%
#   distinct()
# write.csv(harv_B73_top1_JL, "./output/HarvDate/harv_B73_top1_JL.csv", row.names = F, na = "")
# harv_B73_top1_JL %>% filter(Apriori) %>%
#   write.csv(., "./output/HarvDate/harv_B73_top1_apriori_JL.csv", row.names = F, na = "")
#
# harv_PH207_top1_JL <- harv_twas_ph207_sliced.df %>%
#   ungroup() %>%
#   match_JL_to_TWAS(., NAM_JL) %>%
#   distinct()
# write.csv(harv_PH207_top1_JL, "./output/HarvDate/harv_PH207_top1_JL.csv", row.names = F, na = "")
# harv_PH207_top1_JL %>% filter(Apriori) %>%
#   write.csv(., "./output/HarvDate/harv_B73_top1_apriori_JL.csv", row.names = F, na = "")

apriori_not_in_twas <- noharv_twas_b73.df %>% filter(!trait %in% trait_key$trait)
write.csv(apriori_not_in_twas, "output/apriori_not_in_twas.csv", row.names = F)

noharv_twas_b73_JL <- noharv_twas_b73.df %>%
  filter(trait %in% trait_key$trait) %>%
  match_JL_to_TWAS(., NAM_JL) %>%
  distinct()


noharv_B73_top1_JL <- noharv_twas_b73_sliced.df %>%
  ungroup() %>%
  match_JL_to_TWAS(., NAM_JL) %>%
  distinct()
write.csv(noharv_B73_top1_JL, "./output/NoHarvDate/noharv_B73_top1_JL.csv", row.names = F, na = "")
noharv_B73_top1_JL %>% filter(Apriori) %>%
  write.csv(., "./output/NoHarvDate/noharv_B73_top1_apriori_JL.csv", row.names = F, na = "")

noharv_B73_fdrpass_JL <- noharv_twas_b73.df %>%
  filter(p.fdr < 0.05) %>%
  ungroup() %>%
  match_JL_to_TWAS(., NAM_JL) %>%
  distinct()
write.csv(noharv_B73_fdrpass_JL, "./output/NoHarvDate/noharv_B73_FDR_0.05_JL.csv", row.names = F, na = "")


noharv_B73_fdrpass0.1_JL <- noharv_twas_b73.df %>%
  filter(p.fdr < 0.1) %>%
  ungroup() %>%
  match_JL_to_TWAS(., NAM_JL) %>%
  distinct()
write.csv(noharv_B73_fdrpass0.1_JL, "./output/NoHarvDate/noharv_B73_FDR_0.1_JL.csv", row.names = F, na = "")

write.csv(noharv_twas_b73.df, "./output/NoHarvDate/TWAS_all_genes.csv", row.names = F)

# noharv_PH207_top1_JL <- noharv_twas_ph207_sliced.df %>%
#   ungroup() %>%
#   match_JL_to_TWAS(., NAM_JL) %>%
#   distinct()
# write.csv(noharv_PH207_top1_JL, "./output/NoHarvDate/noharv_PH207_top1_JL.csv", row.names = F, na = "")
# noharv_PH207_top1_JL %>% filter(Apriori) %>%
#   write.csv(., "./output/NoHarvDate/noharv_PH207_top1_apriori_JL.csv", row.names = F, na = "")

noharv_B73_top1_JL <- read.csv("./output/NoHarvDate/noharv_B73_top1_JL.csv", na = "")
apriori_top1 <- noharv_B73_top1_JL %>% filter(Apriori)
apriori_top1 %>% dplyr::select(gene_name, description, Common.Support.Interval.ID) %>% distinct()
# Number of a priori genes per trait (not counting JL-QTL support intervals)
apriori_top1 %>% dplyr::select(trait, starts_with("Diepenbrock"), "Wu_21") %>%
  mutate(apriori2 = case_when(
    !is.na(Diepenbrock_17) ~ TRUE,
    !is.na(Diepenbrock_21) ~ TRUE,
    !is.na(Wu_21) ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  filter(apriori2) %>%
  group_by(trait) %>%
  summarize(n()) %>%
  print(n="inf")

# total number a priori genes
apriori.df %>% filter(!is.na(Diepenbrock_21) | !is.na(Diepenbrock_17) | !is.na(Wu_21))%>%
  dplyr::select(RefGen_v4.Gene.ID)  %>% unique()



# new.SC_TWAS %>%
#   write.csv("./output/top1_NAM_JL.csv", row.names = F, na = "")
#
# new.SC_TWAS %>% filter(!is.na(Apriori)) %>%
#   write.csv("./output/top1_apriori_NAM_JL.csv", row.names = F, na = "")
#
# test1 <- new.SC_TWAS %>%
#   #dplyr::select(-Apriori) %>%
#   mutate(Apriori = !is.na(Diepenbrock_17) | !is.na(Diepenbrock_21) | !is.na(Wu_21)) %>% tibble() %>%
#   filter(Apriori) %>%
#   write.csv("./output/toco_carot_top1_apriori_NAM_JL.csv", row.names = F, na = "")


# # twas_apriori <- full_join(twas.df, apriori.df) %>%
# #   mutate(Annotation = as.character(Annotation)) %>% distinct()
# ranked.twas <- twas.df %>% arrange(trait, -neg.log.P) %>% group_by(trait) %>%
#   mutate(rank = rank(-neg.log.P)) #%>% nest()
#
# Zn_Fe_B <- twas.df %>% filter(trait == "Zinc" | trait == "Iron" | trait == "Boron")

# top.twas.apriori <- sliced.twas %>% filter(!is.na(Apriori)) %>% dplyr::select(-Apriori)
# write.csv(top.twas.apriori, "./output/TWAS_top1percent_with_apriori.csv", row.names = F, na="")
# write.csv(sliced.twas, "./output/TWAS_top1percent.csv", row.names = F, na = "")
#
# ranked.twas %>% filter(RefGen_v4.Gene.ID == "Zm00001d044768")
#
# not_in_TWAS <- apriori.df %>% filter(!RefGen_v4.Gene.ID %in% twas.df$RefGen_v4.Gene.ID) %>%
#   drop_na(RefGen_v4.Gene.ID)
# not_in_TWAS %>%
#   write.csv(., "./output/Apriori_genes_not_in_TWAS.csv", row.names = F, na = "")
#
#
# Josh <- read.csv("./data/htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.csv")
#
# not_in_RLOG <- apriori.df %>% filter(!RefGen_v4.Gene.ID %in% Josh$gene_id) %>%
#   drop_na(RefGen_v4.Gene.ID)
# not_in_RLOG %>%
#   write.csv(., "./output/Apriori_genes_not_in_RLOG.csv", row.names = F, na = "")
#
#
# not_in_TWAS %>% filter(!RefGen_v4.Gene.ID %in% not_in_RLOG$RefGen_v4.Gene.ID) %>%
#   write.csv(., "./output/Apriori_genes_lost_in_filtering.csv", row.names = F, na = "")





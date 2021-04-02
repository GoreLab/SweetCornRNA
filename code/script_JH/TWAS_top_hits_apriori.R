library(tidyverse)
library(fs)

read_csv_adapted <- function(path){
  file.df <- read.csv(file = path)
  pat <- "(?<=B73_)(.*)(?=.csv)"
  trait.name <- str_extract(path, pat)
  file.df$trait <- trait.name
  return(file.df)
}

twas_results_dir <- "./output/7.1-TWAS"
twas_files_b73 <- list.files(twas_results_dir, pattern = "B73", full.names = T)
twas.df <- map_dfr(twas_files_b73, read_csv_adapted) %>%
  rename(RefGen_v4.Gene.ID = Gene)

v4_annotations <- read.csv("./data/Zea_mays.B73_RefGen_v4.59_anno (1).csv") %>%
  rename(RefGen_v4.Gene.ID = ID)

apriori.df <- read.csv("./output/apriori_geneIDs_v4.csv") %>%
  distinct() %>%
  left_join(v4_annotations, by = "RefGen_v4.Gene.ID") %>%
  mutate(Apriori = TRUE,
         Annotation = as.character(Annotation),
         RefGen_v4.Gene.ID = as.character(RefGen_v4.Gene.ID)) %>%
  dplyr::select(RefGen_v4.Gene.ID, chr:logic_name, Apriori, Paper, Annotation, -Trait_group) %>%
  pivot_wider(id_cols = c(RefGen_v4.Gene.ID:Apriori),
              names_from = Paper,
              values_from = Annotation)

sliced.twas <- twas.df %>%
  arrange(trait, -neg.log.P) %>%
  group_by(trait) %>%
  mutate(rank = rank(-neg.log.P),
         rank_percent = rank/18831*100) %>% #
  full_join(apriori.df) %>%
  dplyr::select(RefGen_v4.Gene.ID:rank, Apriori:Ziegler_17) %>%
  left_join(v4_annotations) %>%
  slice_max(order_by = neg.log.P, n = 188) %>%
  dplyr::select(RefGen_v4.Gene.ID, neg.log.P:rank, chr:logic_name, Apriori:Ziegler_17)

# top 1% = 188 genes
0.01 * 18831

# twas_apriori <- full_join(twas.df, apriori.df) %>%
#   mutate(Annotation = as.character(Annotation)) %>% distinct()
ranked.twas <- twas.df %>% arrange(trait, -neg.log.P) %>% group_by(trait) %>%
  mutate(rank = rank(-neg.log.P)) #%>% nest()

Zn_Fe_B <- twas.df %>% filter(trait == "Zinc" | trait == "Iron" | trait == "Boron")

top.twas.apriori <- sliced.twas %>% filter(!is.na(Apriori)) %>% dplyr::select(-Apriori)
write.csv(top.twas.apriori, "./output/TWAS_top1percent_with_apriori.csv", row.names = F, na="")
write.csv(sliced.twas, "./output/TWAS_top1percent.csv", row.names = F, na = "")

write.csv(apriori.df, "./output/Apriori_list.csv", row.names = F, na = "")


ranked.twas %>% filter(RefGen_v4.Gene.ID == "Zm00001d044768")

not_in_TWAS <- apriori.df %>% filter(!RefGen_v4.Gene.ID %in% twas.df$RefGen_v4.Gene.ID) %>%
  drop_na(RefGen_v4.Gene.ID)
not_in_TWAS %>%
  write.csv(., "./output/Apriori_genes_not_in_TWAS.csv", row.names = F, na = "")


Josh <- read.csv("./data/htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.csv")

not_in_RLOG <- apriori.df %>% filter(!RefGen_v4.Gene.ID %in% Josh$gene_id) %>%
  drop_na(RefGen_v4.Gene.ID)
not_in_RLOG %>%
  write.csv(., "./output/Apriori_genes_not_in_RLOG.csv", row.names = F, na = "")


not_in_TWAS %>% filter(!RefGen_v4.Gene.ID %in% not_in_RLOG$RefGen_v4.Gene.ID) %>%
  write.csv(., "./output/Apriori_genes_lost_in_filtering.csv", row.names = F, na = "")





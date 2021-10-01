# Compile a priori candidate list and match with B73 v4 genome annotations

library(tidyverse)
library(readxl)
source("code/script_JH/002.TWAS_top_hits_apriori_FUNCTIONS.R")

apriori_colnames <- c("RefGen_v4.Gene.ID", "Paper", "Trait_group", "Annotation")

apriori_christine21 <- read_excel("./data/apriori/koab032_supplementary_data/tpc.01048.2020-s02.xlsx",
                                  skip = 4) %>%
  rename(Chr_v2 = RefGen_v2,
         Chr_v4 = RefGen_v4,
         `v2 Start Open Reading Frame (bp)` = ...7,
         `v2 Stop Open Reading Frame (bp)` = ...8,
         `v4 Start Open Reading Frame (bp)` = ...10,
         `v4 Stop Open Reading Frame (bp)` = ...11,
         Annotation = `Gene Name`,
         RefGen_v4.Gene.ID = `RefGen_v4 Gene ID`) %>%
  drop_na("A Priori Candidate Gene Pathway") %>%
  mutate(Paper = "Diepenbrock_21",
         Trait_group = "Carotenoid",
         RefGen_v4.Gene.ID = case_when(
           `RefGen_v2 Gene ID`  == "GRMZM2G127139" ~ "Zm00001d003513",
           `RefGen_v2 Gene ID`  == "GRMZM5G848550" ~ "Zm00001d025544",
           TRUE ~ RefGen_v4.Gene.ID
         )) %>%
  group_by(RefGen_v4.Gene.ID) %>%
  summarize_at(vars(-group_cols()), ~str_c(unique(.), collapse = ", "))

apriori_di21 <- read_excel("./data/apriori/Supplemental Table Sx. a priori gene list_DW_20210503.xlsx", skip = 2) %>%
  rename(RefGen_v4.Gene.ID = `Gene ID (RefGen v4)`) %>%
  mutate(Paper = "Wu_21",
         Trait_group = "Tocochromanol",
         Annotation = paste(`Gene name (MaizeGDB)`, `Gene annotation`, sep = "; "))

apriori_all <- rbind(apriori_christine21[, apriori_colnames], apriori_di21[,apriori_colnames])


# NAM JL QTL intervals
NAM_JL <- read.csv("./output/NAM_JL_toco_carot_joined.csv", check.names = F)

v4_annotations <- read.csv("./data/B73_v4/Zea_mays.B73_RefGen_v4.59_anno (1).csv") %>%
  rename(RefGen_v4.Gene.ID = ID)

annotations_JL <- match_JL_to_TWAS(SC_TWAS = v4_annotations, NAM_JL = NAM_JL)

write.csv(annotations_JL, "./output/annotations_with_JL.csv", row.names = F)



apriori.df <- apriori_all %>%
  distinct() %>%
  left_join(annotations_JL) %>%
  # left_join(v4_annotations) %>%
  drop_na(RefGen_v4.Gene.ID) %>%
  mutate(Apriori = TRUE,
         Annotation = as.character(Annotation),
         RefGen_v4.Gene.ID = as.character(RefGen_v4.Gene.ID)) %>%
  dplyr::select(RefGen_v4.Gene.ID, chr:logic_name,
                QTL_ID:`NAM JL-QTL traits`, `Marker at the Peak`:Overlap,
                Apriori, Paper, Annotation, -Trait_group) %>%
  pivot_wider(id_cols = c(RefGen_v4.Gene.ID:Apriori),
              names_from = Paper,
              values_from = Annotation)

write.csv(apriori.df, "./output/Apriori_annotated_20210518.csv", row.names = F, na = "")











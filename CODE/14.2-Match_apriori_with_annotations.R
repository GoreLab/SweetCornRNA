# Compile a priori candidate list and match with B73 v4 and Ia453 genome annotations

library(tidyverse)
library(readxl)
source("14.3-TWAS_top_hits_apriori_FUNCTIONS.R")

#### B73 ####

apriori_colnames <- c("RefGen_v4.Gene.ID", "Paper", "Trait_group", "Annotation")

apriori_christine21 <- read_excel("./RAWDATA/Annotation/tpc.01048.2020-s02.xlsx",
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

apriori_di21 <- read_excel("./RAWDATA/Candidate_Genes/a priori gene list_DW_20210503.xlsx", skip = 2) %>%
  rename(RefGen_v4.Gene.ID = `Gene ID (RefGen v4)`) %>%
  mutate(Paper = "Wu_21",
         Trait_group = "Tocochromanol",
         Annotation = paste(`Gene name (MaizeGDB)`, `Gene annotation`, sep = "; "))

apriori_all <- rbind(apriori_christine21[, apriori_colnames], apriori_di21[,apriori_colnames])


# NAM JL QTL intervals
NAM_JL <- read.csv("./RESULT/14.1-TWAS_NAM_JL_QTL/NAM_JL_toco_carot_joined.csv", check.names = F)

v4_annotations <- read.csv("./RAWDATA/Annotation/Zea_mays.B73_RefGen_v4.59_anno.csv") %>%
  rename(RefGen_v4.Gene.ID = ID)

annotations_JL <- match_JL_to_TWAS(SC_TWAS = v4_annotations, NAM_JL = NAM_JL)

write.csv(annotations_JL, "./RESULT/14.2-Match_apriori_with_annotations/annotations_with_JL.csv", row.names = F)



apriori.df <- apriori_all %>%
  distinct() %>%
  left_join(annotations_JL) %>%
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

write.csv(apriori.df, "./RESULT/14.2-Match_apriori_with_annotations/Apriori_annotated_20210518.csv", row.names = F, na = "")


#### Ia453 ####
library(tidyverse)
library(janitor)

# Read in files: B73v4 annotated, Ia453 annotations, B73v5-Ia453 key, B73v4-v5 key
# apriori.df <- read.csv("./output/Apriori_annotated_20210518.csv", check.names = F)
B73.apriori.df <- apriori.df %>%
  rename_with(~paste0("B73_", .x), .cols = chr:strand) %>%
  rename(B73v4_id = RefGen_v4.Gene.ID)
Ia453_annotations <- read.csv("./data/Ia453_annotations.csv") %>%
  rename(Ia453_id = gene_id,
         Ia453_chr = chr,
         Ia453_start = pos_left,
         Ia453_end = pos_right,
         Ia453_arabidopsis_hit = best_arabidopsis_hit,
         Ia453_arabidopsis_hit_annotation = best_arabidopsis_hit_annotation)

# Split B73v4-Ia453 key file so each matching pair has its own line
key.B73v5.Ia453 <- read.table("RAWDATA/Annotation/gffcmp.liftoff_B73.txt") %>%
  janitor::row_to_names(row_number = 1) %>%
  separate(col = qry_id_list,
           into = c("match1", "match2", "match3", "match4", "match5", "match6", "match7"),
           sep = ",") %>%
  pivot_longer(cols = starts_with("match"),
               names_to = "match_order",
               values_to = "B73v5_long",
               values_drop_na = TRUE) %>%
  separate(col = "B73v5_long",
           into = c("B73v5_gene", "B73v5_id_temp"),
           sep = "\\|") %>%
  dplyr::select(ref_gene, B73v5_gene) %>%
  distinct() %>%
  rename(Ia453_id = ref_gene,
         B73v5_id = B73v5_gene)

# Match B73v4-v5 and B73v5-Ia453 key files
B73_v4_v5_key <- read.csv("RAWDATA/Annotation/B73v4_to_B73v5_genes_all.csv", na.strings = "") %>%
  rename(B73v5_id = v5.Gene.Model.ID,
         B73v4_id = v4.Gene.Model.ID,
         B73v3_id = v3.Gene.Model.ID)

key.B73v4.Ia453 <- full_join(key.B73v5.Ia453, B73_v4_v5_key)

# Join a priori annotations with key and then Ia453 annotation file
joined_annotations_all <- B73.apriori.df %>%
  full_join(key.B73v4.Ia453) %>%
  full_join(Ia453_annotations) %>%
  dplyr::select(starts_with("B73"), starts_with("Ia453"), everything()) %>%
  # replace all " , " with ", " and remove all commas and spaces before first letter
  mutate(`NAM Families with QTL` = as.character(`NAM Families with QTL`)) %>%
  mutate(`NAM Families with QTL` = str_replace_all(`NAM Families with QTL`, pattern = " , ", replacement =  ", ")) %>%
  mutate(`NAM Families with QTL` = str_replace(string = `NAM Families with QTL`, pattern = "^, ", replacement =  ""))

joined_annotations <- joined_annotations_all %>% filter(Apriori)

write.csv(joined_annotations_all, "RESULT/14.2-Match_apriori_with_annotations/Ia453_annotations_all.csv", row.names = F)
write.csv(joined_annotations, "RESULT/14.2-Match_apriori_with_annotations/Ia453_annotations_apriori.csv", row.names = F)

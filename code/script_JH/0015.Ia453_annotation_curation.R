# 8/20/21
library(tidyverse)
library(janitor)

# Read in files: B73v4 annotated, Ia453 annotations, B73v5-Ia453 key, B73v4-v5 key

B73.apriori.df <- read.csv("./output/Apriori_annotated_20210518.csv", check.names = F) %>%
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
key.B73v5.Ia453 <- read.table("output/gffcmp.liftoff_B73 copy.txt") %>%
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
B73_v4_v5_key <- read.csv("data/B73v4_to_B73v5_genes_all.csv", na.strings = "") %>%
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

write.csv(joined_annotations_all, "output/Ia453_annotations_all.csv", row.names = F)
write.csv(joined_annotations, "output/Ia453_annotations_apriori.csv", row.names = F)


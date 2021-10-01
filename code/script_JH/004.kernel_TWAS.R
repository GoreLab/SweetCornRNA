library(tidyverse)

# Annotations
ia453_rlog <- read.csv("./data/Ia453_annotations.csv")
annotations <- read.csv("./output/annotations_with_JL.csv") %>% dplyr::select(-Chr)

# TWAS results
ia453_kernelTWAS <- read.csv("./output/NoHarvDate/9.7-TWAS_mutant_type/TwasResult_toco_Ia453_endosperm.mutation.csv") %>%
  rename(gene_id = Gene) %>%
  mutate(overall.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr"),
         rank = rank(-neg.log.P),
         rank_percent = rank/18477*100) %>%
  left_join(ia453_rlog) %>%
  arrange(rank) %>%
  dplyr::select(gene_id, chr:pos_right, neg.log.P:rank_percent, best_arabidopsis_hit_annotation)
b73_kernelTWAS <- read.csv("./output/NoHarvDate/9.7-TWAS_mutant_type/TwasResult_toco_B73_endosperm.mutation.csv") %>%
  rename(RefGen_v4.Gene.ID = Gene) %>%
  mutate(overall.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr"),
         rank = rank(-neg.log.P),
         rank_percent = rank/18765*100) %>%
  left_join(annotations) %>%
  dplyr::select(RefGen_v4.Gene.ID:logic_name) %>%
  distinct() %>%
  arrange(rank)

write.csv(ia453_kernelTWAS, "output/NoHarvDate/kernel_mutant_Ia453_annotated.csv", na = "", row.names = F)
write.csv(b73_kernelTWAS, "output/NoHarvDate/kernel_mutant_B73_annotated.csv", na = "", row.names = F)

ia453_kernelTWAS %>% filter(overall.fdr < 0.05)
b73_kernelTWAS %>% filter(overall.fdr < 0.05)


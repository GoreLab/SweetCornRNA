# Pathway-level FDR
library(tidyverse)
trait_key <- read.csv("data/trait_key.csv")
twas_all <- read.csv("output/NoHarvDate/TWAS_all_genes.csv") %>%
  left_join(trait_key, by = "trait")

twas_ap_carot <- twas_all %>%
  filter(!is.na(Diepenbrock_21),
         pathway == "carot") %>%
  mutate(pathway.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr")) %>%
  dplyr::select(RefGen_v4.Gene.ID:p.fdr, pathway.fdr, everything())

twas_ap_toco <- twas_all %>%
  filter(!is.na(Wu_21),
         pathway == "toco") %>%
  mutate(pathway.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr")) %>%
  dplyr::select(RefGen_v4.Gene.ID:p.fdr, pathway.fdr, everything())


carot_ap_passing <- twas_ap_carot %>% filter(pathway.fdr < 0.05)
toco_ap_passing <- twas_ap_toco %>% filter(pathway.fdr < 0.05)

passing <- rbind(carot_ap_passing, toco_ap_passing)


write.csv(passing, "./output/NoHarvDate/noharv_B73_pathway_ap_FDR5.csv", row.names = F)

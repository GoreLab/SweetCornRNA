library(tidyverse)
library(readxl)
baseggio_2019_sig_SNPs <- read_excel("data/apriori/Baseggio_2019_significant_SNPs-converted.xlsx", sheet = "Table 1", skip = 1)
baseggio_chr2 <- baseggio_2019_sig_SNPs %>% filter(Chr. == 2)


v4_matrix <- read.table("data/Ion_163K_401_v4.hmp.txt", sep = "\t", header = T)
v4_matrix.cut <- v4_matrix %>%
  dplyr::select(rs:pos) %>% rename(`SNP ID` = rs)

twas_hits <- read.csv("output/NoHarvDate/noharv_B73_overallFDR_0.05.csv")
twas_chr2 <- twas_hits %>% filter(chr == 2)
snp_matches <- baseggio_chr2 %>% left_join(v4_matrix.cut) %>%
  mutate(dist_start_Zm00001d007121 = abs(222781921 - pos), # CW-type Zinc Finger
         dist_end_Zm00001d007121 = abs(222806658 - pos),
         dist_start_Zm00001d007489 = abs(232843576 - pos), # Zinc-finger domain of monoamine-oxidase A repressor R1
         dist_end_Zm00001d007489 = abs(232847278 - pos)) %>%
  mutate(kb_start_Zm00001d007121 = abs(222781921 - pos)/1000, # CW-type Zinc Finger
         kb_end_Zm00001d007121 = abs(222806658 - pos)/1000,
         kb_start_Zm00001d007489 = abs(232843576 - pos)/1000, # Zinc-finger domain of monoamine-oxidase A repressor R1
         kb_end_Zm00001d007489 = abs(232847278 - pos)/1000)

remove(v4_matrix)

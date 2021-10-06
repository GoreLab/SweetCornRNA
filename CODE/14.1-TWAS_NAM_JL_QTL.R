# Join NAM JL QTL results from Diepenbrock et al 2017 and Diepenbrock et al 2021

library(tidyverse)
library(readxl)

JL_colnames <- c("Pathway", "Phenotype", "Common Support Interval ID", "Chr", "Marker at the Peak", "Peak Marker Position, cM",
                 "Peak Marker Position (RefGen v2)",
                 "Peak Marker Position (RefGen v4)",
                 "Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)",
                 "Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)",
                 "Common Support Interval Left Bound (bp) (RefGen v4)",
                 "Common Support Interval Right Bound (bp) (RefGen v4)",
                 "Support Interval (α = 0.01), Left Bound (bp) (RefGen v2)",
                 "Support Interval (α = 0.01), Right Bound (bp) (RefGen v2)",
                 "Common Support Interval Left Bound (bp) (RefGen v2)",
                 "Common Support Interval Right Bound (bp) (RefGen v2)",
                 "P-value in JL model",
                 "PVE (%)",
                 "Peak Marker LOD Score",
                 "A Priori Gene in Support Interval",
                 "Number of NAM Families with QTL",
                 "NAM Families with QTL",
                 "Minimum Allelic Effect Estimate (AEE)",
                 "NAM Family with Minimum AEE",
                 "Maximum AEE",
                 "NAM Family with Maximum AEE",
                 "Number of marker-trait associations (RefGen v4)")

#### TOCOCHROMANOLS ####

NAM_GWAS_uplifted <- read_excel("./RAWDATA/Annotation/NAM_GWAS_uplifted.xlsx", skip = 1)
NAM_JL_uplifted.1 <- read_excel("./RAWDATA/Annotation/NAM_JL_uplifted.xlsx", skip = 1)
colnames(NAM_JL_uplifted.1)[16] <- "A Priori Gene in Support Interval" # having weird issues with this col name
NAM_JL_uplifted <- NAM_JL_uplifted.1 %>%
  mutate(Pathway = "toco",
         Phenotype = case_when(
    Phenotype == "γΤ" ~ "gT",
    Phenotype == "αΤ" ~ "aT",
    Phenotype == "δΤ" ~ "dT",
    Phenotype == "γΤ3" ~ "gT3",
    Phenotype == "αΤ3" ~ "aT3",
    Phenotype == "δΤ3" ~ "dT3",
    Phenotype == "ΣTT3" ~ "Total.T3_plus_T",
    Phenotype == "ΣT3" ~ "Total.T3",
    Phenotype == "ΣT" ~ "Total.T"
  )) %>%
  rename(`Peak Marker Position (RefGen v2)` = `Peak marker pos v2 (bp)`,
         `Peak Marker Position (RefGen v4)` = `Peak marker pos v4 (bp)`,
         `Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)` = `Left SI pos v4 (bp)`,
         `Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)` = `Right SI pos v4 (bp)`,
         `Common Support Interval Left Bound (bp) (RefGen v4)` = `Left CSI pos v4 (bp)`,
         `Common Support Interval Right Bound (bp) (RefGen v4)` = `Right CSI pos v4 (bp)`,
         `Support Interval (α = 0.01), Left Bound (bp) (RefGen v2)` = `Support Interval (α = 0.01), Left Bound (bp)`,
         `Support Interval (α = 0.01), Right Bound (bp) (RefGen v2)` = `Support Interval (α = 0.01), Right Bound (bp)`,
         `Common Support Interval Left Bound (bp) (RefGen v2)` = `Common Support Interval Left Bound (bp)`,
         `Common Support Interval Right Bound (bp) (RefGen v2)` = `Common Support Interval Right Bound (bp)`,
         `Number of marker-trait associations (RefGen v2)` = `Number of SNP-trait associations`,
         `Chr` = `Peak marker chr v2`) %>%
  mutate(`Number of marker-trait associations (RefGen v4)` = NA) %>%
  dplyr::select(all_of(JL_colnames))


#### CAROTENOIDS ####

NAM_JL_GWAS_carotenoids.1 <- read_excel("./RAWDATA/Annotation/tpc.01048.2020-s04.xlsx", skip = 13)
colnames(NAM_JL_GWAS_carotenoids.1)[7] <- "Peak Marker Position (RefGen v4)"
colnames(NAM_JL_GWAS_carotenoids.1)[6] <- "Peak Marker Position (RefGen v2)"
NAM_JL_GWAS_carotenoids <- NAM_JL_GWAS_carotenoids.1%>%
  drop_na(Traita) %>%
  filter(Traita != "phytofluene") %>%
  mutate(Phenotype = as.character(Traita),
         Pathway = "carot") %>%
  mutate(Phenotype = case_when(
    Phenotype == "a_carotene" ~ "Other.carotenes",
    Phenotype == "b_carotene" ~ "beta.Carotene",
    Phenotype == "b_cryp" ~ "beta.Cryptoxanthin",
    Phenotype == "total_carot" ~ "Total.carotenoids",
    Phenotype == "zeino" ~ "Zeinoxanthin",
    Phenotype == "zeaxanthin" ~ "Zeaxanthin",
    Phenotype == "lutein" ~ "Lutein"
  )) %>%
  drop_na(Phenotype) %>%
  rename(`Peak Marker Position, cM` = `Peak Marker Position, cMc`,
         `PVE (%)` = `PVE (%)e`,
         Chr = Chr.b,
         `Minimum Allelic Effect Estimate (AEE)` = `Minimum Allelic Effect Estimate (AEE)h`,
         `P-value in JL model` = `P-value in JL modeld`,
         `Peak Marker LOD Score` = `Peak Marker LOD Scoref`,
         `A Priori Gene in Support Interval` = `A Priori Gene in Support Intervalg`,
         `Number of marker-trait associations (RefGen v4)` = `Number of marker-trait associationsi (RefGen v4)`,
         `Number of marker-trait associations (RefGen v2)` = `Number of marker-trait associationsi (RefGen v2)`) %>%
  dplyr::select(all_of(JL_colnames))




NAM_JL_rbound <- rbind(NAM_JL_uplifted, NAM_JL_GWAS_carotenoids) %>%
  mutate(QTL_ID = paste(Pathway, `Common Support Interval ID`, sep = "_")) %>%
  dplyr::select(QTL_ID, everything())

NAM_JL_traitlist <- NAM_JL_rbound %>%
  aggregate(Phenotype ~ QTL_ID, data = ., paste, collapse = ", ")

NAM_JL <- NAM_JL_rbound %>%
  left_join(NAM_JL_traitlist, by = "QTL_ID") %>%
  rename(`NAM JL-QTL traits` = Phenotype.y) %>%
  dplyr::select(QTL_ID,
                `Common Support Interval ID`, Phenotype,
                `NAM JL-QTL traits`,
                Chr:`Peak Marker Position, cM`,
                `Peak Marker Position (RefGen v4)`:`Common Support Interval Right Bound (bp) (RefGen v4)`,
                `P-value in JL model`:`NAM Families with QTL`)

write.csv(NAM_JL, "./RESULT/14.1-TWAS_NAM_JL_QTL/NAM_JL_toco_carot_joined.csv", row.names = F)

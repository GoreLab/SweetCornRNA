library(tidyverse)
library(readxl)

SC_TWAS <- read.csv("./output/TWAS_top1percent.csv", na.strings = "")

JL_colnames <- c("Phenotype", "Common Support Interval ID", "Chr", "Marker at the Peak", "Peak Marker Position, cM",
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

NAM_GWAS_uplifted <- read_excel("./data/Supplemental Table Sx. NAM_GWAS_uplifted.xlsx", skip = 1)
NAM_JL_uplifted.1 <- read_excel("./data/Supplemental Table Sx. NAM_JL_uplifted.xlsx", skip = 1)
colnames(NAM_JL_uplifted.1)[16] <- "A Priori Gene in Support Interval" # having weird issues with this col name
NAM_JL_uplifted <- NAM_JL_uplifted.1 %>%
  mutate(Phenotype = case_when(
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

NAM_JL_GWAS_carotenoids.1 <- read_excel("./data/koab032_supplementary_data/tpc.01048.2020-s04.xlsx", skip = 13)
colnames(NAM_JL_GWAS_carotenoids.1)[7] <- "Peak Marker Position (RefGen v4)"
colnames(NAM_JL_GWAS_carotenoids.1)[6] <- "Peak Marker Position (RefGen v2)"
NAM_JL_GWAS_carotenoids <- NAM_JL_GWAS_carotenoids.1%>%
  drop_na(Traita) %>%
  filter(Traita != "phytofluene") %>%
  mutate(Phenotype = as.character(Traita)) %>%
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
         #`Peak Marker Position (RefGen v4)` = `Peak Marker Position (RefGen v4) `,
         #`Peak Marker Position (RefGen v2)` = `Peak Marker Position (RefGen v2) `,
         `PVE (%)` = `PVE (%)e`,
         Chr = Chr.b,
         `Minimum Allelic Effect Estimate (AEE)` = `Minimum Allelic Effect Estimate (AEE)h`,
         `P-value in JL model` = `P-value in JL modeld`,
         `Peak Marker LOD Score` = `Peak Marker LOD Scoref`,
         `A Priori Gene in Support Interval` = `A Priori Gene in Support Intervalg`,
         `Number of marker-trait associations (RefGen v4)` = `Number of marker-trait associationsi (RefGen v4)`,
         `Number of marker-trait associations (RefGen v2)` = `Number of marker-trait associationsi (RefGen v2)`) %>%
  dplyr::select(all_of(JL_colnames))


NAM_JL <- rbind(NAM_JL_uplifted, NAM_JL_GWAS_carotenoids)

for(i in 1:nrow(SC_TWAS)){
  keep.NAM_JL <- NAM_JL %>%
    filter(Phenotype == SC_TWAS[i, "trait"],
           `Chr` == SC_TWAS[i,"chr"],
           # account for partial overlaps:  start < start and end > start or start < end and end > end
           # a complete overlap should be true for both
           (`Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)` <= SC_TWAS[i,"start"]) & (`Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)` >= SC_TWAS[i,"start"]) |
            (`Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)` <= SC_TWAS[i,"end"]) & (`Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)` >= SC_TWAS[i,"end"])) %>%
    mutate(Overlap = ifelse((`Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)` <= SC_TWAS[i,"start"]) & (`Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)` >= SC_TWAS[i,"end"]),
           "Complete", "Partial"))
  if(nrow(keep.NAM_JL) < 1){
    keep.NAM_JL <- t(as.data.frame(rep(NA, ncol(NAM_JL_uplifted)+1)))
    colnames(keep.NAM_JL) <- c(colnames(NAM_JL_uplifted), "Overlap")
  }

  new.i <- cbind(SC_TWAS[i,], keep.NAM_JL)

  if(i == 1){
    new.SC_TWAS <- new.i
  } else{
    new.SC_TWAS <- rbind(new.SC_TWAS, new.i)
  }

}


new.SC_TWAS %>%
  write.csv("./output/top1_NAM_JL.csv", row.names = F, na = "")

new.SC_TWAS %>% filter(!is.na(Apriori)) %>%
  write.csv("./output/top1_apriori_NAM_JL.csv", row.names = F, na = "")

element.traits <- c("Boron", "Cadmium", "Calcium", "Copper", "Magnesium",
                    "Manganese", "Molybdenum", "Nickel", "Phosphorus", "Potassium",
                    "Rubidium", "Sulfur", "Zinc",
                    "Iron", "Strontium")

test1 <- new.SC_TWAS %>%
  #dplyr::select(-Apriori) %>%
  mutate(Apriori = !is.na(Diepenbrock_17) | !is.na(Diepenbrock_21) | !is.na(Wu_21)) %>% tibble() %>%
  filter(!trait %in% element.traits,
         Apriori) %>%
  dplyr::select(-Baxter_20, -Ziegler_21, -Baseggio_21) %>%
  write.csv("./output/toco_carot_top1_apriori_NAM_JL.csv", row.names = F, na = "")







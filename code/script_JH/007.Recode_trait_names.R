library(tidyverse)
library(readxl)
library(r2symbols)
S3 <- read_excel("~/Google Drive/Work/Sweetcorn_RNA/Supplemental Files for SCRNA/Supplemental Table 3.xlsx", skip = 1)

traits_recoded <- S3 %>%
  mutate(Trait = recode(Trait,
                        "beta.Carotene" = "β-carotene",
                        "beta.Cryptoxanthin" = "β-cryptoxanthin",
                        "Other.carotenes" = "Other carotenes",
                        "alpha.Xanthophylls" = "α-xanthophylls",
                        "beta.Xanthophylls" = "β-xanthophylls",
                        "Total.xanthophylls" = "Total xanthophylls",
                        "Total.carotenes" = "Total carotenes",
                        "Total.carotenoids" = "Total carotenoids",
                        "beta.Carotene_over_beta.cryptoxanthin" = "Ratio of β-carotene to β-cryptoxanthin",
                        "beta.Carotene_over_sum_beta.cryptoxanthin_and_zeaxanthin" = "Ratio of β-carotene to β-cryptoxanthin + zeaxanthin",
                        "beta.Cryptoxanthin_over_zeaxanthin" = "Ratio of β-cryptoxanthin to zeaxanthin",
                        "Zeinoxanthin_over_lutein" = "Ratio of zeinoxanthin to lutein",
                        "beta.Xanthophylls_over_alpha.xanthophylls" = "Ratio of β-xanthophylls to α-xanthophylls",
                        "Total.carotenes_over_total.xanthophylls" = "Ratio of total carotenes to total xanthophylls",
                        "aT" = "αT",
                        "aT3" = "αT3",
                        "dT" = "δT",
                        "dT3" = "δT3",
                        "gT" = "γT",
                        "gT3" = "γT3",
                        "Total.T3" = "ΣT3",
                        "Total.T" = "ΣT",
                        "Total.T3_plus_T" = "ΣT3 + ΣT",
                        "ratio_aT_gT" = "αT/γT",
                        "ratio_dT_aT" = "δT/αT",
                        "ratio_dt_gT" = "δT/γT",
                        "dT_over_sum_gT_aT" = "δT/(γT + αT)",
                        "gt_over_sum_gT_aT" = "γT/(γT + αT)",
                        "ratio_aT3_gT3" = "αT3/γT3",
                        "ratio_dT3_aT3" = "δT3/αT3",
                        "ratio_dT3_gT3" = "δT3/γT3",
                        "dT3_over_sum_gT3_aT3" = "δT3/(γT3 + αT3)",
                        "ratio_aT3_gT3" = "αT3/γT3",
                        "dT3_over_sum_gT3_aT3" = "δT3/(γT3 + αT3)",
                        "gT3_over_sum_gT3_aT3" = "γT3/(γT3 + αT3)",
                        "ratio_TotalT_TotalT3" = "ΣT/ΣT3"
                        ))
write.csv(traits_recoded, "output/traits_recoded.csv", row.names = F)

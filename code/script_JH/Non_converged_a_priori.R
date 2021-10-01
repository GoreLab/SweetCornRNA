library(tidyverse)
library(readxl)

convergence_results <- read.csv("./output/TmpFile_Convergence_Status_of_B73_genes.csv")
summary(convergence_results)

v4 <- read_csv("./data/B73_v4/gene_model_xref_v4.csv", skip = 4, skip_empty_rows = T, col_names = T) %>%
  dplyr::select(v4_gene_model, v2_gene_model) %>%
  rename(`RefGen_v2 Gene ID` = v2_gene_model,
         `RefGen_v4 Gene ID` = v4_gene_model)

apriori_colnames <- c("RefGen_v4 Gene ID", "Paper", "Trait_group", "Annotation")

nonconv_asreml <- convergence_results %>% filter(asreml != "LogLikelihood Converged")

apriori_christine21 <- read_excel("./data/apriori/koab032_supplementary_data/tpc.01048.2020-s02.xlsx",
                                  skip = 4) %>%
  rename(Chr_v2 = RefGen_v2,
         Chr_v4 = RefGen_v4,
         `v2 Start Open Reading Frame (bp)` = ...7,
         `v2 Stop Open Reading Frame (bp)` = ...8,
         `v4 Start Open Reading Frame (bp)` = ...10,
         `v4 Stop Open Reading Frame (bp)` = ...11,
         Annotation = `Gene Name`) %>%
  drop_na("A Priori Candidate Gene Pathway") %>%
  mutate(Paper = "Diepenbrock_21",
         Trait_group = "Carotenoid",
         `RefGen_v4 Gene ID` = case_when(
         `RefGen_v2 Gene ID`  == "GRMZM2G127139" ~ "Zm00001d003513",
         `RefGen_v2 Gene ID`  == "GRMZM5G848550" ~ "Zm00001d025544",
         TRUE ~ `RefGen_v4 Gene ID`
         )) %>%
  group_by(`RefGen_v4 Gene ID`) %>%
  summarize_at(vars(-group_cols()), ~str_c(unique(.), collapse = ", "))
# https://stackoverflow.com/questions/61877676/how-to-pivot-table-in-r-with-multiple-or-duplicate-text-data
#
# apriori_christine17 <- read_excel("./data/apriori/TPC2017-LSB-00475R1_Supplemental_Data_Set_1.xlsx",
#                                   skip = 3) %>%
#   mutate(`RefGen_v2 Gene ID` = str_replace(`RefGen_v2 Gene ID`, " ", "_")) %>% # Fix the one with a space
#   drop_na("A Priori Candidate Gene Pathway") %>%
#   rename(Annotation = "RefGen_v2 Annotated Gene Function") %>%
#   left_join(v4) %>%
#   mutate(Paper = "Diepenbrock_17",
#          Trait_group = "Tocochromanol")
#
# apriori_matt21 <- read_excel("./data/apriori/Baseggio2021_supplmenal_material/Tables/Table S4.xlsx", skip = 1) %>%
#   rename(`RefGen_v4 Gene ID` = `Maize gene ID`,
#          Annotation = "Rice/Arabidopsis homolog gene descriptiona") %>%
#   drop_na(`RefGen_v4 Gene ID`) %>%
#   filter(str_detect(`RefGen_v4 Gene ID`, "Zm")) %>%
#   mutate(Paper = "Baseggio_21",
#          Trait_group = "Ionomics")
#
# apriori_whitt20 <- read_excel("./data/apriori/pld3272-sup-0001-tables1.xlsx", sheet = "Z.mays") %>%
#   rename(`RefGen_v2 Gene ID` = GeneID) %>%
#   mutate(Annotation = paste0(GeneName, " (", Elements, ")")) %>%
#   left_join(v4) %>%
#   dplyr::select(`RefGen_v4 Gene ID`, Annotation) %>%
#   drop_na() %>%
#   mutate(Paper = "Whitt_20",
#          Trait_group = "Ionomics")
#
# apriori_di21 <- read_excel("./data/apriori/Supplemental Table Sx. a priori gene list_DW_20210212.xlsx", skip = 2) %>%
#   rename(`RefGen_v4 Gene ID` = `Gene ID (RefGen v4)`) %>%
#   mutate(Paper = "Wu_21",
#          Trait_group = "Tocochromanol",
#          Annotation = paste(`Gene name (MaizeGDB)`, `Gene annotation`, sep = "; "))

apriori_di21 <- read_excel("./data/apriori/Supplemental Table Sx. a priori gene list_DW_20210503.xlsx", skip = 2) %>%
  rename(`RefGen_v4 Gene ID` = `Gene ID (RefGen v4)`) %>%
  mutate(Paper = "Wu_21",
         Trait_group = "Tocochromanol",
         Annotation = paste(`Gene name (MaizeGDB)`, `Gene annotation`, sep = "; "))

# apriori_ziegler17 <- read_excel("./data/apriori/Baseggio2021_supplmenal_material/Tables/Table S7.xlsx", skip = 1) %>%
#   rename(`RefGen_v4 Gene ID` = `Gene ID`,
#          Annotation = `Annotated gene function`) %>%
#   mutate(Paper = "Ziegler_17",
#          Trait_group = "Ionomics") %>%
#   drop_na(`RefGen_v4 Gene ID`) %>%
#   distinct()

intersect(as.character(nonconv_asreml$GeneID), apriori_christine21$`RefGen_v4 Gene ID`)
# None of the filtered genes in the a priori candidates from Christine's 2021 Plant Cell paper

intersect(as.character(nonconv_asreml$GeneID), apriori_christine17$`RefGen_v4 Gene ID`)
# One is filtered out: "Zm00001d036439"

intersect(as.character(nonconv_asreml$GeneID), apriori_matt21$`RefGen_v4 Gene ID`)

intersect(as.character(nonconv_asreml$GeneID), apriori_whitt20$`RefGen_v4 Gene ID`)

intersect(as.character(nonconv_asreml$GeneID), apriori_di21$`Gene ID (RefGen v4)`)


apriori_all <- rbind(apriori_christine21[,apriori_colnames],
                     # apriori_christine17[,apriori_colnames],
                     # apriori_matt21[,apriori_colnames],
                     # apriori_whitt20[,apriori_colnames],
                     apriori_di21[,apriori_colnames]#,
                     # apriori_ziegler17[,apriori_colnames]
                     )

write.csv(apriori_all, "./output/apriori_geneIDs_v4_20210510.csv", row.names = F)
#write.csv(apriori_all, "./output/apriori_geneIDs_v4.csv", row.names = F)


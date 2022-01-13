# Compile prediction results and generate plot (Figure 3)

library(tidyverse)
library(fs)
library(ggcorrplot)

#### Compile results into one df ####

read_csv_adapted <- function(path){
  pat_trait <- "(?<=PredictionResult_)(.*)(?=.csv)"
  trait.string <- str_extract(path, pat_trait)
  trait.vec <- unlist(str_split(trait.string, pattern = "_"))
  # [1] = trait category
  # [2] = trait
  # [3] = kernel mutant type included
  # [4] = relationship matrices and candidates
  if(length(trait.vec))

  if(length(trait.vec) > 4){
    # some traits split so must rejoin here
    # beta.Xanthophylls_over_alpha.Xanthophylls
    # beta.Cryptoxanthin_over_zeaxanthin
    trait <- paste(trait.vec[2:(length(trait.vec)-2)], collapse = "_")
    } else{
      trait <- trait.vec[2]
      }

  file.df <- read.csv(file = path) %>%
    mutate(TraitCategory = trait.vec[1],
           Trait = trait,
           MutantTypeIncluded = trait.vec[(length(trait.vec) - 1)] == "UseMu",
           Model = trait.vec[length(trait.vec)])
  return(file.df)
}


# Harvest date NOT in model
# 13.1 Prediction stratified
noharv_pred_results_dir_13.1 <- "./RESULT/13.1-Prediction_stratified/"
noharv_pred_files_13.1 <- list.files(noharv_pred_results_dir_13.1, pattern = "PredictionResult", full.names = T)
noharv_pred.df_13.1 <- map_dfr(noharv_pred_files_13.1, read_csv_adapted) %>%
  filter(TraitCategory != "ion") %>% mutate(harvest.date.in.model = FALSE)
write.csv(noharv_pred.df_13.1, "./RESULT/15.1-Prediction_plots/13.1-Prediction_results_NoHarvDate_compiled.csv", row.names = F)

# 13.2 prediction stratified with candidate genes
noharv_pred_results_dir_13.2 <- "./RESULT/13.2-Prediction_stratified_useCand/"
noharv_pred_files_13.2 <- list.files(noharv_pred_results_dir_13.2, pattern = "PredictionResult", full.names = T)
noharv_pred.df_13.2 <- map_dfr(noharv_pred_files_13.2, read_csv_adapted) %>%
  filter(TraitCategory != "ion") %>% mutate(harvest.date.in.model = FALSE)
write.csv(noharv_pred.df_13.2, "./RESULT/15.1-Prediction_plots/13.2-Prediction_results_NoHarvDate_compiled.csv", row.names = F)

# all
noharv_pred.df <- rbind(noharv_pred.df_13.1, noharv_pred.df_13.2)
#### Compute accuracies from predictions and HPLC BLUPs ####

blups_carot <- read.csv("./RESULT/4.1-MakeBlueDatasets/PhenoData_carot.csv")
blups_toco <- read.csv("./RESULT/4.1-MakeBlueDatasets/PhenoData_toco.csv")

blups_carot_long <- blups_carot %>%
  pivot_longer(cols = Antheraxanthin:Total.carotenes_over_total.xanthophylls,
               names_to = "Trait", values_to = "Observed")

blups_toco_long <- blups_toco %>%
  pivot_longer(cols = aT:ratio_TotalT_TotalT3,
               names_to = "Trait", values_to = "Observed")

blups_all <- rbind(blups_toco_long, blups_carot_long)

full <- noharv_pred.df %>%
  pivot_longer(cols = starts_with("Rep"), names_to = "Fold", values_to = "Predicted") %>%
  full_join(blups_all, by = c("Sample.ID", "Trait")) %>%
  drop_na(Predicted, Observed) %>%
  mutate(Predicted = as.numeric(Predicted),
         Observed = as.numeric(Observed)) %>%
  group_by(Trait, Model, Fold, MutantTypeIncluded) %>%
  mutate(PearsonCor = cor(Predicted, Observed, method = "pearson"),
         SpearmanRankCor = cor(Predicted, Observed, method = "spearman"))

toco.summary.traits <- c("ratio_TotalT_TotalT3", "Total.T3_plus_T")
tocotrienols <- unique(full$Trait)[str_detect(unique(full$Trait), "T3")][c(1:8, 10)]
tocopherols <- unique(full$Trait)[str_detect(unique(full$Trait), "T3", negate = T)][20:28]
carotenes <- c("Other.carotenes", "beta.Carotene", "Total.carotenes")
xanthophylls <- c("Lutein", "beta.Cryptoxanthin", "Zeaxanthin", "Violaxanthin", "Zeinoxanthin",
                  "Antheraxanthin", "Total.xanthophylls", "alpha.Xanthophylls", "beta.Xanthophylls")
carot.summary.traits <- c("beta.Carotene_over_beta.cryptoxanthin",
                          "beta.Carotene_over_sum_beta.cryptoxanthin_and_zeaxanthin",
                          "beta.Cryptoxanthin_over_zeaxanthin",
                          "beta.Xanthophylls_over_alpha.xanthophylls",
                          "Total.carotenes_over_total.xanthophylls",
                          "Total.carotenoids", "Zeinoxanthin_over_lutein")

accuracy.df <- full %>%
  mutate(Model = as.character(Model)) %>%
  mutate(ModelType = case_when(
    Model == "mGRM" ~ "GRM",
    str_detect(Model, "mGRM") & str_detect(Model, "tGRM") ~ "GRM + TRM",
    !str_detect(Model, "mGRM") & str_detect(Model, "tGRM") ~ "TRM"
  )) %>%
  mutate(Trait = as.character(Trait),
         TraitCategory = as.character(TraitCategory)) %>%
  mutate(SubCategory = case_when(
    # TraitCategory == "carot" ~ NA,
    Trait %in% tocotrienols ~ "Tocotrienol",
    Trait %in% tocopherols ~ "Tocopherol",
    Trait %in% toco.summary.traits ~ "Tocochromanol summary",
    Trait %in% carotenes ~ "Carotene",
    Trait %in% xanthophylls ~ "Xanthophyll",
    Trait %in% carot.summary.traits ~ "Carotenoid summary"
  )) %>%
  dplyr::select(TraitCategory, SubCategory, Trait:Model, ModelType, Fold, PearsonCor:SpearmanRankCor) %>%
  distinct()

summarized.accuracy.df <- accuracy.df %>%
  ungroup() %>%
  dplyr::select(-Fold) %>%
  group_by(TraitCategory, SubCategory, Trait, MutantTypeIncluded, Model, ModelType) %>%
  filter(!str_detect(Model, "PH207")) %>%
  summarize(MeanPearson = mean(PearsonCor),
            MeanSpearman = mean(SpearmanRankCor))

accuracy.df.man <- summarized.accuracy.df %>%
  mutate(Model = recode_factor(Model,
                               "mGRM" = "GRM",
                               "tGRM.B73" = "TRM.B73",
                               "tGRM.Ia453" = "TRM.Ia453",
                               "tGRM.cand" = "TRM.cand",
                               "tGRM.non.cand" = "TRM.non.cand",
                               "tGRM.both" = "TRM.both",
                               "mGRM.tGRM.B73" = "GRM + TRM.B73",
                               "mGRM.tGRM.Ia453" = "GRM + TRM.Ia453",
                               "mGRM.tGRM.cand" = "GRM + TRM.cand",
                               "mGRM.tGRM.non.cand" = "GRM + TRM.non.cand",
                               "mGRM.tGRM.both" = "GRM + TRM.both"),
         SubCategory = recode(SubCategory,
                              "Carotenoid summary" = "Carotenoid ratio",
                              "Tocochromanol summary" = "Tocochromanol ratio or sum"),
         TraitCategory = recode(TraitCategory,
                                "carot" = "Carotenoid",
                                "toco" = "Tocochromanol"))

write.csv(accuracy.df.man, "./RESULT/15.1-Prediction_plots/prediction_accuracies_all_stratified_manuscript.csv", row.names = F)
write.csv(summarized.accuracy.df, "./RESULT/15.1-Prediction_plots/prediction_accuracies_strat_summarized.csv", row.names = F)

accuracy.df %>% ungroup() %>%
  filter(MutantTypeIncluded) %>%
  filter(Model %in% c("tGRM.Ia453", "tGRM.B73", "mGRM.tGRM.Ia453", "mGRM.tGRM.B73", "mGRM")) %>%
  dplyr::select(TraitCategory, SubCategory, MutantTypeIncluded,
                Trait,
                Model, ModelType, PearsonCor, SpearmanRankCor) %>%
  group_by(TraitCategory,
           MutantTypeIncluded,
           SubCategory,
           Model, ModelType) %>%
  summarize(MeanPearsonCategory = mean(PearsonCor),
            MeanSpearmanCategory = mean(SpearmanRankCor)) %>%
  arrange(TraitCategory, SubCategory, Model) %>%
  print(n = 60)

# pivot wider based on kernel type
accuracy.df %>% ungroup() %>%
  dplyr::select(TraitCategory, SubCategory, MutantTypeIncluded,
                Trait, Model, ModelType, PearsonCor) %>%
   filter(Model %in% c("tGRM.non.cand", "tGRM.cand", "tGRM.both", "mGRM.tGRM.non.cand", "mGRM.tGRM.cand", "mGRM.tGRM.both", "tGRM.B73")) %>%
  group_by(TraitCategory,
           MutantTypeIncluded,
           Model,
           ModelType) %>%
  summarize(MeanPearsonCategory = mean(PearsonCor)) %>%
  pivot_wider(id_cols = c(TraitCategory,
                          Model, ModelType),
              names_from = c(MutantTypeIncluded),
              values_from = c(MeanPearsonCategory), names_prefix = "Pearson_MutantIncluded_") %>%
  mutate(difMutantIncluded = `Pearson_MutantIncluded_TRUE` - `Pearson_MutantIncluded_FALSE`) %>%
  arrange(TraitCategory,
          -`Pearson_MutantIncluded_TRUE`
          ) %>%
  print(n = 60)

accuracy.df %>% ungroup() %>%
  dplyr::select(TraitCategory,
                Trait, MutantTypeIncluded,
                Model, ModelType, PearsonCor) %>%
  filter(Model %in% c("tGRM.non.cand", "tGRM.cand", "tGRM.both", "mGRM.tGRM.non.cand", "mGRM.tGRM.cand", "mGRM.tGRM.both", "tGRM.B73")) %>%
  group_by(TraitCategory,
           MutantTypeIncluded,
           Trait,
           Model,
           ModelType) %>%
  summarize(MeanPearsonCategory = mean(PearsonCor)) %>%
  pivot_wider(id_cols = c(TraitCategory,  Trait, MutantTypeIncluded),
              names_from = Model,
              values_from = MeanPearsonCategory) %>%
  arrange(TraitCategory, Trait) %>%
  print(n = 60)

kernel.covar.TRM.cand <- accuracy.df %>% ungroup() %>%
  dplyr::select(TraitCategory,
                Trait, MutantTypeIncluded,
                Model, ModelType, PearsonCor) %>%
  filter(Model == "tGRM.cand") %>%
  group_by(Trait, MutantTypeIncluded) %>%
  summarize(MeanPearson = mean(PearsonCor)) %>%
  pivot_wider(id_cols = Trait,
              names_from = MutantTypeIncluded,
              values_from = MeanPearson) %>%
  mutate(dif.mut = `TRUE` - `FALSE`) %>%
  print(n = 60)

write.csv(kernel.covar.TRM.cand, "RESULT/15.1-Prediction_plots/TRM.cand_kernel_type_effect.csv", row.names = F)


summarized.accuracy.df %>% filter(MutantTypeIncluded, TraitCategory == "carot") %>%
  pivot_wider(id_cols = Trait, names_from = Model, values_from = MeanPearson) %>%
  column_to_rownames(var = "Trait") %>%
  as.matrix() %>%
  cor() %>% round(2) %>%
  ggcorrplot(., hc.order = FALSE, type = "lower",
           lab = TRUE)

summarized.accuracy.df %>% filter(MutantTypeIncluded, TraitCategory == "toco") %>%
  pivot_wider(id_cols = Trait, names_from = Model, values_from = MeanPearson) %>%
  column_to_rownames(var = "Trait") %>%
  as.matrix() %>%
  cor(method = "pearson") %>% round(2) %>%
  ggcorrplot(., hc.order = FALSE, type = "upper",
             lab = TRUE)


# Zeinoxanthin test
accuracy.df %>% filter(Trait == "Zeinoxanthin") %>% filter(MutantTypeIncluded) %>%
  ggplot(aes(x = Model, y = PearsonCor, fill = ModelType)) +
  geom_boxplot() +
  labs(title = "Zeinoxanthin")


#### MANUSCRIPT ####
accuracy.df <- read.csv("./RESULT/15.1-Prediction_plots/prediction_accuracies_all_stratified_manuscript.csv")
library(ggsci)
library(tidyverse)
library(ggtext)

plotmanuscript <- accuracy.df %>%
  filter(MutantTypeIncluded) %>%
  mutate(TraitCategory = as.character(TraitCategory),
         SubCategory = as.character(SubCategory),
         plot_categories = case_when(
           TraitCategory == "Carotenoid" ~ "Carotenoids",
           SubCategory == "Tocochromanol summary" ~ NA_character_,
           SubCategory == "Tocopherol" ~ "Tocopherols",
           SubCategory == "Tocotrienol" ~ "Tocotrienols"),
         ModelType = factor(ModelType, levels = c("GRM", "TRM", "GRM + TRM"))) %>%
  drop_na(plot_categories) %>%
  mutate(plot_categories = factor(plot_categories)) %>%
  ggplot(aes(x = Model, y = MeanPearson, fill = ModelType)) +
  geom_boxplot(position = "dodge2") +
  labs(y = "Predictive Ability (*r*)")+
  theme_bw() +
  lims(y = c(0,1)) +
  facet_grid(~plot_categories) +
  scale_fill_manual(name = "Model Type", labels = c("GRM", "GRM + TRM", "TRM"),
                    values = c("#F79F1F", "#A3CB38", "#1289A7")) +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.title.y = element_markdown())

plotmanuscript
ggsave(plotmanuscript, filename = "./RESULT/Figure3_predictions_bytype_v2.png", device = "png", units = "in", width = 11, height = 7)

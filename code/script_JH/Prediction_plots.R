library(tidyverse)
library(fs)

#### Compile results into one df ####

read_csv_adapted <- function(path){
  pat_trait <- "(?<=PredictionResult_)(.*)(?=.csv)"
  trait.string <- str_extract(path, pat_trait)
  trait.vec <- unlist(str_split(trait.string, pattern = "_"))
  # [1] = trait category
  # [2] = trait
  # [3] = kernel mutant type included
  # [4] = relationship matrices
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

pred_results_dir <- "./output/8.1-Prediction"
pred_files <- list.files(pred_results_dir, pattern = "PredictionResult", full.names = T)
pred.df <- map_dfr(pred_files, read_csv_adapted)

path = "./output/8.1-Prediction/PredictionResult_carot_beta.Cryptoxanthin_over_zeaxanthin_NoMu_mGRM.tGRM.both.csv"

pred.no.ion.df <- pred.df %>%
  filter(TraitCategory != "ion")

write.csv(pred.df, "./output/8.1-Prediction_results_compiled.csv", row.names = F)
write.csv(pred.no.ion.df, "./output/8.1-Prediction_results_no_ionomics.csv", row.names = F)

#### Compute accuracies from predictions and HPLC BLUPs ####

blups_carot <- read.csv("./output/4.1-MakeBlueDatasets/PhenoData_carot.csv")
blups_toco <- read.csv("./output/4.1-MakeBlueDatasets/PhenoData_toco.csv")
#blups_ion <- read.csv("./output/output/4.1-MakeBlueDatasets/PhenoData_ion.csv")

blups_carot_long <- blups_carot %>%
  pivot_longer(cols = Antheraxanthin:Total.carotenes_over_total.xanthophylls,
               names_to = "Trait", values_to = "Observed")

blups_toco_long <- blups_toco %>%
  pivot_longer(cols = aT:ratio_TotalT_TotalT3,
               names_to = "Trait", values_to = "Observed")

blups_all <- rbind(blups_toco_long, blups_carot_long)

full <- pred.no.ion.df %>%
  pivot_longer(cols = starts_with("Rep"), names_to = "Fold", values_to = "Predicted") %>%
  full_join(blups_all, by = c("Sample.ID", "Trait")) %>%
  drop_na(Predicted, Observed) %>%
  mutate(Predicted = as.numeric(Predicted),
         Observed = as.numeric(Observed)) %>%
  group_by(Trait, Model, Fold, MutantTypeIncluded) %>%
  mutate(PearsonCor = cor(Predicted, Observed, method = "pearson"),
         SpearmanRankCor = cor(Predicted, Observed, method = "spearman"))

accuracy.df <- full %>%
  dplyr::select(TraitCategory:Fold, PearsonCor:SpearmanRankCor) %>%
  distinct()

summarized.accuracy.df <- accuracy.df %>%
  ungroup() %>%
  dplyr::select(-Fold) %>%
  group_by(Trait, MutantTypeIncluded, Model) %>%
  summarize(MeanPearson = mean(PearsonCor),
            MeanSpearman = mean(SpearmanRankCor))

toco.plot <- accuracy.df %>% filter(TraitCategory == "toco") %>%
  mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
  ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
  facet_grid("MutantTypeIncluded") +
  geom_boxplot(position = "dodge2") +
  labs(title = "Sweet Corn Tocochromanol Predictions") +
  theme_bw() +
  lims(y = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

carot.plot <- accuracy.df %>% filter(TraitCategory == "carot") %>%
  mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
  ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
  facet_grid("MutantTypeIncluded") +
  geom_boxplot(position = "dodge2") +
  labs(title = "Sweet Corn Carotenoid Predictions") +
  theme_bw() +
  lims(y = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toco.plot
carot.plot

ggsave(toco.plot, filename = "./output/Tocochromanol_predictions.png", device = "png", units = "in", width = 12, height = 9)
ggsave(carot.plot, filename = "./output/Carotenoid_predictions.png", device = "png", units = "in", width = 12, height = 9)

PH207.models <- c("mGRM.tGRM.both", "mGRM.tGRM.PH207", "tGRM.both", "tGRM.PH207")

toco.noPH207.plot <- accuracy.df %>%
  filter(TraitCategory == "toco",
         !Model %in% PH207.models) %>%
  mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
  ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
  facet_grid("MutantTypeIncluded") +
  geom_boxplot(position = "dodge2") +
  labs(title = "Sweet Corn Tocochromanol Predictions") +
  theme_bw() +
  lims(y = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

carot.noPH207.plot <- accuracy.df %>%
  filter(TraitCategory == "carot",
         !Model %in% PH207.models) %>%
  mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
  ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
  facet_grid("MutantTypeIncluded") +
  geom_boxplot(position = "dodge2") +
  labs(title = "Sweet Corn Carotenoid Predictions") +
  theme_bw() +
  lims(y = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toco.noPH207.plot
carot.noPH207.plot

ggsave(toco.noPH207.plot, filename = "./output/Tocochromanol_predictions_noPH207.png", device = "png", units = "in", width = 12, height = 9)
ggsave(carot.noPH207.plot, filename = "./output/Carotenoid_predictions_noPH207.png", device = "png", units = "in", width = 12, height = 9)



summarized.accuracy.df <- accuracy.df %>%
  filter(!Model %in% PH207.models) %>%
  ungroup() %>%
  dplyr::select(-Fold) %>%
  group_by(Trait, MutantTypeIncluded, Model) %>%
  summarize(MeanPearson = mean(PearsonCor),
            MeanSpearman = mean(SpearmanRankCor))

summarized.accuracy.df %>%
  arrange(MeanPearson) %>%
  filter(Model == "mGRM",
         MutantTypeIncluded) %>%
  head()

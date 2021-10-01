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

# # Harvest date in model
# harv_pred_results_dir <- "./output/HarvDate/8.1-Prediction"
# harv_pred_files <- list.files(harv_pred_results_dir, pattern = "PredictionResult", full.names = T)
# harv_pred.df <- map_dfr(harv_pred_files, read_csv_adapted) %>%
#   filter(TraitCategory != "ion") %>% mutate(harvest.date.in.model = TRUE)
# write.csv(harv_pred.df, "./output/HarvDate/8.1-Prediction_results_HarvDate_compiled.csv", row.names = F)

# Harvest date NOT in model
noharv_pred_results_dir <- "./output/NoHarvDate/8.1-Prediction"
noharv_pred_files <- list.files(noharv_pred_results_dir, pattern = "PredictionResult", full.names = T)
noharv_pred.df <- map_dfr(noharv_pred_files, read_csv_adapted) %>%
  filter(TraitCategory != "ion") %>% mutate(harvest.date.in.model = FALSE)
write.csv(noharv_pred.df, "./output/NoHarvDate/8.1-Prediction_results_NoHarvDate_compiled.csv", row.names = F)
#noharv_pred.df <- read.csv("./output/NoHarvDate/8.1-Prediction_results_NoHarvDate_compiled.csv")

#### Compute accuracies from predictions and HPLC BLUPs ####

blups_carot <- read.csv("./output/4.1-MakeBlueDatasets/PhenoData_carot.csv")
blups_toco <- read.csv("./output/4.1-MakeBlueDatasets/PhenoData_toco.csv")
#blups_ion <- read.csv("./output/4.1-MakeBlueDatasets/PhenoData_ion.csv")

blups_carot_long <- blups_carot %>%
  pivot_longer(cols = Antheraxanthin:Total.carotenes_over_total.xanthophylls,
               names_to = "Trait", values_to = "Observed")

blups_toco_long <- blups_toco %>%
  pivot_longer(cols = aT:ratio_TotalT_TotalT3,
               names_to = "Trait", values_to = "Observed")

blups_all <- rbind(blups_toco_long, blups_carot_long)

PH207.models <- c("mGRM.tGRM.both", "mGRM.tGRM.PH207", "tGRM.both", "tGRM.PH207")

full <- noharv_pred.df %>% #rbind(harv_pred.df, noharv_pred.df) %>%
  filter(!Model %in% PH207.models) %>%
  pivot_longer(cols = starts_with("Rep"), names_to = "Fold", values_to = "Predicted") %>%
  full_join(blups_all, by = c("Sample.ID", "Trait")) %>%
  drop_na(Predicted, Observed) %>%
  mutate(Predicted = as.numeric(Predicted),
         Observed = as.numeric(Observed)) %>%
  group_by(Trait, Model, Fold, MutantTypeIncluded) %>% #, harvest.date.in.model) %>%
  mutate(PearsonCor = cor(Predicted, Observed, method = "pearson"),
         SpearmanRankCor = cor(Predicted, Observed, method = "spearman"))

accuracy.df <- full %>%
  dplyr::select(#harvest.date.in.model,
                TraitCategory:Fold, PearsonCor:SpearmanRankCor) %>%
  distinct()

summarized.accuracy.df <- accuracy.df %>%
  ungroup() %>%
  dplyr::select(-Fold) %>%
  group_by(Trait, MutantTypeIncluded, Model#, harvest.date.in.model
           ) %>%
  summarize(MeanPearson = mean(PearsonCor),
            MeanSpearman = mean(SpearmanRankCor))

write.csv(accuracy.df, "./output/NoHarvDate/prediction_accuracies_all.csv", row.names = F)

toco.plot <- accuracy.df %>% filter(TraitCategory == "toco") %>%
  mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
  #mutate(harvest.date.in.model = ifelse(harvest.date.in.model, "Harvest date included", "Harvest date NOT included")) %>%
  ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
  #facet_grid(MutantTypeIncluded ~ harvest.date.in.model) +
  facet_grid("MutantTypeIncluded") +
  geom_boxplot(position = "dodge2") +
  labs(title = "Sweet Corn Tocochromanol Predictions") +
  theme_bw() +
  lims(y = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

carot.plot <- accuracy.df %>% filter(TraitCategory == "carot") %>%
  mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
  ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
  #facet_grid(MutantTypeIncluded ~ harvest.date.in.model) +
  facet_grid("MutantTypeIncluded") +
  geom_boxplot(position = "dodge2") +
  labs(title = "Sweet Corn Carotenoid Predictions") +
  theme_bw() +
  lims(y = c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toco.plot
carot.plot

ggsave(toco.plot, filename = "./output/NoHarvDate/Tocochromanol_predictions.png", device = "png", units = "in", width = 12, height = 12)
ggsave(carot.plot, filename = "./output/NoHarvDate/Carotenoid_predictions.png", device = "png", units = "in", width = 12, height = 12)


#### EXIT SEMINAR ####

traits.carot <- accuracy.df %>% ungroup() %>%  filter(TraitCategory == "carot") %>% dplyr::select(Trait) %>% unique() %>% unlist()
traits.toco <- accuracy.df %>% ungroup() %>%  filter(TraitCategory == "toco") %>% dplyr::select(Trait) %>% unique() %>% unlist()
ratios <- c(traits.carot[c(4:5, 7,9,13,19)], traits.toco[c(4,6,8,10,11,12,13:17)])

accuracy.df

toco.plot.seminar <- accuracy.df %>%
  filter(TraitCategory == "toco",
         !Trait %in% ratios,
         MutantTypeIncluded) %>%
  # mutate(Trait = recode(Trait, "Total.T3_plus_T" = "Total Tocochromanols",
  #                       "Total.T" = "Total Tocopherols",
  #                       "Total.T3" = "Total Tocotrienols")) %>%
  #mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
  #mutate(harvest.date.in.model = ifelse(harvest.date.in.model, "Harvest date included", "Harvest date NOT included")) %>%
  ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
  #facet_grid(MutantTypeIncluded ~ harvest.date.in.model) +
  #facet_grid("MutantTypeIncluded") +
  geom_boxplot(position = "dodge2") +
  labs(y = "Accuracy")+#title = "Sweet Corn Tocochromanol Predictions") +
  theme_bw() +
  lims(y = c(0,1)) +
  scale_fill_discrete(name = "Relationship matrices", labels = c("Genomic only", "Genomic + Transcriptomic", "Transcriptomic only"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

carot.plot.seminar <- accuracy.df %>%
  filter(TraitCategory == "carot",
         !Trait %in% ratios,
         MutantTypeIncluded) %>%
  filter(MutantTypeIncluded) %>%
  #mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
  ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
  #facet_grid(MutantTypeIncluded ~ harvest.date.in.model) +
  #facet_grid("MutantTypeIncluded") +
  geom_boxplot(position = "dodge2") +
  labs(y = "Accuracy")+#title = "Sweet Corn Carotenoid Predictions") +
  theme_bw() +
  lims(y = c(0,1)) +
  scale_fill_discrete(name = "Relationship matrices", labels = c("Genomic only", "Genomic + Transcriptomic", "Transcriptomic only"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

toco.plot.seminar
carot.plot.seminar

ggsave(toco.plot.seminar, filename = "./output/NoHarvDate/Tocochromanol_predictions_seminar.png", device = "png", units = "in", width = 9, height = 5)
ggsave(carot.plot.seminar, filename = "./output/NoHarvDate/Carotenoid_predictions_seminar.png", device = "png", units = "in", width = 9, height = 5)












# PH207.models <- c("mGRM.tGRM.both", "mGRM.tGRM.PH207", "tGRM.both", "tGRM.PH207")
#
# toco.noPH207.plot <- accuracy.df %>%
#   filter(TraitCategory == "toco",
#          !Model %in% PH207.models) %>%
#   mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
#   ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
#   facet_grid(MutantTypeIncluded~harvest.date.in.model) +
#   geom_boxplot(position = "dodge2") +
#   labs(title = "Sweet Corn Tocochromanol Predictions") +
#   theme_bw() +
#   lims(y = c(0,1)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
# toco.noPH207.plot_color <- accuracy.df %>%
#   filter(TraitCategory == "toco",
#          !Model %in% PH207.models) %>%
#   mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
#   ggplot(aes(x = Trait, y = PearsonCor, fill = Model, color = harvest.date.in.model)) +
#   facet_grid("MutantTypeIncluded") +
#   geom_boxplot(position = "dodge2") +
#   scale_color_manual(values = c("#999999", "#E69F00"))+
#   labs(title = "Sweet Corn Tocochromanol Predictions") +
#   theme_bw() +
#   lims(y = c(0,1)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
#
# carot.noPH207.plot <- accuracy.df %>%
#   filter(TraitCategory == "carot",
#          !Model %in% PH207.models) %>%
#   mutate(MutantTypeIncluded = ifelse(MutantTypeIncluded, "Kernel mutant type included", "Kernel mutant type NOT included")) %>%
#   ggplot(aes(x = Trait, y = PearsonCor, fill = Model)) +
#   facet_grid(MutantTypeIncluded~harvest.date.in.model) +
#   geom_boxplot(position = "dodge2") +
#   labs(title = "Sweet Corn Carotenoid Predictions") +
#   theme_bw() +
#   lims(y = c(0,1)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
# toco.noPH207.plot
# carot.noPH207.plot
#
# ggsave(toco.noPH207.plot, filename = "./output/Tocochromanol_predictions_noPH207.png", device = "png", units = "in", width = 12, height = 12)
# ggsave(carot.noPH207.plot, filename = "./output/Carotenoid_predictions_noPH207.png", device = "png", units = "in", width = 12, height = 12)
#
#
#
# summarized.accuracy.df <- accuracy.df %>%
#   filter(!Model %in% PH207.models) %>%
#   ungroup() %>%
#   dplyr::select(-Fold) %>%
#   group_by(Trait, MutantTypeIncluded, Model) %>%
#   summarize(MeanPearson = mean(PearsonCor),
#             MeanSpearman = mean(SpearmanRankCor))
#
# summarized.accuracy.df %>%
#   arrange(MeanPearson) %>%
#   filter(Model == "mGRM",
#          MutantTypeIncluded) %>%
#   head()

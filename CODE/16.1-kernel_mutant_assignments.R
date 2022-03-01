# Assign kernel mutant types based on line name

library(tidyverse)
library(readxl)
library(readr)
#remotes::install_github("clauswilke/ggtext")
library(ggtext)

Baseggio_assignments_2020 <- read_excel("RAWDATA/Carotenoids/carotenoid_transformed_BLUPs-converted.xlsx", skip = 1) %>%
  janitor::clean_names() %>%
  mutate(Baseggio_assignment = as.character(endosperm_mutation)) %>%
  dplyr::select(sample_id:gbs_id, Baseggio_assignment) %>%
  rename_at(vars(-gbs_id), function(x) paste0(x,"_carot"))

Baseggio_assignments <- read_csv("RAWDATA/Tocochromanols/toco_back_transformed_blups.csv") %>%
  janitor::clean_names() %>%
  mutate(Baseggio_assignment = as.character(endosperm_mutation)) %>%
  dplyr::select(sample_id:gbs_id, classification_method, Baseggio_assignment) %>%
  rename_at(vars(-gbs_id, -classification_method), function(x) paste0(x,"_toco")) %>%
  full_join(Baseggio_assignments_2020) %>%
  mutate(sample_id = case_when(
    !is.na(sample_id_toco) ~ sample_id_toco,
    is.na(sample_id_toco) & !is.na(sample_id_carot) ~ sample_id_carot),
    Baseggio_assignment = case_when(
      !is.na(Baseggio_assignment_toco) ~ Baseggio_assignment_toco,
      is.na(Baseggio_assignment_toco) & !is.na(Baseggio_assignment_carot) ~ Baseggio_assignment_carot)) %>%
  dplyr::select(sample_id, Baseggio_assignment, classification_method) %>%
  drop_na() %>%
  distinct()


Hersh_assignments <- read_csv("RAWDATA/Metadata/Hershberger_assignments.csv") %>%
  drop_na(Hersh_assignment) %>%
  rename(sample_id = VarietyName) %>%
  dplyr::select(sample_id, Hersh_assignment)

assignments <- Baseggio_assignments %>%
  left_join(Hersh_assignments) %>%
  mutate(Baseggio_assignment = case_when(
    sample_id == "304A" ~ "su1",
    sample_id == "675A" ~ "su1",
    TRUE ~ Baseggio_assignment)) %>%
  drop_na(Baseggio_assignment) %>%
  dplyr::select(sample_id, Hersh_assignment, Baseggio_assignment,
                classification_method) %>%
  mutate(sample_id = recode(sample_id,
                            C8_pseudo = "C8p",
                            C23_su = "C23",
                            IaEv191 = "IaEV191",
                            IL101T = "Il101t",
                            IL677a = "Il677a",
                            Il731a = "IL731a",
                            Il767b = "IL767b",
                            Il779a = "IL779a",
                            IL802b_Ht2 = "IL802b",
                            Luther_Hill = "LUTHER_HILL",
                            Me244_Wb = "Me244_wb",
                            Olcott = "OLCOTT",
                            P39A = "P39a",
                            P39_Goodman_Buckler = "P39Goodman_Buckler",
                            P39Le_253_60 = "P39Le_253",
                            Strain_T20_2_68B = "T20",
                            Strain_T24_395_68B = "T24",
                            Strain_T32_397_68r = "T32",
                            Strain_T33_399_68B = "T33",
                            Strain_T35_388_68A = "T35",
                            T62S = "T62s"
                            ))

write.csv(assignments, "RESULT/16.1-kernel_mutant_assignments/kernel_assignments_20210923.csv", row.names = F, na = "")

rlog_b73 <- read.csv("RAWDATA/SweetCorn_TagSeq/htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.csv") %>%
  filter(gene_id %in% c("Zm00001d044129", "Zm00001d049753")) %>%
  dplyr::select(gene_id, starts_with("SC")) %>%
  pivot_longer(cols = starts_with("SC"), names_to = "sample", values_to = "rlog_count") %>%
  mutate(plot_sample_id = str_extract(sample, pattern = "(?<=19A).*")) %>%
  mutate(sample_id = str_extract(plot_sample_id, pattern = "(?<=_).*" )) %>%
  dplyr::select(-plot_sample_id) %>%
  pivot_wider(id_cols = c(sample, sample_id), names_from = gene_id, values_from = rlog_count) %>%
  mutate(sample_id = str_replace(sample_id, pattern = "\\.", replacement = "_"))

rlog_b73_assigned <- rlog_b73 %>%
  full_join(assignments)

write.csv(rlog_b73_assigned, "RESULT/16.1-kernel_mutant_assignments/kernel_mutants_rlog_B73.csv", row.names = F)
#rlog_b73_assigned <- read.csv("RESULT/16.1-kernel_mutant_assignments/kernel_mutants_rlog_B73.csv")

sh2i_b73 <- rlog_b73_assigned %>%
  drop_na(Baseggio_assignment, classification_method) %>%
  rename(`Classification Method` = classification_method) %>%
  mutate(`Name-based\nclassification` = ifelse(is.na(Hersh_assignment), "Sh2 or unknown",
                                               as.character(Hersh_assignment))) %>%
  pivot_longer(cols = starts_with("Z"), names_to = "Gene", values_to = "rlog_count") %>%
  filter(Gene %in% c("Zm00001d044129", "Zm00001d049753")) %>%
  mutate(GeneLab = case_when(
    Gene == "Zm00001d044129" ~ "*shrunken2* (Zm00001d044129)",
    Gene == "Zm00001d049753" ~ "*sugary1* (Zm00001d049753)",
  )) %>%
  mutate(Baseggio_assignment = case_when(
    Baseggio_assignment == "sh2" ~ "Su1sh2",
    Baseggio_assignment == "su1" ~ "su1Sh2",
    Baseggio_assignment == "su1sh2" ~ "su1sh2",
  )) %>%
  ggplot(aes(y = rlog_count, x = Baseggio_assignment)) +
  geom_violin() +
  geom_jitter(aes(color = `Name-based\nclassification`,
                  shape = `Name-based\nclassification`),
             alpha = 0.7) +
  labs(x = "Visual and marker-based classification", y = "Normalized transcript abundance") +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     labels = c("*su1Sh2* or unknown", "*Su1sh2-R*", "*su1sh2-i*", "*su1sh2-R*"))+
  scale_shape_manual(values = 15:18,
                     labels = c("*su1Sh2* or unknown", "*Su1sh2-R*", "*su1sh2-i*", "*su1sh2-R*"))+
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.text = element_markdown(),
        strip.text.x = element_markdown())+
  facet_grid(~GeneLab)
sh2i_b73
ggsave(plot = sh2i_b73, filename = "RESULT/16.1-kernel_mutant_assignments/Figure2_highres.png",
       device = "png", units = "in", width = 7, height = 5, dpi = 400)


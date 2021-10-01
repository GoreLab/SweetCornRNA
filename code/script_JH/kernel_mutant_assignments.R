library(tidyverse)
library(readxl)
library(readr)
#remotes::install_github("clauswilke/ggtext")
library(ggtext)

# read.csv not working with this file for some reason so manually imported as df
Hislop_assignments <- read_csv("data/Hislop_assignments.csv") %>%
  dplyr::select(VarietyName, GenoName:Program) %>%
  rename(sample_id = VarietyName,
         Hislop_assignment = Endospermtype) %>%
  drop_na(Hislop_assignment) %>%
  mutate(sample_id = str_replace(sample_id, "-", "_")) %>%
  mutate(sample_id = recode(sample_id,
                            C8pseudo = "C8_pseudo",
                            C68R = "C68Rx552_63",
                            `Fa56A` = "Fa56a",
                            IaEV191 = "IaEv191",
                            IL101t = "IL101T",
                            IL14H = "Il14H",
                            IL731 = "Il731a",
                            IL767b = "Il767b",
                            IL779a = "Il779a",
                            `IL802b Ht2` = "IL802b_Ht2",
                            `LUTHER HILL` = "Luther_Hill",
                            Me244_wb = "Me244_Wb",
                            OLCOTT = "Olcott",
                            `P39 Goodman_Buckler` = "P39_Goodman_Buckler",
                            `P39(New)` = "P39_New",
                            P39a = "P39A",
                            `P39Le_253:60` = "P39Le_253_60",
                            `P39LE_261:61` = "P39LE_261_61",
                            `STRAIN 304A_408-68B` = "Strain_304A_408_68B",
                            `STRAIN 675A_415-68(B)` = "Strain_675A_415_68B",
                            `STRAIN T20_2-68B` = "Strain_T20_2_68B",
                            `STRAIN T24_395-68B` = "Strain_T24_395_68B",
                            `STRAIN T32_397-68B` = "Strain_T32_397_68r",
                            `STRAIN T33_399-68B` = "Strain_T33_399_68B",
                            `STRAIN T35_388-68A` = "Strain_T35_388_68A",
                            T62s = "T62S"
                            ))

Baseggio_assignments_2020 <- read_excel("data/phenotypes/carotenoid_transformed_BLUPs-converted.xlsx", skip = 1) %>%
  janitor::clean_names() %>%
  mutate(Baseggio_assignment = as.character(endosperm_mutation)) %>%
  dplyr::select(sample_id:gbs_id, Baseggio_assignment) %>%
  rename_at(vars(-gbs_id), function(x) paste0(x,"_carot"))

Baseggio_assignments <- read_csv("data/phenotypes/toco_back_transformed_blups.csv") %>%
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


Hersh_assignments <- read_csv("data/Hershberger_assignments.csv") %>%
  drop_na(Hersh_assignment) %>%
  rename(sample_id = VarietyName) %>%
  dplyr::select(sample_id, Hersh_assignment)

assignments <- full_join(Hislop_assignments, Baseggio_assignments) %>%
  left_join(Hersh_assignments) %>%
  mutate(Baseggio_assignment = case_when(
    sample_id == "304A" ~ "su1",
    sample_id == "675A" ~ "su1",
    TRUE ~ Baseggio_assignment)) %>%
  drop_na(Baseggio_assignment) %>%
  mutate(assignment_matches = Hislop_assignment == Baseggio_assignment) %>%
  dplyr::select(sample_id, Hersh_assignment, Hislop_assignment, Baseggio_assignment,
                assignment_matches, classification_method) %>%
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


mismatch_assignments <- assignments %>%
  filter(!assignment_matches & Hislop_assignment != "sh2i" & Hislop_assignment != "se") %>% print(n = 100)

write.csv(assignments, "output/kernel_assignments_20210923.csv", row.names = F, na = "")

rlog_b73 <- read.csv("data/htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.csv") %>%
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

write.csv(rlog_b73_assigned, "output/kernel_mutants_rlog_B73.csv", row.names = F)
rlog_b73_assigned <- read.csv("output/kernel_mutants_rlog_B73.csv")

sh2i_b73 <- rlog_b73_assigned %>%
  drop_na(Baseggio_assignment, classification_method) %>%
  rename(`Classification Method` = classification_method,
         `Hislop assignment` = Hislop_assignment) %>%
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
  #geom_jitter(aes(color = `Classification Method`))+
  geom_jitter(aes(color = `Name-based\nclassification`,
                  shape = `Name-based\nclassification`),
             alpha = 0.7
              )+
  #geom_jitter(aes(color = `Hislop assignment`, shape = `Classification Method`))+
  #labs(title = "Zm00001d044129", subtitle = "shrunken2 gene, B73 alignment") +
  labs(x = "Visual and marker-based classification", y = "Normalized transcript abundance") +
  # labs(y = "Baseggio kernel endosperm mutation type assignment", x = "Normalized transcript abundance") +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     labels = c("*su1Sh2* or unknown", "*Su1sh2-R*", "*su1sh2-i*", "*su1sh2-R*")
                     )+
  scale_shape_manual(values = 15:18,
                     labels = c("*su1Sh2* or unknown", "*Su1sh2-R*", "*su1sh2-i*", "*su1sh2-R*")
                     )+
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.text = element_markdown(),
        strip.text.x = element_markdown())+
  facet_grid(~GeneLab)
sh2i_b73
ggsave(plot = sh2i_b73, filename = "output/kernel_assignments/sh2i_B73alignment_vertical.png",
       device = "png", units = "in", width = 7, height = 5)



sh2i_b73_seminar <- rlog_b73_assigned %>%
  drop_na(Baseggio_assignment, classification_method) %>%
  rename(`Classification Method` = classification_method,
         `Hislop assignment` = Hislop_assignment) %>%
  #mutate("*shrunken2* allele" = ifelse(is.na(Hersh_assignment), "Sh2 or unknown", Hersh_assignment)) %>%
  pivot_longer(cols = starts_with("Z"), names_to = "Gene", values_to = "rlog_count") %>%
  filter(Gene %in% c("Zm00001d044129", "Zm00001d049753")) %>%
  mutate(GeneLab = case_when(
    Gene == "Zm00001d044129" ~ "*shrunken2*",
    Gene == "Zm00001d049753" ~ "*sugary1*",
  )) %>%
  mutate(`*shrunken2* allele` = case_when(
    str_detect(Hersh_assignment, "sh2-i") ~ "sh2-i",
    str_detect(Hersh_assignment, "sh2-R") ~ "sh2-R",
    is.na(Hersh_assignment) ~ "Sh2 or unknown"
  )) %>%
  ggplot(aes(y = rlog_count, x = Baseggio_assignment)) +
  geom_violin() +
  geom_jitter(aes(color = `*shrunken2* allele`), alpha = 0.7)+
  labs(x = "Visual and marker-based classification", y = "Normalized transcript abundance") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"),
                     labels = c("*Sh2* or unknown", "*sh2-i*", "*sh2-R*"))+
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.text = element_markdown(),
        legend.title = element_markdown(),
        # legend.position = "bottom",
        strip.text.y = element_markdown())+
  facet_grid(~GeneLab)
sh2i_b73_seminar
ggsave(plot = sh2i_b73_seminar, filename = "output/kernel_assignments/sh2i_B73alignment_seminar.png",
       device = "png", units = "in", width = 6, height = 4)




rlog_b73_assigned <- read.csv("output/kernel_mutants_rlog_B73.csv")
sh2_b73 <- rlog_b73_assigned %>%
  drop_na(Baseggio_assignment, classification_method) %>%
  rename(`Classification Method` = classification_method) %>%
  pivot_longer(cols = starts_with("Z"), names_to = "Gene", values_to = "rlog_count") %>%
  filter(Gene == "Zm00001d044129") %>%
  ggplot(aes(x = rlog_count, y = Baseggio_assignment)) +
  geom_violin() +
  geom_jitter(aes(color = `Classification Method`))+
  #geom_jitter(aes(color = Hislop_assignments, shape = `Baseggio Classification Method`))+
  #labs(title = "Zm00001d044129", subtitle = "shrunken2 gene, B73 alignment") +
  labs(y = "Kernel endosperm mutation type assignment", x = "Normalized transcript abundance") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
#facet_grid(~Gene)
sh2_b73
ggsave(plot = sh2_b73, filename = "output/kernel_assignments/sh2_B73alignment.png",
       device = "png", units = "in", width = 6, height = 5)


sh2su1_b73 <- rlog_b73_assigned %>%
  drop_na(Baseggio_assignment, classification_method) %>%
  rename(`Classification Method` = classification_method) %>%
  pivot_longer(cols = starts_with("Z"), names_to = "Gene", values_to = "rlog_count") %>%
  filter(Gene %in% c("Zm00001d044129", "Zm00001d049753")) %>%
  ggplot(aes(x = rlog_count, y = Baseggio_assignment)) +
  geom_violin() +
  geom_jitter(aes(color = `Classification Method`), alpha = 0.7) +
  labs(y = "Kernel endosperm mutation type assignment", x = "Normalized transcript abundance") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic")) +
  facet_grid(~Gene)
sh2su1_b73
ggsave(plot = sh2su1_b73, filename = "output/kernel_assignments/su1sh2_B73alignment.png",
       device = "png", units = "in", width = 7, height = 4)







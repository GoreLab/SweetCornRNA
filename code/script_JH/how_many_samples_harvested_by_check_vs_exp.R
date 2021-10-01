crosses_2019 <- read_excel("data/metadata/crosses_2019.xlsx", na = "NA")

harvested_crosses <- crosses_2019 %>% drop_na(harvested)
length(unique(harvested_crosses$plot_id))

fieldbook <- read.csv("data/metadata/rnaseq_trial_2019_upload.csv") %>% rename(plot_id = plot_name)

keyfile <- harvested_crosses %>% left_join(fieldbook) %>% dplyr::select(plot_id, is_a_control)
table(keyfile$is_a_control)

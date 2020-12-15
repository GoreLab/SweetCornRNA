#!/usr/bin/env Rscript

# RNA extraction QC
# Jenna Hershberger
# jmh579@cornell.edu
# 07/14/2020

# required input:
# args[1] = working directory ex// "/Users/Jenna/Documents/GitHub/SweetCornRNA"


# example:
# Rscript rna_qc_app.R "/Users/Jenna/Documents/GitHub/SweetCornRNA"
library(tidyverse)
library(googledrive)
library(readxl)
library(lubridate)

# for debugging
#args = "/Users/Jenna/Documents/GitHub/SweetCornRNA"

args = commandArgs(trailingOnly=TRUE)

# test if there are three arguments: if not, return an error
if (length(args)<1) {
  stop("An input argument is required: Working directory", call.=FALSE)
}

#working.directory <- args[1]
working.directory <- "/Users/Jenna/Documents/GitHub/SweetCornRNA"
#working.directory <- getwd()
setwd(working.directory) #"/Users/Jenna/Documents/GitHub/SweetCornRNA"


#### GRINDING LOG ####
drive_download(file = "SC_grinding_2020", overwrite = T,
               path = paste0(working.directory, "/data/SC_grinding_2020.xlsx"))
grinding_log <- read_excel("./data/SC_grinding_2020.xlsx", na = "NA") %>%
  separate(Polycon_barcode, into = c("Sample_barcode", "rep"), sep = "_", remove = FALSE) %>%
  mutate(Sample_barcode = case_when(Sample_barcode == "19A0012" & Notes == "from sample harvested 8/22" ~ Polycon_barcode,
                                    Sample_barcode == "19A0012" & Notes != "from sample harvested 8/22" ~ "19A0012_2",
                                    Sample_barcode == "19A0096" & Notes == "from sample harvested 9/6/2019" ~ Polycon_barcode,
                                    Sample_barcode == "19A0096" & Notes == "from sample harvested 9/11/2019" ~ "19A0096_2",
                                    Sample_barcode == "19A0059" & Notes != "sample harvested 9/20" ~ Polycon_barcode,
                                    Sample_barcode == "19A0059" & Notes == "sample harvested 9/20" ~ "19A0059_2",
                                    TRUE ~ Sample_barcode)) %>%
  rename(grind_notes = Notes) %>%
  filter(!is.na(`5ml_Rack_Number`)) %>%
  dplyr::select(Sample_barcode, Grind_Date, `5ml_Rack_Number`, grind_notes)

#### EXTRACTIONS LOG ####
drive_download(file = "rna_extractions_log.xlsx", overwrite = T,
               path = paste0(working.directory, "/data/rna_extractions_log.xlsx"))
raw_firstpass <- read_excel("./data/rna_extractions_log.xlsx", sheet = "FirstPass", col_names = TRUE,
                            range = cell_cols("A:H"), na = "NA") %>%
  rename(extraction_notes = Notes)
raw_QC <- read_excel("./data/rna_extractions_log.xlsx", sheet = "QC", col_names = TRUE,
                     range = cell_cols("A:P"), na = "NA") %>% rename(Sample_barcode = Name)
raw_gel <- read_excel("./data/rna_extractions_log.xlsx", sheet = "Gels", col_names = TRUE,
                     range = cell_cols("A:G"), na = "NA") %>%
  filter(Gel_num != 1 & Gel_num != 3)


#### MATCH AND SUMMARIZE ####
matched <- raw_firstpass %>%
  full_join(raw_QC, by = "Sample_barcode") %>%
  full_join(raw_gel, by = "Sample_barcode") %>%
  full_join(grinding_log, by = "Sample_barcode")

summary_df <- matched %>%
  dplyr::select(Sample_barcode, Extraction_day_1_date, done_box_number, extraction_notes,
                QC_date, `260/280`:Gel_date, gel_pass, gel_notes, Grind_Date, `5ml_Rack_Number`,
                grind_notes)
redos <- summary_df %>%
  filter(gel_pass == FALSE | QC_pass == FALSE) %>%
  mutate(qc_matches = ifelse(gel_pass == QC_pass, TRUE, FALSE))

redos.sorted <- redos %>%
  arrange(`5ml_Rack_Number`) %>%
  dplyr::select(Sample_barcode, `5ml_Rack_Number`, grind_notes, extraction_notes) %>%
  filter(Sample_barcode != "19A0012_1") %>%
  filter(Sample_barcode != "19A0096_1")
#write.csv(redos.sorted)
#write.csv(redos.sorted, paste0(working.directory, "/output/rna_qc_sorted_redos", today(), ".csv"), row.names = F)

redos.done <- matched %>%
  filter(extraction_notes == "REDO") %>%
  rename(Redo_barcode = Sample_barcode,
         Redo_extraction_day_1 = Extraction_day_1_date,
         Redo_tube_number = Tube_number,
         Redo_box_number = Box_number,
         Redo_notes = extraction_notes,
         Redo_extraction_day_2 = Extraction_day_2_date,
         Redo_dnase_date = dnase_date,
         Redo_done_box_number = done_box_number,
         Redo_QC_date = QC_date,
         Redo_QC_pass = QC_pass) %>%
  dplyr::select(starts_with("Redo"))
redos.done.report <- redos %>%
  mutate(Redo_barcode = paste0(Sample_barcode, "_2")) %>%
  full_join(redos.done) %>%
  filter(Sample_barcode != "19A0012_1") %>%
  filter(Sample_barcode != "19A0096_1") %>%
  mutate(Redo_tissue_to_tube = "7/16/2020") %>%
  dplyr::select(Sample_barcode:qc_matches, Redo_tissue_to_tube, Redo_extraction_day_1,
                Redo_barcode:Redo_QC_pass)
#Redo_tissue_to_tube	Redo_extraction_day_1	Redo_barcode	Redo_tube_number	Redo_box_number	Redo_notes	Redo_extraction_day_2	Redo_dnase_date	Redo_done_box_number	Redo_QC_date	Redo_QC_pass
#write.csv(redos, "redos.csv", row.names = F)
write.csv(redos.done.report, paste0(working.directory, "/output/rna_qc_redos", today(), ".csv"), row.names = F,
          na = "")
write.csv(summary_df, paste0(working.directory, "/output/rna_qc_summary_", today(), ".csv"), row.names = F)

plating <- matched %>% separate(Sample_barcode, into = c("Sample_barcode", "rep"), sep = "_",
                        remove = TRUE) %>% mutate(rep = ifelse(is.na(rep), 1, rep)) %>%
  dplyr::select(Sample_barcode, rep, done_box_number, Extraction_day_1_date, Tube_number, QC_pass, starts_with("ng"))

plating.nondups <- plating %>% group_by(Sample_barcode) %>%
  filter(n()==1)

plating.dups <- plating %>% group_by(Sample_barcode) %>%
  filter(n()>1) %>%
  # all duplicates should just be rep 2 except for 59, which both reps passed but rep 2 looks better QC-wise
  filter(rep == 1 & Sample_barcode == "19A0059" | rep == 2 & Sample_barcode != "19A0059") %>%
  filter(Sample_barcode != "19A0254" & Sample_barcode != "19A0276") # these redos had no RNA

get_dilution <- function(concentration){
  if(concentration < 250 & concentration > 100){
    return(0) # dilution.H2O <- 0
  } else if(concentration > 250){
    sample.vol <- 4000 / concentration
    dilution.H2O <- 20 - sample.vol
    return(floor(dilution.H2O)) #
  } else{
    return(NA)
  }
}


set.seed(438)
plate.order.vec <- sample(x = 438, replace = FALSE)
plating.complete <- rbind(plating.nondups, plating.dups) %>%
  arrange(Extraction_day_1_date, Tube_number) %>%
  rename(initial.conc = starts_with("ng")) %>%
  mutate(dilution.H2O = unlist(map(.f = get_dilution, . = initial.conc)),
         dilution.sample = 20 - dilution.H2O,
         final.conc = initial.conc * dilution.sample / 20) %>%
  rename(Sample.Name = Sample_barcode,
         Concentration.ng.ul. = final.conc)

plating.complete$plate.order <- plate.order.vec


write.csv(plating.complete, paste0("./output/plating_", today(), ".csv"), row.names = F)

#### PLATES ####
# Need 2 blanks and 2 positive controls per plate

# One negative control (blank) specified by sequencing facility
# Need to select 3/96 wells for other controls
# Templates are mismatched from online form because they were generated before putting in 96 wells and it won't let me regenerate.
# Need to fix with updated blank wells below
plate.template <- read.csv("~/Desktop/RNA_Extraction/submission/SC_RNA_01Template.csv",
                             skip = 1, header = T, colClasses = "character")
plate.template$Within.plate.num <- 1:96

plate.template.master <- rbind(plate.template, plate.template, plate.template, plate.template, plate.template) %>%
  rownames_to_column() %>% dplyr::select(everything(), rowname)
plate.template.master$Plate.Name <- rep(c("SC_RNA_01", "SC_RNA_02", "SC_RNA_03", "SC_RNA_04", "SC_RNA_05"), c(96,96,96,96,96))


#plate.template$row.num <- 1:96

blank.row <- c(1, "BLANK", 0, 0, "Neg_Control", NA)
neg.control.row <- c(1, "Neg_Control", 0, 20, "Neg_Control", NA)
pos.control.row <- c(4577, "Pos_Control", 200, 20, "Pos_Control", "Kernel")

set.seed(1)
control.wells <- sample(1:96, 15, replace = F)

# controls.df <- as.data.frame(matrix(nrow = 20, ncol = 6))
# colnames(controls.df) <- colnames(plate.template.master)[4:9]
# controls.df$Plate.Name <- rep(c("SC_RNA_01", "SC_RNA_02", "SC_RNA_03", "SC_RNA_04", "SC_RNA_05"), c(4,4,4,4,4))
#
# controls.df[1, 6:9]



# Plate 1
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_01" &
                              plate.template.master$Well == "C04"), 6:11] <- blank.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_01" &
                              plate.template.master$Within.plate.num == control.wells[1]), 6:11] <- neg.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_01" &
                              plate.template.master$Within.plate.num == control.wells[2]), 6:11] <- pos.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_01" &
                              plate.template.master$Within.plate.num == control.wells[3]), 6:11] <- pos.control.row

# Plate 2
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_02" &
                              plate.template.master$Well == "G08"), 6:11] <- blank.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_02" &
                              plate.template.master$Within.plate.num == control.wells[4]), 6:11] <- neg.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_02" &
                              plate.template.master$Within.plate.num == control.wells[5]), 6:11] <- pos.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_02" &
                              plate.template.master$Within.plate.num == control.wells[6]), 6:11] <- pos.control.row

# Plate 3
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_03" &
                              plate.template.master$Well == "D08"), 6:11] <- blank.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_03" &
                              plate.template.master$Within.plate.num == control.wells[7]), 6:11] <- neg.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_03" &
                              plate.template.master$Within.plate.num == control.wells[8]), 6:11] <- pos.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_03" &
                              plate.template.master$Within.plate.num == control.wells[9]), 6:11] <- pos.control.row


# Plate 4
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_04" &
                              plate.template.master$Well == "D07"), 6:11] <- blank.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_04" &
                              plate.template.master$Within.plate.num == control.wells[10]), 6:11] <- neg.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_04" &
                              plate.template.master$Within.plate.num == control.wells[11]), 6:11] <- pos.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_04" &
                              plate.template.master$Within.plate.num == control.wells[12]), 6:11] <- pos.control.row


# Plate 5
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_05" &
                              plate.template.master$Well == "G04"), 6:11] <- blank.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_05" &
                              plate.template.master$Within.plate.num == control.wells[13]), 6:11] <- neg.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_05" &
                              plate.template.master$Within.plate.num == control.wells[14]), 6:11] <- pos.control.row
plate.template.master[which(plate.template.master$Plate.Name == "SC_RNA_05" &
                              plate.template.master$Within.plate.num == control.wells[15]), 6:11] <- pos.control.row

#### Populate with samples! ####
samples.master <- plate.template.master %>%
  filter(Sample.Name != "BLANK" & Sample.Name != "Pos_Control" & Sample.Name != "Neg_Control")
samples.master$plate.order <- 1:nrow(samples.master)

controls.master <- plate.template.master %>%
  filter(Sample.Name == "BLANK" | Sample.Name == "Pos_Control" | Sample.Name == "Neg_Control") %>%
  mutate(plate.order = NA,
         rep = NA,
         done_box_number = NA,
         Extraction_day_1_date = NA,
         Tube_number = NA,
         QC_pass = NA,
         initial.conc = NA,
         dilution.H2O = NA,
         dilution.sample = NA)

plate.master <- samples.master %>% dplyr::select(-Sample.Name, -Concentration.ng.ul.) %>%
  full_join(plating.complete, by = "plate.order") %>%
  mutate(Volume.ul. = 20,
            Experimental.Condition = "Experimental",
            Tissue.Type = "Kernel") %>%
  dplyr::select(Project.Name, User.Email, Plate.Name, Well, NCBITaxId,
                Sample.Name, Concentration.ng.ul., Volume.ul., Experimental.Condition, Tissue.Type,
                Within.plate.num, everything()) %>%
  rbind(controls.master) %>%
  arrange(as.numeric(rowname)) %>%
  drop_na(Sample.Name) %>%
  mutate(Concentration.ng.ul. = floor(as.numeric(Concentration.ng.ul.)))

write.csv(plate.master, "./output/plate_master_v2.csv", row.names = F)

plate.toprint.extractionorder <- plate.master %>% mutate(plate.col.num = parse_number(as.character(Well)),
                                         plate.row.letter = substr(Well, 1,1)) %>%
  mutate(plate.row.number = match(plate.row.letter, toupper(letters))) %>%
  dplyr::select(Sample.Name, done_box_number, Tube_number, Plate.Name, Well, dilution.H2O, dilution.sample) %>%
  arrange(done_box_number, Tube_number)

plate.toprint.plateorder <- plate.toprint.extractionorder %>%
  arrange(Plate.Name, Well)

write.csv(plate.toprint.extractionorder, "./output/master_extraction_order.csv", row.names = F)
write.csv(plate.toprint.plateorder, "./output/master_plate_order.csv", row.names = F)


#### Visualize ####
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}
plate.master %>% mutate(plate.col.num = parse_number(as.character(Well)),
                                 plate.row.letter = substr(Well, 1,1)) %>%
  mutate(plate.row.number = match(plate.row.letter, toupper(letters))) %>%
  ggplot(aes(x = plate.col.num, y = plate.row.number, fill = Experimental.Condition)) +
  geom_tile(color = "black") + facet_grid(~Plate.Name) + scale_y_reverse(breaks = 1:8,labels = toupper(letters[1:8])) +
  scale_x_continuous(breaks = integer_breaks()) + geom_text(aes(label = Sample.Name), size = 1, angle = 45) +
  theme(legend.position="bottom") + labs(title = "Sweet corn RNA plate layout", x = "Plate Column", y = "Plate Row")


# Plate 1
# blank = C04
plate.master %>% filter(Plate.Name == "SC_RNA_01") %>%
  dplyr::select(colnames(plate.template), -Within.plate.num) %>%
  write.table(x = ., file = "./output/SC_RNA_01.txt", sep = "\t", row.names = F)

# Plate 2
# blank = G08
plate.master %>% filter(Plate.Name == "SC_RNA_02") %>%
  dplyr::select(colnames(plate.template), -Within.plate.num) %>%
  write.table(x = ., file = "./output/SC_RNA_02.txt", sep = "\t", row.names = F)

# Plate 3
# blank = D08
plate.master %>% filter(Plate.Name == "SC_RNA_03") %>%
  dplyr::select(colnames(plate.template), -Within.plate.num) %>%
  write.table(x = ., file = "./output/SC_RNA_03.txt", sep = "\t", row.names = F)

# Plate 4
# blank = D07
plate.master %>% filter(Plate.Name == "SC_RNA_04") %>%
  dplyr::select(colnames(plate.template), -Within.plate.num) %>%
  write.table(x = ., file = "./output/SC_RNA_04.txt", sep = "\t", row.names = F)

# Plate 5
# blank = D03
plate.master %>% filter(Plate.Name == "SC_RNA_05") %>%
  dplyr::select(colnames(plate.template), -Within.plate.num) %>%
  write.table(x = ., file = "./output/SC_RNA_05_v2.txt", sep = "\t", row.names = F)

plate.master %>% dplyr::select(colnames(plate.template), -Within.plate.num) %>%
  write.csv(x = ., file = "./output/SC_RNA_plate_master_for_submission.csv", row.names = F)


#### Match to genotype to select positive control ####
field_layout <- read.csv("../GoreLabBase/2019/formatted_for_upload/trials/rnaseq_trial_2019_upload.csv")
field_layout %>%
  dplyr::select(plot_name, accession_name) %>%
  rename(Sample.Name = plot_name) %>%
  full_join(plate.master) %>%
  filter(accession_name == "fill_c2" | accession_name == "CHECK2") %>%
  filter(initial.conc > 1000) %>%
  dplyr::select(Sample.Name, accession_name, Plate.Name, Well, initial.conc:dilution.sample, Extraction_day_1_date) %>%
  arrange(Plate.Name, Well)



#### plate 1 ####
plate.1.fig <- plate.master %>% filter(Plate.Name == "SC_RNA_01") %>%
  mutate(plate.col.num = parse_number(as.character(Well)),
                        plate.row.letter = substr(Well, 1,1)) %>%
  mutate(plate.row.number = match(plate.row.letter, toupper(letters))) %>%
  ggplot(aes(x = plate.col.num, y = plate.row.number))+#, fill = Experimental.Condition)) +
  geom_tile(color = "black", alpha = 0) +
  #facet_grid(~Plate.Name) +
  scale_y_reverse(breaks = 1:8,labels = toupper(letters[1:8])) +
  scale_x_continuous(breaks = integer_breaks()) +
  geom_text(aes(label = paste0(Sample.Name, "\n", dilution.sample)), size = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Sweet corn RNA plate 1 layout", x = "Plate Column", y = "Plate Row")

ggsave(plate.1.fig, filename = "./output/plate.1.layout.png",device = "png",width = 10, height = 7.5,units = "in",
       bg =  "transparent")

#### plate 2 ####
plate.2.fig <- plate.master %>% filter(Plate.Name == "SC_RNA_02") %>%
  mutate(plate.col.num = parse_number(as.character(Well)),
         plate.row.letter = substr(Well, 1,1)) %>%
  mutate(plate.row.number = match(plate.row.letter, toupper(letters))) %>%
  ggplot(aes(x = plate.col.num, y = plate.row.number))+#, fill = Experimental.Condition)) +
  geom_tile(color = "black", alpha = 0) +
  #facet_grid(~Plate.Name) +
  scale_y_reverse(breaks = 1:8,labels = toupper(letters[1:8])) +
  scale_x_continuous(breaks = integer_breaks()) +
  geom_text(aes(label = paste0(Sample.Name, "\n", dilution.sample)), size = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Sweet corn RNA plate 2 layout", x = "Plate Column", y = "Plate Row")

ggsave(plate.2.fig, filename = "./output/plate.2.layout.png",device = "png",width = 10, height = 7.5,units = "in",
       bg =  "transparent")


#### plate 3 ####
plate.3.fig <- plate.master %>% filter(Plate.Name == "SC_RNA_03") %>%
  mutate(plate.col.num = parse_number(as.character(Well)),
         plate.row.letter = substr(Well, 1,1)) %>%
  mutate(plate.row.number = match(plate.row.letter, toupper(letters))) %>%
  ggplot(aes(x = plate.col.num, y = plate.row.number))+#, fill = Experimental.Condition)) +
  geom_tile(color = "black", alpha = 0) +
  #facet_grid(~Plate.Name) +
  scale_y_reverse(breaks = 1:8,labels = toupper(letters[1:8])) +
  scale_x_continuous(breaks = integer_breaks()) +
  geom_text(aes(label = paste0(Sample.Name, "\n", dilution.sample)), size = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Sweet corn RNA plate 3 layout", x = "Plate Column", y = "Plate Row")

ggsave(plate.3.fig, filename = "./output/plate.3.layout.png",device = "png",width = 10, height = 7.5,units = "in",
       bg =  "transparent")


#### plate 4 ####
plate.4.fig <- plate.master %>% filter(Plate.Name == "SC_RNA_04") %>%
  mutate(plate.col.num = parse_number(as.character(Well)),
         plate.row.letter = substr(Well, 1,1)) %>%
  mutate(plate.row.number = match(plate.row.letter, toupper(letters))) %>%
  ggplot(aes(x = plate.col.num, y = plate.row.number))+#, fill = Experimental.Condition)) +
  geom_tile(color = "black", alpha = 0) +
  #facet_grid(~Plate.Name) +
  scale_y_reverse(breaks = 1:8,labels = toupper(letters[1:8])) +
  scale_x_continuous(breaks = integer_breaks()) +
  geom_text(aes(label = paste0(Sample.Name, "\n", dilution.sample)), size = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Sweet corn RNA plate 4 layout", x = "Plate Column", y = "Plate Row")

ggsave(plate.4.fig, filename = "./output/plate.4.layout.png",device = "png",width = 10, height = 7.5,units = "in",
       bg =  "transparent")



#### plate 5 ####
plate.5.fig <- plate.master %>% filter(Plate.Name == "SC_RNA_05") %>%
  mutate(plate.col.num = parse_number(as.character(Well)),
         plate.row.letter = substr(Well, 1,1)) %>%
  mutate(plate.row.number = match(plate.row.letter, toupper(letters))) %>%
  ggplot(aes(x = plate.col.num, y = plate.row.number))+#, fill = Experimental.Condition)) +
  geom_tile(color = "black", alpha = 0) +
  #facet_grid(~Plate.Name) +
  scale_y_reverse(breaks = 1:8,labels = toupper(letters[1:8])) +
  scale_x_continuous(breaks = integer_breaks()) +
  geom_text(aes(label = paste0(Sample.Name, "\n", dilution.sample)), size = 3) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Sweet corn RNA plate 5 layout", x = "Plate Column", y = "Plate Row")

ggsave(plate.5.fig, filename = "./output/plate.5.layout.png",device = "png",width = 10, height = 7.5,units = "in",
       bg =  "transparent")

positive.control.conc = (1266.154+1237.101)/2
pos.control.H2O = get_dilution(positive.control.conc)
pos.control.sample = 20 - pos.control.H2O


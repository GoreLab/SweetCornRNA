# calculate BLUE via asreml

# get argument from shell
args <- commandArgs(trailingOnly = T)
s <- as.numeric(args[1])
ref <- args[2] # "B73" or "PH207"
rand.seed <- 2018

# source
library(asreml)
library(data.table)
library(doParallel)

# file I
filename.data <- paste0("RESULT/1.3-OutlierRemoval_and_Imputation/ExpressionData_For_BLUE_", ref, ".csv")
file.in.key <- "RAWDATA/master_key.csv"
file.in.expinfo <- "RAWDATA/Metadata/rnaseq_trial_2019_upload.csv"

# file O
foldername.save <- paste0("RESULT_NoHarvDate/3.1-BLUE/ResEach_", ref)
filename.save.1 <- paste0("RESULT_NoHarvDate/3.1-BLUE/ResEach_", ref, "/exp_BLUE_all_", s, ".Rdata")
filename.save.2 <- paste0("RESULT_NoHarvDate/3.1-BLUE/ResEach_", ref, "/exp_BLUE_pred_", s, ".Rdata")
filename.save.3 <- paste0("RESULT_NoHarvDate/3.1-BLUE/ResEach_", ref, "/exp_BLUE_lastm_", s, ".Rdata")
filename.save.4 <- paste0("RESULT_NoHarvDate/3.1-BLUE/ResEach_", ref, "/exp_BLUE_seed_", s, ".txt")

# make folder to save result
dir.create(foldername.save, recursive = TRUE)

# load data
exp.data <- fread(filename.data)
sample.name <- exp.data[, 1]
rlog.mat <- as.matrix(exp.data[, -1])
rownames(rlog.mat) <- unlist(sample.name)

# load master data
key.df <- read.csv(file = file.in.key)
setdiff(rownames(rlog.mat), key.df$rlog_ids) # OK

# load experiment info data
expinfo.df <- read.csv(file = file.in.expinfo)
setdiff(key.df$sample_name_metadata, expinfo.df$plot_name) # OK -- "Pos" is positive control

# add rlog value of the s-th gene to the key file
df.model <- data.frame("rlog.id" = rownames(rlog.mat),
											 "Expr.rlog" = rlog.mat[, s],
											 row.names = NULL)
df.model <- merge.data.frame(x = df.model, y = key.df, by.x = "rlog.id", by.y = "rlog_ids")
df.model$accession_name <- as.factor(as.character(df.model$accession_name))

# merge metadata
length(expinfo.df$plot_name) == length(unique(expinfo.df$plot_name)) # uniqueness: OK
length(df.model$sample_name_metadata) == length(unique(df.model$sample_name_metadata)) # uniqueness: OK
setdiff(df.model$sample_name_metadata, expinfo.df$plot_name) # check inclusion: OK
df.model <- merge.data.frame(x = df.model, y = expinfo.df,
														 by.x = "sample_name_metadata", by.y = "plot_name",
														 all.x = TRUE, all.y = FALSE)

# rename columns for simplicity
colnames(df.model)[colnames(df.model) == "accession_name.y"] <- "accession_name_meta"
colnames(df.model)[colnames(df.model) == "accession_name.x"] <- "accession_name"

# Harvest.Date
days.harvest <- as.integer(as.Date(df.model$tissue_harvest_date) - min(as.Date(df.model$tissue_harvest_date)))
days.harvest.sc <- scale(days.harvest)

# days.rna.extract should be two-level factor
rna.extract <- rep("spring", nrow(df.model))
rna.extract["2020-05-01" < as.Date(df.model$nucleic_acid_extraction_date)] <- "summer"
rna.extract <- as.factor(rna.extract)

# final data
df.fit <- data.frame("y" = df.model$Expr.rlog,
										 "Genotype" = df.model$accession_name,
										 "Range" = as.factor(df.model$range_number),
										 "Block" = as.factor(df.model$block_number),
										 "Column" = as.factor(df.model$col_number),
										 "Row" = as.factor(df.model$row_number),
										 "Plate" = as.factor(df.model$plate_number),
										 "Harvest.Date" = days.harvest.sc,
										 "RnaExtract" = as.factor(rna.extract))

# save the data frame to build BLUE/BLUP model
if ( s == 1 ) {
	df.save <- df.fit[, -1]
	write.csv(df.save, file = paste0(foldername.save, "/DataFrame_BLUE.csv"), row.names = F)
}

# run asreml
set.seed(rand.seed)
asr.mod <- asreml(fixed = y ~ -1 + Genotype,
									random = ~ Range + Range:Block + Range:Block:Column + Row + Plate + RnaExtract,
									maxiter = 100,
									data = df.fit,
									trace = F)
last.m <- asr.mod$last.message
pred <- predict(asr.mod, classify = "Genotype", data = df.fit)
pred.val <- pred$predictions$pvals$predicted.value
names(pred.val) <- pred$predictions$pvals$Genotype

# make object for the output
coef.all <- coef(asr.mod)
coef.fixed <- setNames(coef.all$fixed, rownames(coef.all$fixed))
coef.random <- setNames(coef.all$random, rownames(coef.all$random))
gamma.rand <- asr.mod$gammas
asr.summary.all <- c(coef.fixed, coef.random, gamma.rand)

# write
saveRDS(asr.summary.all, file = filename.save.1)
saveRDS(pred.val, file = filename.save.2)
saveRDS(last.m, file = filename.save.3)
write(rand.seed, file = filename.save.4)

# print
if ( asr.mod$converge == FALSE ) { print(s); print("Not converged") }

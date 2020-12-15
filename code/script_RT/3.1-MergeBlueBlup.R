# merge BLUP results

# source
library(e1071)
library(ggplot2)
library(reshape2)

# file I/O
dir.in.rlog <- "RAWDATA/Seetcorn_TagSeq"
file.in.rlog <- "htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.txt"
dir.in.key <- "RAWDATA"
file.in.key <- "master_key.csv"
dir.blue.in.lmer <- "RESULT/2.1-BLUE_lmer"
dir.blup.in.lmer <- "RESULT/2.2-BLUP_lmer"
dir.save <- "RESULT/3.1-MergeBlueBlup_lmer"

# objects
names.var.BLUE <- c("var.range", "var.block", "var.column", "var.row", 
                    "var.plate", "var.extract", "var.resid")
names.var.BLUP <- c("var.geno", "var.range", "var.block", "var.column", "var.row", 
                     "var.plate", "var.extract", "var.resid")
names.check <- c("CHECK1", "CHECK2", "CHECK3", "CHECK4")


# make dir to save result
dir.create(dir.save, recursive = TRUE)


# --------------------------------------------------------------------------- #
# ----- 1. BLUE
# --------------------------------------------------------------------------- #
# load rlog data
rlog.df <- read.delim(file = paste0(dir.in.rlog, "/", file.in.rlog)) # a few min

# rlog data matrix
gene.info.df <- rlog.df[, 1:5]
rlog.mat <- as.matrix(rlog.df[, 6:ncol(rlog.df)])
rownames(rlog.mat) <- gene.info.df$gene_id

# load key file
data.key <- read.csv(paste0(dir.in.key, "/", file.in.key))

# get numbers for which we could not get BLUE
filenames.all <- dir(dir.blue.in.lmer)
error.files <- filenames.all[grep("num_error", filenames.all)]
num.error <- c()
for ( i in 1:length(error.files) ) {
  file.i <- error.files[i]
  res.try <- try(
    df.error <- read.csv(file = paste0(dir.blue.in.lmer, "/", file.i))
  )
  if ( class(res.try) != "try-error" ) {
    num.error <- c(num.error, df.error$num.error)
  }
}

# get BLUE
num.get <- setdiff(1:nrow(rlog.mat), num.error)
df.BLUE.all <- NULL
for (i in num.get) {
  num.i <- formatC(i, digits = 5, flag = "0")
  file.i <- paste0(dir.blue.in.lmer, "/BLUE_", num.i, ".csv")
  dat.i <- read.csv(file.i)
  if ( i == num.get[1] ) {
    df.BLUE.all <- dat.i
  } else {
    df.BLUE.all <- cbind.data.frame(df.BLUE.all, dat.i$Est.effect)
  }
}
colnames(df.BLUE.all)[2:ncol(df.BLUE.all)] <- rownames(rlog.mat)[num.get]

# make a data frame of variance components
df.var <- df.BLUE.all[df.BLUE.all$Model.Term %in% names.var.BLUE, ]
write.csv(df.var, file = paste0(dir.save, "/EstVarComp_BLUE_lmer.csv"), row.names = F)

# make a data fram of BLUEs
geno.common <- setdiff(intersect(df.BLUE.all$Model.Term, 
                                 data.key$accession_name), 
                       names.check)
df.BLUE <- df.BLUE.all[df.BLUE.all$Model.Term %in% geno.common, ]
write.csv(df.BLUE, file = paste0(dir.save, "/BLUE_lmer.csv"), row.names = F)

# make a data fram of other effects
name.others <- setdiff(df.BLUE.all$Model.Term, c(df.var$Model.Term, df.BLUE$Model.Term))
df.others <- df.BLUE.all[df.BLUE.all$Model.Term %in% name.others, ]
write.csv(df.others, file = paste0(dir.save, "/OtherEffects_BLUE_lmer.csv"), row.names = F)




# --------------------------------------------------------------------------- #
# ----- 1. BLUP
# --------------------------------------------------------------------------- #
# get numbers for which we could not get BLUP
filenames.all <- dir(dir.blup.in.lmer)
error.files <- filenames.all[grep("num_error", filenames.all)]
num.error <- c()
for ( i in 1:length(error.files) ) {
  file.i <- error.files[i]
  res.try <- try(
    df.error <- read.csv(file = paste0(dir.blup.in.lmer, "/", file.i))
  )
  if ( class(res.try) != "try-error" ) {
    num.error <- c(num.error, df.error$num.error)
  }
}

# get BLUP
num.get <- setdiff(1:nrow(rlog.mat), num.error)
df.BLUP.all <- NULL
for (i in num.get) {
  num.i <- formatC(i, digits = 5, flag = "0")
  file.i <- paste0(dir.blup.in.lmer, "/BLUP_", num.i, ".csv")
  dat.i <- read.csv(file.i)
  if ( i == num.get[1] ) {
    df.BLUP.all <- dat.i
  } else {
    df.BLUP.all <- cbind.data.frame(df.BLUP.all, dat.i$Est.effect)
  }
}
colnames(df.BLUP.all)[2:ncol(df.BLUP.all)] <- rownames(rlog.mat)[num.get]

# make a data frame of variance components
df.var <- df.BLUP.all[df.BLUP.all$Model.Term %in% names.var.BLUP, ]
write.csv(df.var, file = paste0(dir.save, "/EstVarComp_BLUP_lmer.csv"), row.names = F)

# make a data frame of BLUPs
geno.common <- setdiff(intersect(df.BLUP.all$Model.Term, 
                                 data.key$accession_name), 
                       names.check)
df.BLUP <- df.BLUP.all[df.BLUP.all$Model.Term %in% geno.common, ]
write.csv(df.BLUP, file = paste0(dir.save, "/BLUP_lmer.csv"), row.names = F)

# make a data frame of other effects
name.others <- setdiff(df.BLUP.all$Model.Term, c(df.var$Model.Term, df.BLUP$Model.Term))
df.others <- df.BLUP.all[df.BLUP.all$Model.Term %in% name.others, ]
write.csv(df.others, file = paste0(dir.save, "/OtherEffects_BLUP_lmer.csv"), row.names = F)






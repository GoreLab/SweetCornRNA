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
dir.in.lmer <- "RESULT/2.2-BLUP_lmer"


dir.save <- "RESULT/3.2-MergeBLUP"

# objects

# make dir to save result
dir.create(dir.save, recursive = TRUE)

# --------------------------------------------------------------------------- #
# ----- 1. Load data & make a few objects
# --------------------------------------------------------------------------- #
# load rlog data
rlog.df <- read.delim(file = paste0(dir.in.rlog, "/", file.in.rlog)) # a few min

# rlog data matrix
gene.info.df <- rlog.df[, 1:5]
rlog.mat <- as.matrix(rlog.df[, 6:ncol(rlog.df)])
rownames(rlog.mat) <- gene.info.df$gene_id

# get numbers for which we could not get BLUP
filenames.all <- dir("RESULT/2.2-BLUP_lmer")
error.files <- filenames.all[grep("num_error", filenames.all)]
num.error <- c()
for ( i in 1:length(error.files) ) {
  file.i <- error.files[i]
  df.error <- read.csv(file = paste0(dir.in.lmer, "/", file.i))
  num.error <- c(num.error, df.error$num.error)
}

# look at the ones we could not get BLUP
n.zero <- apply(rlog.mat == 0, 1, sum, na.rm = T)
df.fig <- data.frame("gene" = rownames(rlog.mat),
                     "n.zero" = n.zero,
                     "pr.zero" = n.zero / ncol(rlog.mat),
                     "error.BLUP" = "No Error")
df.fig$error.BLUP[num.error] <- "Error"
p <- ggplot(df.fig, aes(x = pr.zero))
p <- p + geom_histogram()
p <- p + facet_wrap(~ error.BLUP, scales = "free_y", ncol = 1)
p <- p + xlab("% of zero rlog2 value")
p <- p + ggtitle("When there are many zero rlog2, we may get error")
p

# get BLUP
num.get <- setdiff(1:nrow(rlog.mat), num.error)
df.BLUP.all <- NULL
for (i in num.get) {
  num.i <- formatC(i, digits = 5, flag = "0")
  file.i <- paste0(dir.in.lmer, "/BLUP_", num.i, ".csv")
  dat.i <- read.csv(file.i)
  if ( i == num.get[1] ) {
    df.BLUP.all <- dat.i
  } else {
    df.BLUP.all <- cbind.data.frame(df.BLUP.all, dat.i$Est.effect)
  }
}
colnames(df.BLUP.all)[2:ncol(df.BLUP.all)] <- rownames(rlog.mat)[num.get]

# Plate effects
name.plate <- c("plate1", "plate2", "plate3", "plate4", "plate5")
sd.vec <- apply(rlog.mat, 1, sd, na.rm = T)
df.plate.effect <- df.BLUP.all[df.BLUP.all$Model.Term %in% name.plate, ]
for ( i in 2:ncol(df.plate.effect) ) {
  gene.i <- colnames(df.plate.effect)[i]
  df.plate.effect[, i] <- df.plate.effect[, i] / sd.vec[gene.i]
}
df.fig <- melt(df.plate.effect, )
p <- ggplot(df.fig, aes(x = Model.Term, y = value))
p <- p + geom_boxplot()
p <- p + ylab("Scaled plate effect = est.effect / sd(Y)")
p <- p + ggtitle("Estimated plate effect")
p

# heritability
name.var <- c("var.geno", "var.range", "var.block", "var.column", "var.row",
              "var.plate", "var.extract", "var.resid")




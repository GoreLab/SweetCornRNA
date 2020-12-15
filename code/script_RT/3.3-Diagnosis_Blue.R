# merge blue results

# source
library(e1071)
library(ggplot2)
library(reshape2)
library(data.table)

# file I/O
dir.in.rlog <- "RAWDATA/Seetcorn_TagSeq"
file.in.rlog <- "htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.txt"
dir.blue.in.lmer <- "RESULT/3.1-MergeBlueBlup_lmer"
file.in.blue <- "BLUE_lmer.csv"
file.in.var.blue <- "EstVarComp_BLUE_lmer.csv"
file.in.others.blue <- "OtherEffects_BLUE_lmer.csv"
dir.save <- "RESULT/3.3-Diagnosis_Blue"

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

# get blue results
blue <- fread(file = paste0(dir.blue.in.lmer, "/", file.in.blue), data.table = F)
blue.mat <- as.matrix(blue[, -1])
rownames(blue.mat) <- blue$Model.Term

# get var comp
df.var <- fread(file = paste0(dir.blue.in.lmer, "/", file.in.var.blue), data.table = F)
var.mat <- as.matrix(df.var[, -1])
rownames(var.mat) <- df.var$Model.Term

# get other effects
df.others <- fread(file = paste0(dir.blue.in.lmer, "/", file.in.others.blue), data.table = F)

# --------------------------------------------------------------------------- #
# ----- 2. Make a summary
# --------------------------------------------------------------------------- #
# list of genes for which we could not get blue
genes.out <- setdiff(rownames(rlog.mat), colnames(blue.mat))
write(genes.out, file = paste0(dir.save, "/Genes_Noblue.txt")) # 3 genes

# variance comp
m <- match(colnames(blue.mat), rownames(rlog.mat))
rlog.mat.for.blue <- t(rlog.mat[m, ])
sd.obs <- apply(rlog.mat.for.blue, 2, sd, na.rm = T)
var.mat.scaled <- t(var.mat) / sd.obs
colnames(var.mat.scaled) <- c("Range", "Block", "Column", "Row",
                              "Plate", "RNA.Extract", "Residual")
df.fig <- reshape2::melt(var.mat.scaled)
colnames(df.fig) <- c("Gene", "Var.Comp", "Est.Var")

# make a figure
p <- ggplot(df.fig, aes(x = Var.Comp, y = Est.Var))
p <- p + geom_boxplot()
p <- p + xlab("Model Term") + ylab("Estimated Variance (scaled by var(obs))")
p <- p + ggtitle("Estimated variance components in BLUE model")
ggsave(filename = paste0(dir.save, "/EstVarComp_v1.png"), p, width = 9, height = 7)

# make a figure
var.mat.scaled <- t(var.mat) / apply(var.mat, 2, sum)
colnames(var.mat.scaled) <- c("Range", "Block", "Column", "Row",
                              "Plate", "RNA.Extract", "Residual")
df.fig <- reshape2::melt(var.mat.scaled)
colnames(df.fig) <- c("Gene", "Var.Comp", "Est.Var")
p <- ggplot(df.fig, aes(x = Var.Comp, y = Est.Var))
p <- p + geom_boxplot()
p <- p + xlab("Model Term") + ylab("Estimated Variance (scaled by sum(var.comp))")
p <- p + ggtitle("Estimated variance components in BLUE model (v2)")
ggsave(filename = paste0(dir.save, "/EstVarComp_v2.png"), p, width = 9, height = 7)

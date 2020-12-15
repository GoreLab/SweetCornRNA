# merge BLUP results

# source
library(e1071)
library(ggplot2)
library(reshape2)
library(data.table)

# file I/O
dir.in.rlog <- "RAWDATA/Seetcorn_TagSeq"
file.in.rlog <- "htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.txt"
dir.blup.in.lmer <- "RESULT/3.1-MergeBlueBlup_lmer"
file.in.blup <- "BLUP_lmer.csv"
file.in.var.blup <- "EstVarComp_BLUP_lmer.csv"
file.in.others.blup <- "OtherEffects_BLUP_lmer.csv"
dir.save <- "RESULT/3.2-Diagnosis_Blup"

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

# get BLUP results
BLUP <- fread(file = paste0(dir.blup.in.lmer, "/", file.in.blup), data.table = F)
BLUP.mat <- as.matrix(BLUP[, -1])
rownames(BLUP.mat) <- BLUP$Model.Term

# get var comp
df.var <- fread(file = paste0(dir.blup.in.lmer, "/", file.in.var.blup), data.table = F)
var.mat <- as.matrix(df.var[, -1])
rownames(var.mat) <- df.var$Model.Term

# get other effects
df.others <- fread(file = paste0(dir.blup.in.lmer, "/", file.in.others.blup), data.table = F)

# --------------------------------------------------------------------------- #
# ----- 2. Make a summary
# --------------------------------------------------------------------------- #
# list of genes for which we could not get BLUP
genes.out <- setdiff(rownames(rlog.mat), colnames(BLUP.mat))
write(genes.out, file = paste0(dir.save, "/Genes_NoBlup.txt")) # 112 genes

# variance comp
m <- match(colnames(BLUP.mat), rownames(rlog.mat))
rlog.mat.for.BLUP <- t(rlog.mat[m, ])
sd.obs <- apply(rlog.mat.for.BLUP, 2, sd, na.rm = T)
var.mat.scaled <- t(var.mat) / sd.obs
colnames(var.mat.scaled) <- c("Genotype", "Range", "Block", "Column", "Row",
                              "Plate", "RNA.Extract", "Residual")
df.fig <- reshape2::melt(var.mat.scaled)
colnames(df.fig) <- c("Gene", "Var.Comp", "Est.Var")

# make a figure
p <- ggplot(df.fig, aes(x = Var.Comp, y = Est.Var))
p <- p + geom_boxplot()
p <- p + xlab("Model Term") + ylab("Estimated Variance (scaled by var(obs))")
p <- p + ggtitle("Estimated variance components in BLUP model")
ggsave(filename = paste0(dir.save, "/EstVarComp_v1.png"), p, width = 9, height = 7)

# make a figure
var.mat.scaled <- t(var.mat) / apply(var.mat, 2, sum)
colnames(var.mat.scaled) <- c("Genotype", "Range", "Block", "Column", "Row",
                              "Plate", "RNA.Extract", "Residual")
df.fig <- reshape2::melt(var.mat.scaled)
colnames(df.fig) <- c("Gene", "Var.Comp", "Est.Var")
p <- ggplot(df.fig, aes(x = Var.Comp, y = Est.Var))
p <- p + geom_boxplot()
p <- p + xlab("Model Term") + ylab("Estimated Variance (scaled by sum(var.comp))")
p <- p + ggtitle("Estimated variance components in BLUP model (v2)")
ggsave(filename = paste0(dir.save, "/EstVarComp_v2.png"), p, width = 9, height = 7)

# heritability
h2 <- var.mat["var.geno", ] / apply(var.mat, 2, sum)
df.fig <- data.frame("Heritability" = h2)
p <- ggplot(df.fig, aes(x = Heritability))
p <- p + geom_histogram(color = "white")
p <- p + ggtitle("Model-based heritability of transcriptome")
ggsave(filename = paste0(dir.save, "/Heritability.png"), p, width = 9, height = 7)

# heritability
sum(is.na(h2)) # residual variance = 0 in this case (i.e., h2 = +inf)
sum(h2 == 0, na.rm = T) # 316 genes, h2 = 0
sum(h2 < 0.01, na.rm = T) # 900 genes, h2 < 1%








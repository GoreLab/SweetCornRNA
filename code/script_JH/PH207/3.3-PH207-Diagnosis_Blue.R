# merge blue results

# source
library(e1071)
library(ggplot2)
library(reshape2)
library(data.table)

# file I/O
dir.in.rlog <- "data"
file.in.rlog <- "htseq_count_matrix_sweetcorn_PH207_RLOG_all_info_v1.txt"
dir.blue.in.lmer <- "output/PH207_BLUEs"
file.in.blue <- "PH207_BLUE_lmer.csv"
file.in.var.blue <- "PH207_EstVarComp_BLUE_lmer.csv"
file.in.others.blue <- "PH207_OtherEffects_BLUE_lmer.csv"
dir.save <- "output/PH207_BLUEs"

# objects

# make dir to save result
dir.create(dir.save, recursive = TRUE)

# --------------------------------------------------------------------------- #
# ----- 1. Load data & make a few objects
# --------------------------------------------------------------------------- #
# load rlog data
rlog.df <- read.delim(file = paste0(dir.in.rlog, "/", file.in.rlog)) # a few min

# rlog data matrix
gene.info.df <- rlog.df[, 1:6]
rlog.mat <- as.matrix(rlog.df[, 7:ncol(rlog.df)])
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
write(genes.out, file = paste0(dir.save, "/PH207_Genes_Noblue.txt")) # 1 gene: Zm00008a003145

# variance comp
m <- match(colnames(blue.mat), rownames(rlog.mat))
rlog.mat.for.blue <- t(rlog.mat[m, ])
sd.obs <- apply(rlog.mat.for.blue, 2, sd, na.rm = T)

var.mat.scaled.var.obs <- t(var.mat) / sd.obs
colnames(var.mat.scaled.var.obs) <- c("Range", "Block", "Column", "Row",
                              "Plate", "RNA.Extract", "Residual")
df.fig.1 <- reshape2::melt(var.mat.scaled.var.obs)
colnames(df.fig.1) <- c("Gene", "Var.Comp", "Est.Var")

# make a figure
p1 <- ggplot(df.fig.1, aes(x = Var.Comp, y = Est.Var)) +
  geom_boxplot() +
  labs(x = "Model Term",
       y = "Estimated Variance (scaled by var(obs))",
       title = "Estimated variance components in BLUE model (v2)")
p1
ggsave(filename = paste0(dir.save, "/PH207_EstVarComp_v1.png"), p1, width = 9, height = 7)

# make a figure
var.mat.scaled.sum.var.comp <- t(var.mat) / apply(var.mat, 2, sum)
colnames(var.mat.scaled.sum.var.comp) <- c("Range", "Block", "Column", "Row",
                              "Plate", "RNA.Extract", "Residual")
df.fig.2 <- reshape2::melt(var.mat.scaled.sum.var.comp)
colnames(df.fig.2) <- c("Gene", "Var.Comp", "Est.Var")
p2 <- ggplot(df.fig.2, aes(x = Var.Comp, y = Est.Var)) +
  geom_boxplot() +
  labs(x = "Model Term",
       y = "Estimated Variance (scaled by sum(var.comp))",
       title = "Estimated variance components in BLUE model (v2)")
p2
ggsave(filename = paste0(dir.save, "/PH207_EstVarComp_v2.png"), p2, width = 9, height = 7)


# Check on genes identified as potential outliers
# Zm00008a023627 for SC_RNA_05_G06_19A0208_MDM
# Zm00008a022414 for SC_RNA_02_D10_19A0088_Wh10140R
# Zm00008a027192 for SC_RNA_02_C07_19A0094_Wuh07466i

which.max(blue.mat[,"Zm00008a023627"])
hist(blue.mat[,"Zm00008a023627"], breaks = 20,
     xlab = "BLUE value",
     main = "Histogram of Zm00008a023627 BLUEs for 432 samples",
     col = "grey")

which.max(blue.mat[,"Zm00008a022414"])
hist(blue.mat[,"Zm00008a022414"], breaks = 20,
     xlab = "BLUE value",
     main = "Histogram of Zm00008a022414 BLUEs for 432 samples",
     col = "grey")

which.max(blue.mat[,"Zm00008a027192"])
hist(blue.mat[,"Zm00008a027192"], breaks = 20,
     xlab = "BLUE value",
     main = "Histogram of Zm00008a027192 BLUEs for 432 samples",
     col = "grey")






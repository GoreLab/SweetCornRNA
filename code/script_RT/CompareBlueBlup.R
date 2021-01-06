# Example to interpret the BLUE/BLUP results

# library
library(ggplot2)

# load results
dat.BLUE <- read.csv("RESULT/2.1-BLUE_lmer/BLUE_005332.csv")
dat.BLUP <- read.csv("RESULT/2.2-BLUP_lmer/BLUP_005332.csv")

# e.g., BLUE vs BLUP
name.exp.gen <- dat.BLUP$Model.Term[14:371]

m1 <- match(name.exp.gen, dat.BLUP$Model.Term)
blup <- dat.BLUP$Est.effect[m1] + dat.BLUP$Est.effect[dat.BLUP$Model.Term == "Intercept"]

m2 <- match(name.exp.gen, dat.BLUE$Model.Term)
blue <- dat.BLUE$Est.effect[m2]

lim <- range(c(blue, blup))
plot(x = blue, y = blup, 
     xlab = "BLUE", ylab = "BLUP",
     main = "BLUE vs BLUP",
     xlim = lim, ylim = lim, pch = 20)
abline(0, 1, lty = 2)

# Estimated variance components in BLUP
name.varcomp <- c("var.geno", "var.range", "var.block", "var.column", "var.row",
                  "var.plate", "var.extract", "var.resid")
rename.varcomp <- c("Genotype", "Range", "Block", "Col", "Row", "Plate", "Extract.Date", "Resid")
m <- match(name.varcomp, dat.BLUP$Model.Term)
df.fig <- dat.BLUP[m, ]
df.fig$Model.Term <- factor(rename.varcomp, levels = rename.varcomp)

p <- ggplot(df.fig, aes(x = Model.Term, y = Est.effect))
p <- p + geom_bar(stat = "identity")
p <- p + ylab("Estimated Variance")
p <- p + ggtitle("Estimated Variacne in BLUP model")
p

# e.g., visualize plate effect in BLUE model
name.plate <- c("plate1", "plate2", "plate3", "plate4", "plate5")
m <- match(name.plate, dat.BLUE$Model.Term)
df.fig <- dat.BLUE[m, ]

p <- ggplot(df.fig, aes(x = Model.Term, y = Est.effect))
p <- p + geom_bar(stat = "identity")
p <- p + ylab("Estimated Effect")
p <- p + ggtitle("Plate effect in BLUE model")
p



# ---------------------------------------------------------------------------- #
# get raw data
dir.in.rlog <- "RAWDATA/Seetcorn_TagSeq"
file.in.rlog <- "htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.txt"
rlog.df <- read.delim(file = paste0(dir.in.rlog, "/", file.in.rlog)) # a few min
gene.info.df <- rlog.df[, 1:5]
rlog.mat <- as.matrix(rlog.df[, 6:ncol(rlog.df)])
rownames(rlog.mat) <- gene.info.df$gene_id

# get BLUP results
dir.blup.in.lmer <- "RESULT/3.1-MergeBlueBlup_lmer"
file.in.blup <- "BLUP_lmer.csv"
BLUP <- fread(file = paste0(dir.blup.in.lmer, "/", file.in.blup), data.table = F)
BLUP.mat <- as.matrix(BLUP[, -1])
rownames(BLUP.mat) <- BLUP$Model.Term

# a <- setdiff(rownames(rlog.mat), colnames(BLUP.mat))
# b <- setdiff(rownames(rlog.mat), colnames(BLUE.mat))
# setdiff(a, b)
# setdiff(b, a)

# get BLUE results
dir.blue.in.lmer <- "RESULT/3.1-MergeBlueBlup_lmer"
file.in.blue <- "BLUE_lmer.csv"
BLUE <- fread(file = paste0(dir.blue.in.lmer, "/", file.in.blue), data.table = F)
BLUE.mat <- as.matrix(BLUE[, -1])
rownames(BLUE.mat) <- BLUE$Model.Term

# common gene
common.gene <- intersect(colnames(BLUP.mat), colnames(BLUE.mat))
BLUP.common <- BLUP.mat[, common.gene]
BLUE.common <- BLUE.mat[, common.gene]
genes.not.common <- c(setdiff(colnames(BLUP.mat), common.gene), 
                      setdiff(colnames(BLUE.mat), common.gene))

# cor(BLUE, BLUP)
r.vec <- rep(NA, ncol(BLUP.common))
names(r.vec) <- colnames(BLUP.common)
for ( i in 1:ncol(BLUP.common) ) {
  r.vec[i] <- cor(BLUP.common[, i], BLUE.common[, i])
}
hist(r.vec, main = "Histogram of 20292 genes", xlab = "Cor(BLUE, BLUP)")
# r.vec[is.na(r.vec)] <- 0
# genes.low.cor <- names(r.vec)[r.vec < 0.6]

# zero-proportion
n.zero <- apply(rlog.mat == 0, 1, sum, na.rm = T)
pr.zero <- n.zero / ncol(rlog.mat)
genes.good <- names(pr.zero)[pr.zero < 0.90]
sum(pr.zero < 0.90)

# 
genes.no.variation.BLUE <- colnames(BLUE.common)[apply(BLUE.common, 2, sd) == 0]
genes.no.variation.BLUP <- colnames(BLUP.common)[apply(BLUP.common, 2, sd) == 0]
genes.no.variation <- union(genes.no.variation.BLUP, genes.no.variation.BLUE)
length(genes.no.variation.BLUE)
length(genes.no.variation.BLUP)
length(genes.no.variation)
hist(pr.zero[genes.no.variation],
     main = "Histogram of the 334 'no variation' genes",
     xlab = "Proportion of zero rlog2 values")

hist(pr.zero[genes.no.variation.BLUE],
     main = "Histogram of the 20 'no variation in BLUE' genes",
     xlab = "Proportion of zero rlog2 values")


# cor(BLUE, BLUP)
r.vec <- rep(NA, ncol(BLUP.common))
names(r.vec) <- colnames(BLUP.common)
for ( i in 1:ncol(BLUP.common) ) {
  r.vec[i] <- cor(BLUP.common[, i], BLUE.common[, i])
}
length(genes.good)
hist(r.vec[genes.good],
     main = "Histogram of the 19715 'high quality' genes",
     xlab = "Proportion of zero rlog2 values")
Dwhich.min(r.vec[genes.good])
plot(BLUP.common[, "Zm00001d044527"], BLUE.common[, "Zm00001d044527"], pch = 16)
hist(rlog.mat["Zm00001d044527", ])


# library(MASS)
# sk <- apply(rlog.mat, 1, skewness, na.rm = T)
# hist(sk)
# hist(sk[genes.low.cor])


hist(pr.zero[genes.not.common], 
     main = "Histogram of 115 'calculation-failed' genes",
     xlab = "Proportion of zero rlog2 values")

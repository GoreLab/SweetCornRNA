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
dir.save <- "RESULT/3.4-Filtering_BLUE"

# make dir to save result
dir.create(dir.save, recursive = TRUE)

# --------------------------------------------------------------------------- #
# ----- 1. Load data & make a few objects
# --------------------------------------------------------------------------- #
# load rlog data
rlog.df <- fread(file = paste0(dir.in.rlog, "/", file.in.rlog), data.table = F) # a few min

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
# ----- 3. remove genes with many zero rlog2 values
# --------------------------------------------------------------------------- #
# zero-proportion
n.zero <- apply(rlog.mat == 0, 1, sum, na.rm = T)
pr.zero <- n.zero / ncol(rlog.mat)

# check the filter
tf <- pr.zero < 0.9
genes.retain <- rownames(rlog.mat)[tf]
genes.remove <- rownames(rlog.mat)[!tf]

# check
length(genes.remove) # remove 1160 genes

# filtering
blue.mat.filtered <- blue.mat[, genes.retain]

# --------------------------------------------------------------------------- #
# ----- 4. Save
# --------------------------------------------------------------------------- #
df.save <- data.frame("Model.Term" = rownames(blue.mat.filtered), blue.mat.filtered)
fwrite(df.save, row.names = F, file = paste0(dir.save, "/BLUE_lmer_filtered.csv"))

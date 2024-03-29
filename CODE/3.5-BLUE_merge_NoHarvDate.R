# This script simply merge all BLUE calculation result objects 
# JUST MERGE, NOTHING WILL BE DONE

# package
library(data.table)

# params
args <- commandArgs(trailingOnly = T)
ref <- args[1] # "B73" or "PH207"

# load data
exp.data <- fread(paste0("RESULT/1.3-OutlierRemoval_and_Imputation/ExpressionData_For_BLUE_", ref, ".csv"))
gene.names <- colnames(exp.data)[-1]
N <- length(gene.names)

# ----------- Merge BLUE of accessions ----------- 
# 1st one
filename.load <- paste0("RESULT_NoHarvDate/3.1-BLUE/ResEach_", ref, "/exp_BLUE_pred_1.Rdata")
res.s <- readRDS(filename.load)
P <- length(res.s)
# load all & merge
res.all <- matrix(NA, nr = N, nc = P)
rownames(res.all) <- 1:N
for (s in 1:N) {
   filename.load <- paste0("RESULT_NoHarvDate/3.1-BLUE/ResEach_", ref, "/exp_BLUE_pred_", s, ".Rdata")
   res.s <- readRDS(filename.load)
   res.all[s, ] <- res.s
}
colnames(res.all) <- names(res.s)
df.save <- data.frame("GeneID" = gene.names, res.all)
# write
fwrite(df.save, paste0("RESULT_NoHarvDate/3.1-BLUE/expression_BLUE_", ref, "_raw.csv"))



# ----------- Merge last message (error message/convergence) ----------- 
# 1st one
filename.load <- paste0("RESULT_NoHarvDate/3.1-BLUE/ResEach_", ref, "/exp_BLUE_lastm_1.Rdata")
res.s <- readRDS(filename.load)
# load all & merge
res.all <- data.frame("GeneID" = gene.names, "last.m" = NA)
for (s in 1:N) {
   filename.load <- paste0("RESULT_NoHarvDate/3.1-BLUE/ResEach_", ref, "/exp_BLUE_lastm_", s, ".Rdata")
   res.s <- readRDS(filename.load)
   res.all[s, "last.m"] <- res.s
}
# write
fwrite(res.all, paste0("RESULT_NoHarvDate/3.1-BLUE/expression_BLUE_", ref, "_message.csv"))



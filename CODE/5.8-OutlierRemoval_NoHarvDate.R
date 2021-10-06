# library
library(data.table)
library(asreml)

# params
args <- commandArgs(trailingOnly = T)
ref <- args[1] # "B73" or "PH207"
dat <- args[2] # "toco", "ion", "carot"
rand.seed <- 2018

# mkdir
dir.create("RESULT_NoHarvDate/5.4-OutlierRemoval")

# file I/O
file.in <- paste0("RESULT_NoHarvDate/5.3-Peer_UseOptFact/PeerResult_UseOptFact_", dat, "_", ref, "_residuals.txt")
file.out1 <- paste0("RESULT_NoHarvDate/5.4-OutlierRemoval/Bool_Matrix_RmOut_", dat, "_", ref, ".Rdata")
file.out2 <- paste0("RESULT_NoHarvDate/5.4-OutlierRemoval/PeerResiduals_RmOut_", dat, "_", ref, ".txt")
file.out3 <- paste0("RESULT_NoHarvDate/5.4-OutlierRemoval/Stud_resid_", dat, "_", ref, ".txt")
file.out.seed <- paste0("RESULT_NoHarvDate/5.4-OutlierRemoval/RandomSeed_", dat, "_", ref, ".txt")

# load data
RawData <- fread(file.in)
data.mat <- as.matrix(RawData[, -1])

# numbers
N <- nrow(data.mat)
P <- ncol(data.mat)
threshold <- qt(p =  1 - 0.05 / (2 * N), df = (N - 2))

# do outl. detection
signif.mat <- matrix(NA, nr = N, nc = P)
rownames(signif.mat) <- RawData$V1
colnames(signif.mat) <- colnames(data.mat)
std.res.mat <- signif.mat
for (j in 1:P) {
   # j-th data
   dat <- data.frame("y" = data.mat[, j])
   
   # use asreml fit 
   set.seed(rand.seed)
   FM.asr <- asreml(fixed = y ~ 1, data = dat, trace = F,
                    na.method.X = "include", aom = TRUE, max.iter = 100)
   dfFM <- FM.asr$nedf
   
   # std.resid
   stdRes <- resid(FM.asr, type = "stdCond")
   studentizedRes <- stdRes / sqrt((dfFM - stdRes ^ 2) / (dfFM - 1))
   
   # filtering
   significance <- studentizedRes
   significance[which(abs(studentizedRes) <= threshold)] <- 'pass'
   significance[which(abs(studentizedRes) > threshold)] <- 'outlier_removed'
   
   # save
   std.res.mat[, j] <- studentizedRes
   signif.mat[, j] <- significance
}

# write this info
saveRDS(signif.mat, file.out1)

# save df
M <- as.matrix(RawData[, -1])
M[signif.mat == "outlier_removed"] <- NA
df.save <- data.frame("ID" = RawData$V1, M)
fwrite(df.save, file = file.out2)

# save df
df.save <- data.frame("ID" = RawData$V1, std.res.mat)
fwrite(df.save, file = file.out3)

# save random seed
write(rand.seed, file = file.out.seed)



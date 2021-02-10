# lib
library(data.table)
library(MASS)
library(asreml)

# object
alpha <- 0.05
rand.seed <- 2018

# mkdir
dir.save <- "RESULT/4.4-OutlierRemoval_asreml"
dir.create(dir.save, recursive = T)

# file I/O
file.in <- "RESULT/4.3-Peer_UseOptFact/PeerResult_UseOptFact_BLUE_lmer_residuals.txt"
file.out1 <- paste0(dir.save, "/Bool_Matrix_RmOut_UseOptFact_BLUE_lmer_alpha_", alpha * 100,".Rdata")
file.out2 <- paste0(dir.save, "/PeerResiduals_RmOut_UseOptFact_BLUE_lmer_alpha_", alpha * 100,".txt")
file.out3 <- paste0(dir.save, "/StudentizedRes_UseOptFact_BLUE_lmer_alpha_", alpha * 100,".txt")
file.out.seed <- paste0(dir.save, "/RandomSeed_", alpha * 100,".txt")

# load data
RawData <- fread(file.in)
data.mat <- as.matrix(RawData[, -1])
rownames(data.mat) <- RawData$V1

# numbers
N <- nrow(data.mat)
P <- ncol(data.mat)
threshold <- qt(p =  1 - alpha / (2 * N), df = (N - 2)) # Bonferroni correction

# outlier detection for each gene
studentizedRes.mat <- matrix(NA, nr = N, nc = P)
signif.mat <- matrix(NA, nr = N, nc = P)
rownames(signif.mat) <- RawData$V1
colnames(signif.mat) <- colnames(data.mat)
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
   
   # outlier removal
   significance <- studentizedRes
   significance[which(abs(studentizedRes) <= threshold)] <- 'pass'
   significance[which(abs(studentizedRes) > threshold)] <- 'outlier_removed'
   
   # save
   signif.mat[, j] <- significance
   studentizedRes.mat[, j] <- studentizedRes
}

# write this info
saveRDS(signif.mat, file.out1)

# save df
M <- as.matrix(RawData[, -1])
M[signif.mat == "outlier_removed"] <- NA
df.save <- data.frame("ID" = RawData$V1, M)
fwrite(df.save, file = file.out2)

# save df
df.save <- data.frame("ID" = RawData$V1, studentizedRes.mat)
fwrite(df.save, file = file.out3)

# save random seed
write(rand.seed, file = file.out.seed)


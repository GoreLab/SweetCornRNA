# library
library(data.table)
library(peer)

# params
args <- commandArgs(trailingOnly = T)
ref <- args[1] # "B73" or "PH207"
dat <- args[2] # "toco", "ion", "carot"
K <- as.numeric(args[3]) # number of factors
rand.seed <- 2020

# file I/O
file.input <- paste0("RESULT_NoHarvDate/4.1-MakeBlueDatasets/BLUE_matrix_for_", dat, "_", ref, ".csv")
file.out.pdf <- paste0("RESULT_NoHarvDate/5.3-Peer_UseOptFact/PeerResult_UseOptFact_", dat, "_", ref, "_plotModel.pdf")
file.out.factors <- paste0("RESULT_NoHarvDate/5.3-Peer_UseOptFact/PeerResult_UseOptFact_", dat, "_", ref, "_factors.txt")
file.out.weights <- paste0("RESULT_NoHarvDate/5.3-Peer_UseOptFact/PeerResult_UseOptFact_", dat, "_", ref, "_weights.txt")
file.out.precision <- paste0("RESULT_NoHarvDate/5.3-Peer_UseOptFact/PeerResult_UseOptFact_", dat, "_", ref, "_precision.txt")
file.out.residuals <- paste0("RESULT_NoHarvDate/5.3-Peer_UseOptFact/PeerResult_UseOptFact_", dat, "_", ref, "_residuals.txt")
file.out.seed <- paste0("RESULT_NoHarvDate/5.3-Peer_UseOptFact/PeerResult_UseOptFact_", dat, "_", ref, "_seed.txt")

# mkdir
dir.create("RESULT_NoHarvDate/5.3-Peer_UseOptFact")

# load data
ex_clean.raw <- fread(file = file.input, head = TRUE, data.table = F)
ex_clean <- ex_clean.raw[, -1]

### create the model object
set.seed(rand.seed)
model = PEER()

### Set 11 factors to model
PEER_setNk(model, K)

### Set expression data
PEER_setPhenoMean(model, as.matrix(ex_clean))
dim(PEER_getPhenoMean(model))

### Train the model, observing convergence
PEER_update(model)

### Plot the posterior variance of the factor weights and convergence diagnostics
pdf(file.out.pdf, width = 8, height = 8)
PEER_plotModel(model)
dev.off()

# name of genes & accessions
genes <- colnames(ex_clean)
samples <- as.character(ex_clean.raw[, 1])

### write PEER factors
fac <- PEER_getX(model)
fac <- as.data.frame(fac)
row.names(fac) <- samples
colnames(fac) <- paste("PEERfac_", 1:K, sep="")
fwrite(fac, file = file.out.factors, col.names = T, row.names = T, sep = "\t", quote = F)

### write PEER factors weights
wg <- PEER_getW(model)
wg <- as.data.frame(wg)
row.names(wg) <- genes
colnames(wg) <- paste("PEERfac_", 1:K, sep="")
fwrite(wg, file = file.out.weights, col.names = T, row.names = T, sep = "\t", quote = F)

### write precision (inverse variance) of the weights
pre <- PEER_getAlpha(model)
pre <- as.data.frame(pre)
pre <- cbind(paste("PEER", 1:K, sep = ""), pre)
colnames(pre) <- c("factor", "precision")
write.table(pre, file = file.out.precision, col.names = T, row.names = F, sep = "\t", quote = F)

### write residual dataset
res <- PEER_getResiduals(model)
res <- as.data.frame(res)
row.names(res) <- samples
colnames(res) <- genes
fwrite(res, file = file.out.residuals, col.names = T, row.names = T, sep = "\t", quote = F)

# save random seed
write(rand.seed, file = file.out.seed)

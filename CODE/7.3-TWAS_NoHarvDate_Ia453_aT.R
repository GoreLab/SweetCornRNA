# Use PC and kernel mutant type, according to the BIC-based model selection

# source
library(data.table)
library(rrBLUP)

# object
P3D <- FALSE # SHOULD BE "FALSE" FOR TWAS
n.core <- 10

# params
args <- commandArgs(trailingOnly = T)
ref <- "Ia453"
dat <- "toco"

# mkdir
dir.save <- "RESULT_NoHarvDate/7.1-TWAS"
dir.create(dir.save, recursive = T)

# file I/O
file.in.PeerResid <- paste0("RESULT_NoHarvDate/5.4-OutlierRemoval/PeerResiduals_RmOut_", dat, "_", ref, ".txt")
file.in.kinship <- paste0("RESULT/6.1-CalcKinship_GAPIT/KinMat_", dat, "_11K.csv")
file.in.pheno <- paste0("RESULT_NoHarvDate/4.1-MakeBlueDatasets/PhenoData_", dat, ".csv")
file.in.OptModel <- paste0("RESULT/6.2-ModelComparison_BIC/OptimalModel_BIC_", dat, ".csv")
if ( ref == "B73" ) {
   file.in.gff <- "RAWDATA/Annotation/Zea_mays.B73_RefGen_v4.59_anno.csv"
} else if ( ref == "PH207" ) {
   file.in.gff <- "RAWDATA/Annotation/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0.v03eksc_anno.csv"
}

# load data
PeerResid <- fread(file.in.PeerResid, data.table = F)
KinMat <- read.csv(file.in.kinship)
Pheno.All <- read.csv(file.in.pheno)
OptModel <- read.csv(file.in.OptModel)
if ( ref %in% c("B73", "PH207") ) { GeneInfo <- fread(file.in.gff, data.table = F) }

# kinship matrix
K <- as.matrix(KinMat[, -1])
rownames(K) <- colnames(K) <- KinMat[, 1]

# Peer residual matrix
PeerResidMat <- as.matrix(PeerResid[, -1])
myFun.01 <- function(VEC) { MIN <- min(VEC, na.rm = T); MAX <- max(VEC, na.rm = T); VEC.NEW <- (VEC - MIN) / (MAX - MIN); return(VEC.NEW)}
PeerResidMat.sc <- apply(X = PeerResidMat, 2, FUN = myFun.01)
rownames(PeerResidMat.sc) <- as.character(PeerResid$ID)
m <- match(rownames(K), rownames(PeerResidMat.sc))
PeerResidMat.sc <- PeerResidMat.sc[m, ]

# chr & pos
if ( ref %in% c("B73", "PH207") ) {
   m <- match(colnames(PeerResidMat.sc), GeneInfo$ID)
   chr <- GeneInfo$chr[m]
   pos <- GeneInfo$start[m]
   snp <- GeneInfo$ID[m] # snp = gene id
} else {
   chr <- rep(0, ncol(PeerResidMat.sc))
   pos <- rep(0, ncol(PeerResidMat.sc))
   snp <- colnames(PeerResidMat.sc)
}

# mutation
mutation.vec <- Pheno.All$Endosperm.mutation

# PC
PC <- eigen(K)$vector[, 1:10]
rownames(PC) <- rownames(K)
colnames(PC) <- paste0("PC", 1:10)

# trait
trait.all <- colnames(Pheno.All)[3:ncol(Pheno.All)]
trait <- "aT"

# optimal model
n.pc <- OptModel[OptModel$trait == trait, "n.PC"]
use.kin <- OptModel[OptModel$trait == trait, "use.kin"]
use.mu <- OptModel[OptModel$trait == trait, "use.mu"]

# filename to save result
file.TwasRes <- paste0(dir.save, "/TwasResult_", dat, "_", ref, "_", trait, ".csv")

# genotype data
geno <- data.frame("SNP" = snp, "Chr" = chr, "Pos" = pos,
                   t(PeerResidMat.sc),
                   row.names = 1:length(snp))
colnames(geno) <- c("SNP", "Chr", "Pos", rownames(PeerResidMat.sc))

# target phenotype data
if ( use.mu == "yes" ) {
   pheno <- cbind(Pheno.All[, c("Sample.ID", trait)], "mu" = mutation.vec)
} else {
   pheno <- Pheno.All[, c("Sample.ID", trait)]
}

# Run GWAS (TWAS)
b <- Sys.time()
if ( use.kin == "yes" ) {
   if ( use.mu == "yes" ) {
      print("Use kernel mutant type as covariate")
      res.GWAS <- GWAS(pheno = pheno, geno = geno, K = K, n.PC = n.pc, fixed = "mu",
                       min.MAF = -Inf, P3D = P3D, plot = FALSE, n.core = n.core)
   } else {
      res.GWAS <- GWAS(pheno = pheno, geno = geno, K = K, n.PC = n.pc,
                       min.MAF = -Inf, P3D = P3D, plot = FALSE, n.core = n.core)
   }
} else {
   pval.all <- rep(NA, nrow(geno))
   for ( k in 1:nrow(geno) ) {
      y <- pheno[, 2]
      x <- unlist(geno[k, 4:ncol(geno)])
      if ( n.pc == 0 ) {
         lmod <- lm(y ~ x)
      } else {
         lmod <- lm(y ~ x + PC[, 1:n.pc])
      }
      summary.lmod <- summary(lmod)
      pval.all[k] <- summary.lmod$coefficients["x", "Pr(>|t|)"]
   }
   res.GWAS <- data.frame("Gene" = geno$SNP,
                          "Chr" = geno$Chr,
                          "Pos" = geno$Pos,
                          "neg.log.P" = -log10(pval.all))
}
colnames(res.GWAS)[4] <- "neg.log.P"
colnames(res.GWAS)[1] <- "Gene"
a <- Sys.time()
print(a-b)

# save result
fwrite(res.GWAS, file = file.TwasRes)


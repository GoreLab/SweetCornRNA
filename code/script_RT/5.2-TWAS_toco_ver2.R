# Use PC and kernel mutant type, according to the BIC-based model selection

# source
library(data.table)
library(rrBLUP)

# object
P3D <- FALSE # SHOULD BE "FALSE" FOR TWAS
n.core <- 40

# mkdir
dir.save <- "RESULT/5.2-TWAS_toco_ver2"
dir.create(dir.save, recursive = T)

# file I/O
file.in.PeerResid <- paste0("RESULT/4.4-OutlierRemoval/PeerResiduals_RmOut_UseOptFact_BLUE_lmer_alpha_5.txt")
file.in.gff <- "RAWDATA/Annotation/Zea_mays.B73_RefGen_v4.59_anno.csv"
file.in.kinship <- "RAWDATA/SNPs_and_kinship/Kinship_from_GAPIT.csv"
file.in.pheno <- "RAWDATA/Tocochromanols/toco_transformed_blups.csv"
file.in.key <- "RAWDATA/master_key.csv"
file.in.OptModel <- "RESULT/5.1-ModelComparison_BIC/OptimalModel_BIC.csv"
filename.summary <- paste0(dir.save, "/TwasSummary.txt")

# load data
PeerResid <- fread(file.in.PeerResid, data.table = F)
KinMat <- read.csv(file.in.kinship)
Pheno.All <- read.csv(file.in.pheno)
GeneInfo <- fread(file.in.gff, data.table = F)
key.df <- read.csv(file.in.key)
OptModel <- read.csv(file.in.OptModel)

# Peer residual matrix
PeerResidMat <- as.matrix(PeerResid[, -1])
myFun.01 <- function(VEC) { MIN <- min(VEC, na.rm = T); MAX <- max(VEC, na.rm = T); VEC.NEW <- (VEC - MIN) / (MAX - MIN); return(VEC.NEW)}
PeerResidMat.sc <- apply(X = PeerResidMat, 2, FUN = myFun.01)
rownames(PeerResidMat.sc) <- as.character(PeerResid$ID)

# chr & pos
m <- match(colnames(PeerResidMat.sc), GeneInfo$ID)
chr <- GeneInfo$chr[m]
pos <- GeneInfo$start[m]
snp <- GeneInfo$ID[m] # snp = gene id

# kinship matrix
K <- as.matrix(KinMat[, -1])
rownames(K) <- colnames(K) <- KinMat[, 1]

################################################################################
# convert "accession_name" to "accession_name_toco"
setdiff(rownames(PeerResidMat.sc), key.df$accession_name) # OK!
m1 <- match(rownames(PeerResidMat.sc), key.df$accession_name)
rownames(PeerResidMat.sc) <- key.df$accession_name_toco[m1]

# remove NA (this is, indeeed, B73)
tf <- !is.na(rownames(PeerResidMat.sc))
PeerResidMat.sc <- PeerResidMat.sc[tf, ]

# match phenotype data with PEER
setdiff(rownames(PeerResidMat.sc), Pheno.All$Sample.ID) # OK!
m2 <- match(rownames(PeerResidMat.sc), Pheno.All$Sample.ID)
Pheno.Sub <- Pheno.All[m2, ]
all(Pheno.Sub$Sample.ID == rownames(PeerResidMat.sc))

# rename a few genotype of kinship matrix (manually confirm)
rownames(K)[rownames(K) == "IL731a"] <- "Il731a"
rownames(K)[rownames(K) == "IL767b"] <- "Il767b"
rownames(K)[rownames(K) == "IL777a"] <- "Il779a"
rownames(K)[rownames(K) == "IaEV191"] <- "IaEv191"
rownames(K)[rownames(K) == "IL14H"] <- "Il14H"
rownames(K)[rownames(K) == "P39a"] <- "P39A"
rownames(K)[rownames(K) == "P39a"] <- "P39A"
rownames(K)[rownames(K) == "T62s"] <- "T62S"
colnames(K) <- rownames(K)
setdiff(rownames(PeerResidMat.sc), rownames(K)) # "M6388" does not exist (fine)

# intersect of the genotype names among data
common.name <- intersect(rownames(PeerResidMat.sc), rownames(K))
m.peer <- match(common.name, rownames(PeerResidMat.sc))
PeerResidMat.sc.final <- PeerResidMat.sc[m.peer, ]
m.pheno <- match(common.name, Pheno.Sub$Sample.ID)
Pheno.Sub <- Pheno.Sub[m.pheno, ]
m.kinship <- match(common.name, rownames(K))
GRM <- K[m.kinship, m.kinship]

# check
all.equal(rownames(GRM), rownames(PeerResidMat.sc.final))
all.equal(rownames(GRM), Pheno.Sub$Sample.ID)

# Endosperm.mutation
mutation.vec <- Pheno.Sub$Endosperm.mutation

# remove unnecessary columns & modify column name
Pheno.Sub <- Pheno.Sub[, c(1, 6:ncol(Pheno.Sub))]
colnames(Pheno.Sub)[1] <- "Taxa"
################################################################################


# LOOP FOR ALL TRAITS
trait.all <- colnames(Pheno.Sub)[2:ncol(Pheno.Sub)]
for ( i in 1:length(trait.all) ) {
   # trait
   trait <- trait.all[i]
   
   # print 
   print(paste0("Start TWAS for ", trait, ": ", i, "-th trait out of ", length(trait.all), " traits."))
   
   # optimal model
   n.pc <- OptModel[OptModel$trait == trait, "n.PC"]
   use.kin <- OptModel[OptModel$trait == trait, "use.kin"]
   use.mu <- OptModel[OptModel$trait == trait, "use.mu"]
   
   # filename to save result
   file.TwasRes <- paste0(dir.save, "/TwasResult_", trait, ".csv")
   
   # genotype data
   geno <- data.frame("SNP" = snp, "Chr" = chr, "Pos" = pos,
                      t(PeerResidMat.sc.final),
                      row.names = 1:length(snp))
   colnames(geno) <- c("SNP", "Chr", "Pos", rownames(PeerResidMat.sc.final))
   
   # target phenotype data
   if ( use.mu == "yes" ) {
      pheno <- cbind(Pheno.Sub[, c("Taxa", trait)], "mu" = mutation.vec)
   } else {
      pheno <- Pheno.Sub[, c("Taxa", trait)]
   }
   
   # Run GWAS (TWAS)
   b <- Sys.time()
   if ( use.kin == "yes" ) {
      if ( use.mu == "yes" ) {
         print("Use kernel mutant type as covariate")
         res.GWAS <- GWAS(pheno = pheno, geno = geno, K = GRM, n.PC = n.pc, fixed = "mu",
                          min.MAF = -Inf, P3D = P3D, plot = FALSE, n.core = n.core)
      } else {
         res.GWAS <- GWAS(pheno = pheno, geno = geno, K = GRM, n.PC = n.pc,
                          min.MAF = -Inf, P3D = P3D, plot = FALSE, n.core = n.core)
      }
   } else {
      print("Error: GWAS/TWAS without kinship matrix is not implemented now.")
   }
   colnames(res.GWAS)[4] <- "neg.log.P"
   colnames(res.GWAS)[1] <- "Gene"
   a <- Sys.time()
   print(a-b)
   
   # save result
   fwrite(res.GWAS, file = file.TwasRes)
}

# write summary
text.b <- paste0("Phenotype data has ", nrow(Pheno.All), " accessions")
text.c <- paste0("Peer residual matrix has ", nrow(PeerResidMat.sc), " accessions and ", ncol(PeerResidMat.sc), " genes")
text.d <- paste0("Number of overlapped accessions = ", length(common.name))
SummaryText <- c(text.b, text.c, text.d)
write(x = SummaryText, file = filename.summary, sep = "/n")


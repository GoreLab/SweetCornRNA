# source
library(data.table)
library(BGLR)

# params
args <- commandArgs(trailingOnly = T)
model <- args[1]
method <- args[2]
n.fold <- 5
n.rep <- 10
nIter <- 12000
burnIn <- 8000

# mkdir
dir.save <- "RESULT/6.1-Prediction_toco"
dir.create(dir.save, recursive = T)

# file I/O
file.in.PeerResid <- paste0("RESULT/4.3-Peer_UseOptFact/PeerResult_UseOptFact_", model, "_", method, "_residuals.txt")
file.in.gff <- "RAWDATA/Annotation/Zea_mays.B73_RefGen_v4.59_anno.csv"
file.in.kinship <- "RAWDATA/SNPs_and_kinship/Kinship_from_GAPIT.csv"
file.in.pheno <- "RAWDATA/Tocochromanols/toco_transformed_blups.csv"
file.in.key <- "RAWDATA/master_key.csv"
filename.summary <- paste0(dir.save, "/TwasSummary_", model, "_", method, ".txt")

# load data
PeerResid <- fread(file.in.PeerResid, data.table = F)
KinMat <- read.csv(file.in.kinship)
Pheno.All <- read.csv(file.in.pheno)
GeneInfo <- fread(file.in.gff, data.table = F)
key.df <- read.csv(file = file.in.key)

# Peer residual matrix
PeerResidMat <- as.matrix(PeerResid[, -1])
myFun.01 <- function(VEC) { MIN <- min(VEC, na.rm = T); MAX <- max(VEC, na.rm = T); VEC.NEW <- (VEC - MIN) / (MAX - MIN); return(VEC.NEW)}
PeerResidMat.sc <- apply(X = PeerResidMat, 2, FUN = myFun.01)
rownames(PeerResidMat.sc) <- as.character(PeerResid$V1)

# TBD: remove genes without any variation
tf <- apply(is.na(PeerResidMat.sc), 2, sum) == 0
PeerResidMat.sc <- PeerResidMat.sc[, tf]

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

# remove unnecessary columns & modify column name
Pheno.Sub <- Pheno.Sub[, c(1, 6:ncol(Pheno.Sub))]
colnames(Pheno.Sub)[1] <- "Taxa"
################################################################################


# make a tGRM
PeerResidMat.sc.final.sc <- scale(PeerResidMat.sc.final)
tGRM <- tcrossprod(PeerResidMat.sc.final.sc) / ncol(PeerResidMat.sc.final.sc) # tGRM

# set model
ETA.mGRM <- list("G" = list(K = GRM, model = "RKHS"))
ETA.tGRM <- list("T" = list("K" = tGRM, "model" = "RKHS"))
ETA.mGRM.tGRM <- list("G" = list("K" = GRM, "model" = "RKHS"),
                      "T" = list("K" = tGRM, "model" = "RKHS"))

# make a cross validation split
set.seed(2020)
df.CV <- data.frame("accession_name_toco" = Pheno.Sub$Taxa)
mat <- matrix(NA, nr = nrow(df.CV), nc = n.rep)
colnames(mat) <- paste0("Rep", formatC(1:n.rep, width = 2, flag = "0"))
for ( i in 1:n.rep ) {
   mat[, i] <- sample(rep(1:n.fold, length.out = nrow(df.CV)))
}
df.CV <- cbind(df.CV, mat)
write.csv(df.CV, file = paste0(dir.save, "/CrossValidationFold.csv"), row.names = F)

# mkdir to save log file (not so important)
dir.log <- paste0(dir.save, "/LOGFILE")
dir.create(dir.log, recursive = T)

# define a function to do cross validation
myFun.Cv <- function(y.all = y.all, ETA = ETA.mGRM, df.CV = df.CV,
                     nIter = 12000, burnIn = 8000,
                     name.log = paste0(dir.log, "/fm_"), ...) {
   # get numbers
   n.rep <- ncol(df.CV) - 1
   n.fold <- max(df.CV[, 2])
   
   # cross validation: n.rep/n.fold 
   pred.mat <- matrix(NA, nr = nrow(df.CV), nc = n.rep)
   rownames(pred.mat) <- df.CV[, 1]
   colnames(pred.mat) <- colnames(df.CV)[2:ncol(df.CV)]
   for ( j in 1:n.rep ) {
      # make an object to save prediction result
      pred.vec <- rep(NA, length(y.all))
      
      for ( k in 1:n.fold ) {
         # cv numbers
         cv.num <- df.CV[, j + 1]
         
         # cross validation
         y.obs <- y.all
         y.obs[cv.num == k] <- NA # mask phenotype
         
         # run prediction
         fm <- BGLR(y = y.obs, ETA = ETA.mGRM,
                    nIter = nIter, burnIn = burnIn, verbose = F,
                    saveAt = name.log)
         pred.vec[cv.num == k] <- fm$yHat[cv.num == k]
         
         # print
         comment <- paste0("Cross Validation: ", k, "-th fold of ", n.fold, " folds, ",
                           j, "-th rep of ", n.rep, " replications.")
         print(comment)
      }
      
      # save
      pred.mat[, j] <- pred.vec
   }
   
   # output
   return(pred.mat)
}

# LOOP FOR ALL TRAITS with different models
trait.all <- colnames(Pheno.Sub)[2:ncol(Pheno.Sub)]
for ( i in 1:length(trait.all) ) {
   # trait
   trait <- trait.all[i]
   
   # target phenotype data
   y.all <- Pheno.Sub[, trait]
   
   # Run
   pred.mat.mGRM <- myFun.Cv(y.all = y.all, ETA = ETA.mGRM, df.CV = df.CV,
                             nIter = nIter, burnIn = burnIn,
                             name.log = paste0(dir.log, "/fm_"))
   pred.mat.tGRM <- myFun.Cv(y.all = y.all, ETA = ETA.tGRM, df.CV = df.CV,
                             nIter = nIter, burnIn = burnIn,
                             name.log = paste0(dir.log, "/fm_"))
   pred.mat.mGRM.tGRM <- myFun.Cv(y.all = y.all, ETA = ETA.mGRM.tGRM, df.CV = df.CV,
                                  nIter = nIter, burnIn = burnIn,
                                  name.log = paste0(dir.log, "/fm_"))
   
   # make a data frame of CV result & save it
   df.CV.res.mGRM <- data.frame("accession_name_toco" = df.CV[, 1], pred.mat.mGRM)
   df.CV.res.tGRM <- data.frame("accession_name_toco" = df.CV[, 1], pred.mat.tGRM)
   df.CV.res.mGRM.tGRM <- data.frame("accession_name_toco" = df.CV[, 1], pred.mat.mGRM.tGRM)
   write.csv(df.CV.res.mGRM, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", trait, "_mGRM_", 
                           model, "_", method, ".csv"))
   write.csv(df.CV.res.tGRM, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", trait, "_tGRM_", 
                           model, "_", method, ".csv"))
   write.csv(df.CV.res.mGRM.tGRM, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", trait, "_mGRM.tGRM_", 
                           model, "_", method, ".csv"))
}


# df.fig <- data.frame("observed" = Pheno.Sub[, trait], "predicted" = df.CV.res.tGRM[, 2])
# lim <- range(c(df.fig$observed, df.fig$predicted), na.rm = T)
# library(ggplot2)
# p <- ggplot(df.fig, aes(x = observed, y = predicted))
# p <- p + geom_point() + xlim(lim) + ylim(lim)
# p <- p + geom_abline(slope = 1, intercept = 0, lty = 2)
# p <- p + ggtitle(paste0("Predcition for ", trait, " by using tGRM"))
# p
# 
# cor.vec.m <- as.numeric(cor(Pheno.Sub[, trait], df.CV.res.mGRM[, 2:ncol(df.CV.res.mGRM)]))
# cor.vec.t <- as.numeric(cor(Pheno.Sub[, trait], df.CV.res.tGRM[, 2:ncol(df.CV.res.tGRM)]))
# cor.vec.mt <- as.numeric(cor(Pheno.Sub[, trait], df.CV.res.mGRM.tGRM[, 2:ncol(df.CV.res.mGRM.tGRM)]))
# 
# cor.vec.m; cor.vec.t; cor.vec.mt
# 
# t.test(cor.vec.m, cor.vec.t)
# t.test(cor.vec.m, cor.vec.mt)


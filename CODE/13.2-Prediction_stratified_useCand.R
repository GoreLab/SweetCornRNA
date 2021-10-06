# source
library(data.table)
library(BGLR)

# params
n.fold <- 5
n.rep <- 10
nIter <- 12000
burnIn <- 8000

# mkdir
dir.save <- "RESULT_NoHarvDate/13.2-Prediction_stratified_useCand"
dir.create(dir.save, recursive = T)

# arguments
args <- commandArgs(trailingOnly = T)
dat <- args[1] # "toco", "ion", or "carot"
mu <- args[2] # "UseMu" or "NoMu"

# ---------------------------------------------------------------------------- #
# BLUE-B73
f <- paste0("RESULT_NoHarvDate/4.1-MakeBlueDatasets/BLUE_matrix_for_", dat, "_B73.csv")
ExprDat.B73 <- fread(f, data.table = F)
ExprMat.B73 <- as.matrix(ExprDat.B73[, -1])
rownames(ExprMat.B73) <- ExprDat.B73[, 1]

# kinship
f <- paste0("RESULT/6.1-CalcKinship_GAPIT/KinMat_", dat, "_163K.csv")
KinDat <- read.csv(f)
KinMat <- as.matrix(KinDat[, -1])
rownames(KinMat) <- colnames(KinMat) <- KinDat[, 1]

# phenotype
f <- paste0("RESULT_NoHarvDate/4.1-MakeBlueDatasets/PhenoData_", dat, ".csv")
PhenoDat <- fread(f, data.table = F)

# candidate gene
f <- "RAWDATA/Candidate_Genes/PlantCell_genes.csv"
PcGenesDat <- fread(f, data.table = F)


# ---------------------------------------------------------------------------- #
# ad hoc solution to the garbled characters...
colnames(PcGenesDat) <- gsub("﻿", "", colnames(PcGenesDat))
for ( i in 1:ncol(PcGenesDat) ) {
   PcGenesDat[, i] <- gsub("﻿", "", PcGenesDat[, i])
}


# ---------------------------------------------------------------------------- #
# PC genes present in our data
tf <- PcGenesDat$`RefGen_v4 ID` %in% colnames(ExprMat.B73)
PcGenesDat.present <- PcGenesDat[tf, ]

# gene id to be used
pc.genes <- PcGenesDat.present$`RefGen_v4 ID`[PcGenesDat.present$Pathway == dat]
non.pc.genes <- setdiff(colnames(ExprMat.B73), pc.genes)


# ---------------------------------------------------------------------------- #
# make a tGRM 
myFun.Make.tGRM <- function(M) {
   M.sc <- scale(M)
   tf <- is.na(M.sc)
   if ( sum(tf) != 0 ) { M.sc[tf] <- 0 } # missing value = impute by mean (= zero)
   tGRM <- tcrossprod(M.sc) / ncol(M.sc) # tGRM = linear kernel
   return(tGRM)
}
tGRM.PC.genes <- myFun.Make.tGRM(M = ExprMat.B73[, pc.genes])
tGRM.non.PC.genes <- myFun.Make.tGRM(M = ExprMat.B73[, non.pc.genes])


# ---------------------------------------------------------------------------- #
if ( mu == "UseMu" ) {
   # design matrix for endosperm mutation
   mutation.vec <- PhenoDat$Endosperm.mutation
   design.X <- model.matrix(~ -1 + mutation.vec)
   colnames(design.X) <- gsub("mutation.vec", "", colnames(design.X))
   
   # set models
   ETA.tGRM.non.cand <- list("T" = list("K" = tGRM.non.PC.genes, "model" = "RKHS"),
                             "Mu" = list(X = design.X, model = "FIXED"))
   ETA.tGRM.cand <- list("T" = list("K" = tGRM.PC.genes, "model" = "RKHS"),
                         "Mu" = list(X = design.X, model = "FIXED"))
   ETA.tGRM.both <- list("T1" = list("K" = tGRM.PC.genes, "model" = "RKHS"),
                         "T2" = list("K" = tGRM.non.PC.genes, "model" = "RKHS"),
                         "Mu" = list(X = design.X, model = "FIXED"))
   ETA.mGRM.tGRM.non.cand <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                                  "T" = list("K" = tGRM.non.PC.genes, "model" = "RKHS"),
                                  "Mu" = list(X = design.X, model = "FIXED"))
   ETA.mGRM.tGRM.cand <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                              "T" = list("K" = tGRM.PC.genes, "model" = "RKHS"),
                              "Mu" = list(X = design.X, model = "FIXED"))
   ETA.mGRM.tGRM.both <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                              "T1" = list("K" = tGRM.PC.genes, "model" = "RKHS"),
                              "T2" = list("K" = tGRM.non.PC.genes, "model" = "RKHS"),
                              "Mu" = list(X = design.X, model = "FIXED"))
} else if ( mu == "NoMu" ) {
   # set models
   ETA.tGRM.non.cand <- list("T" = list("K" = tGRM.non.PC.genes, "model" = "RKHS"))
   ETA.tGRM.cand <- list("T" = list("K" = tGRM.PC.genes, "model" = "RKHS"))
   ETA.tGRM.both <- list("T1" = list("K" = tGRM.PC.genes, "model" = "RKHS"),
                         "T2" = list("K" = tGRM.non.PC.genes, "model" = "RKHS"))
   ETA.mGRM.tGRM.non.cand <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                                  "T" = list("K" = tGRM.non.PC.genes, "model" = "RKHS"))
   ETA.mGRM.tGRM.cand <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                              "T" = list("K" = tGRM.PC.genes, "model" = "RKHS"))
   ETA.mGRM.tGRM.both <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                              "T1" = list("K" = tGRM.PC.genes, "model" = "RKHS"),
                              "T2" = list("K" = tGRM.non.PC.genes, "model" = "RKHS"))
} else {
   print("Error: Incorrect specification for the use of endosperm mutant")
}


# ---------------------------------------------------------------------------- #
# Use a cross validation split
f <- paste0("RESULT_NoHarvDate/13.1-Prediction_stratified/CrossValidationFold_", dat, "_", mu, ".csv")
df.CV.save <- read.csv(f)
df.CV <- df.CV.save[, -2]

# ---------------------------------------------------------------------------- #
# mkdir to save log file (not so important)
dir.log <- paste0(dir.save, "/LOGFILE")
dir.create(dir.log, recursive = T)


# ---------------------------------------------------------------------------- #
# define a function to do cross validation
myFun.Cv <- function(y.all, ETA, df.CV,
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
         fm <- BGLR(y = y.obs, ETA = ETA,
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


# ---------------------------------------------------------------------------- #
# LOOP FOR ALL TRAITS with different models
trait.all <- colnames(PhenoDat)[3:ncol(PhenoDat)]
for ( i in 1:length(trait.all) ) {
   # trait
   trait <- trait.all[i]
   
   # target phenotype data
   y.all <- PhenoDat[, trait]
   
   # Run
   pred.mat.tGRM.non.cand <- myFun.Cv(y.all = y.all, 
                                      ETA = ETA.tGRM.non.cand, 
                                      df.CV = df.CV,
                                      nIter = nIter, burnIn = burnIn,
                                      name.log = paste0(dir.log, "/fm_"))
   pred.mat.tGRM.cand <- myFun.Cv(y.all = y.all, 
                                  ETA = ETA.tGRM.cand, 
                                  df.CV = df.CV,
                                  nIter = nIter, burnIn = burnIn,
                                  name.log = paste0(dir.log, "/fm_"))
   pred.mat.tGRM.both <- myFun.Cv(y.all = y.all, 
                                  ETA = ETA.tGRM.both, 
                                  df.CV = df.CV,
                                  nIter = nIter, burnIn = burnIn,
                                  name.log = paste0(dir.log, "/fm_"))
   pred.mat.mGRM.tGRM.non.cand <- myFun.Cv(y.all = y.all, 
                                           ETA = ETA.mGRM.tGRM.non.cand, 
                                           df.CV = df.CV,
                                           nIter = nIter, burnIn = burnIn,
                                           name.log = paste0(dir.log, "/fm_"))
   pred.mat.mGRM.tGRM.cand <- myFun.Cv(y.all = y.all, 
                                       ETA = ETA.mGRM.tGRM.cand, 
                                       df.CV = df.CV,
                                       nIter = nIter, burnIn = burnIn,
                                       name.log = paste0(dir.log, "/fm_"))
   pred.mat.mGRM.tGRM.both <- myFun.Cv(y.all = y.all, 
                                       ETA = ETA.mGRM.tGRM.both, 
                                       df.CV = df.CV,
                                       nIter = nIter, burnIn = burnIn,
                                       name.log = paste0(dir.log, "/fm_"))
   
   # make a data frame of CV result & save it
   df.CV.res.tGRM.non.cand <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.tGRM.non.cand)
   df.CV.res.tGRM.cand <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.tGRM.cand)
   df.CV.res.tGRM.both <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.tGRM.both)
   df.CV.res.mGRM.tGRM.non.cand <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.mGRM.tGRM.non.cand)
   df.CV.res.mGRM.tGRM.cand <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.mGRM.tGRM.cand)
   df.CV.res.mGRM.tGRM.both <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.mGRM.tGRM.both)
   write.csv(df.CV.res.tGRM.non.cand, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_tGRM.non.cand.csv"))
   write.csv(df.CV.res.tGRM.cand, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_tGRM.cand.csv"))
   write.csv(df.CV.res.tGRM.both, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_tGRM.both.csv"))
   write.csv(df.CV.res.mGRM.tGRM.non.cand, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_mGRM.tGRM.non.cand.csv"))
   write.csv(df.CV.res.mGRM.tGRM.cand, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_mGRM.tGRM.cand.csv"))
   write.csv(df.CV.res.mGRM.tGRM.both, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_mGRM.tGRM.both.csv"))
   
}

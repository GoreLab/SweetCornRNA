# source
library(data.table)
library(BGLR)

# params
n.fold <- 5
n.rep <- 10
nIter <- 12000
burnIn <- 8000

# mkdir
dir.save <- "RESULT_NoHarvDate/13.1-Prediction_stratified"
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

# BLUE-PH207
f <- paste0("RESULT_NoHarvDate/4.1-MakeBlueDatasets/BLUE_matrix_for_", dat, "_PH207.csv")
ExprDat.PH207 <- fread(f, data.table = F)
ExprMat.PH207 <- as.matrix(ExprDat.PH207[, -1])
rownames(ExprMat.PH207) <- ExprDat.PH207[, 1]

# BLUE-Ia453
f <- paste0("RESULT_NoHarvDate/4.1-MakeBlueDatasets/BLUE_matrix_for_", dat, "_Ia453.csv")
ExprDat.Ia453 <- fread(f, data.table = F)
ExprMat.Ia453 <- as.matrix(ExprDat.Ia453[, -1])
rownames(ExprMat.Ia453) <- ExprDat.Ia453[, 1]
m <- match(rownames(ExprMat.PH207), rownames(ExprMat.Ia453))
ExprMat.Ia453 <- ExprMat.Ia453[m, ] # sort

# kinship
f <- paste0("RESULT/6.1-CalcKinship_GAPIT/KinMat_", dat, "_163K.csv")
KinDat <- read.csv(f)
KinMat <- as.matrix(KinDat[, -1])
rownames(KinMat) <- colnames(KinMat) <- KinDat[, 1]

# phenotype
f <- paste0("RESULT_NoHarvDate/4.1-MakeBlueDatasets/PhenoData_", dat, ".csv")
PhenoDat <- fread(f, data.table = F)


# ---------------------------------------------------------------------------- #
# make a tGRM 
myFun.Make.tGRM <- function(M) {
   M.sc <- scale(M)
   tf <- is.na(M.sc)
   if ( sum(tf) != 0 ) { M.sc[tf] <- 0 } # missing value = impute by mean (= zero)
   tGRM <- tcrossprod(M.sc) / ncol(M.sc) # tGRM = linear kernel
   return(tGRM)
}
tGRM.B73 <- myFun.Make.tGRM(M = ExprMat.B73)
tGRM.PH207 <- myFun.Make.tGRM(M = ExprMat.PH207)
tGRM.Ia453 <- myFun.Make.tGRM(M = ExprMat.Ia453)

# ---------------------------------------------------------------------------- #
if ( mu == "UseMu" ) {
   # design matrix for endosperm mutation
   mutation.vec <- PhenoDat$Endosperm.mutation
   design.X <- model.matrix(~ -1 + mutation.vec)
   colnames(design.X) <- gsub("mutation.vec", "", colnames(design.X))
   
   # set models
   ETA.mGRM <- list("G" = list(K = KinMat, model = "RKHS"),
                    "Mu" = list(X = design.X, model = "FIXED"))
   ETA.tGRM.B73 <- list("T" = list("K" = tGRM.B73, "model" = "RKHS"),
                        "Mu" = list(X = design.X, model = "FIXED"))
   ETA.tGRM.PH207 <- list("T" = list("K" = tGRM.PH207, "model" = "RKHS"),
                          "Mu" = list(X = design.X, model = "FIXED"))
   ETA.tGRM.Ia453 <- list("T" = list("K" = tGRM.Ia453, "model" = "RKHS"),
                          "Mu" = list(X = design.X, model = "FIXED"))
   ETA.mGRM.tGRM.B73 <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                             "T" = list("K" = tGRM.B73, "model" = "RKHS"),
                             "Mu" = list(X = design.X, model = "FIXED"))
   ETA.mGRM.tGRM.PH207 <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                               "T" = list("K" = tGRM.PH207, "model" = "RKHS"),
                               "Mu" = list(X = design.X, model = "FIXED"))
   ETA.mGRM.tGRM.Ia453 <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                               "T" = list("K" = tGRM.Ia453, "model" = "RKHS"),
                               "Mu" = list(X = design.X, model = "FIXED"))
} else if ( mu == "NoMu" ) {
   # set models
   ETA.mGRM <- list("G" = list(K = KinMat, model = "RKHS"))
   ETA.tGRM.B73 <- list("T" = list("K" = tGRM.B73, "model" = "RKHS"))
   ETA.tGRM.PH207 <- list("T" = list("K" = tGRM.PH207, "model" = "RKHS"))
   ETA.tGRM.Ia453 <- list("T" = list("K" = tGRM.Ia453, "model" = "RKHS"))
   ETA.mGRM.tGRM.B73 <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                             "T" = list("K" = tGRM.B73, "model" = "RKHS"))
   ETA.mGRM.tGRM.PH207 <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                               "T" = list("K" = tGRM.PH207, "model" = "RKHS"))
   ETA.mGRM.tGRM.Ia453 <- list("G" = list("K" = KinMat, "model" = "RKHS"),
                               "T" = list("K" = tGRM.Ia453, "model" = "RKHS"))
} else {
   print("Error: Incorrect specification for the use of endosperm mutant")
}

# ---------------------------------------------------------------------------- #
# make a cross validation split
set.seed(2020)
df.CV.raw <- data.frame("Sample.ID" = PhenoDat$Sample.ID)
mat <- matrix(NA, nr = nrow(df.CV.raw), nc = n.rep)
colnames(mat) <- paste0("Rep", formatC(1:n.rep, width = 2, flag = "0"))
for ( r in 1:n.rep ) {
   n.tab <- table(PhenoDat$Endosperm.mutation)
   tmp.n.tab <- cumsum(c(0, n.tab))
   cv.num.with.names.all <- c()
   for ( i in 1:length(n.tab) ) {
      mu.type <- names(n.tab)[i]
      sample.names <- PhenoDat$Sample.ID[PhenoDat$Endosperm.mutation == mu.type]
      seq.num <- (tmp.n.tab[i]+1):(tmp.n.tab[i+1])
      cv.num <- 1 + seq.num %% n.fold
      cv.num.rand <- sample(cv.num)
      cv.num.with.names <- setNames(cv.num.rand, sample.names)
      cv.num.with.names.all <- c(cv.num.with.names.all, cv.num.with.names)
   }
   cv.num.with.names.all.ord <- cv.num.with.names.all[PhenoDat$Sample.ID] # sort
   mat[, r] <- cv.num.with.names.all.ord
}
df.CV <- cbind(df.CV.raw, mat)
df.CV.save <- cbind(df.CV.raw, "Endosperm.mutation" = PhenoDat$Endosperm.mutation, mat)
f <- paste0(dir.save, "/CrossValidationFold_", dat, "_", mu, ".csv")
write.csv(df.CV.save, file = f, row.names = F)

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
   pred.mat.mGRM <- myFun.Cv(y.all = y.all, ETA = ETA.mGRM, df.CV = df.CV,
                             nIter = nIter, burnIn = burnIn,
                             name.log = paste0(dir.log, "/fm_"))
   pred.mat.tGRM.B73 <- myFun.Cv(y.all = y.all, ETA = ETA.tGRM.B73, df.CV = df.CV,
                                 nIter = nIter, burnIn = burnIn,
                                 name.log = paste0(dir.log, "/fm_"))
   pred.mat.tGRM.PH207 <- myFun.Cv(y.all = y.all, ETA = ETA.tGRM.PH207, df.CV = df.CV,
                                   nIter = nIter, burnIn = burnIn,
                                   name.log = paste0(dir.log, "/fm_"))
   pred.mat.tGRM.Ia453 <- myFun.Cv(y.all = y.all, ETA = ETA.tGRM.Ia453, df.CV = df.CV,
                                   nIter = nIter, burnIn = burnIn,
                                   name.log = paste0(dir.log, "/fm_"))
   pred.mat.mGRM.tGRM.B73 <- myFun.Cv(y.all = y.all, ETA = ETA.mGRM.tGRM.B73, df.CV = df.CV,
                                      nIter = nIter, burnIn = burnIn,
                                      name.log = paste0(dir.log, "/fm_"))
   pred.mat.mGRM.tGRM.PH207 <- myFun.Cv(y.all = y.all, ETA = ETA.mGRM.tGRM.PH207, df.CV = df.CV,
                                        nIter = nIter, burnIn = burnIn,
                                        name.log = paste0(dir.log, "/fm_"))
   pred.mat.mGRM.tGRM.Ia453 <- myFun.Cv(y.all = y.all, ETA = ETA.mGRM.tGRM.Ia453, df.CV = df.CV,
                                        nIter = nIter, burnIn = burnIn,
                                        name.log = paste0(dir.log, "/fm_"))
   
   # make a data frame of CV result & save it
   df.CV.res.mGRM <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.mGRM)
   df.CV.res.tGRM.B73 <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.tGRM.B73)
   df.CV.res.tGRM.PH207 <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.tGRM.PH207)
   df.CV.res.tGRM.Ia453 <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.tGRM.Ia453)
   df.CV.res.mGRM.tGRM.B73 <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.mGRM.tGRM.B73)
   df.CV.res.mGRM.tGRM.PH207 <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.mGRM.tGRM.PH207)
   df.CV.res.mGRM.tGRM.Ia453 <- data.frame("Sample.ID" = df.CV[, 1], pred.mat.mGRM.tGRM.Ia453)
   write.csv(df.CV.res.mGRM, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_mGRM.csv"))
   write.csv(df.CV.res.tGRM.B73, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_tGRM.B73.csv"))
   write.csv(df.CV.res.tGRM.PH207, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_tGRM.PH207.csv"))
   write.csv(df.CV.res.tGRM.Ia453, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_tGRM.Ia453.csv"))
   write.csv(df.CV.res.mGRM.tGRM.B73, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_mGRM.tGRM.B73.csv"))
   write.csv(df.CV.res.mGRM.tGRM.PH207, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_mGRM.tGRM.PH207.csv"))
   write.csv(df.CV.res.mGRM.tGRM.Ia453, row.names = F,
             file = paste0(dir.save, "/PredictionResult_", dat, "_", trait, "_", mu, "_mGRM.tGRM.Ia453.csv"))
}

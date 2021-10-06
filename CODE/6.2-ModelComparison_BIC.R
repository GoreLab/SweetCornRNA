# source
library(data.table)
library(rrBLUP)
library(ggplot2)

# params
args <- commandArgs(trailingOnly = T)
dat <- args[1] # "toco", "ion", "carot"

# mkdir
dir.save <- "RESULT/6.2-ModelComparison_BIC"
dir.create(dir.save, recursive = T)

# file I/O
file.in.kinship <- paste0("RESULT/6.1-CalcKinship_GAPIT/KinMat_", dat, "_11K.csv")
file.in.pheno <- paste0("RESULT/4.1-MakeBlueDatasets/PhenoData_", dat, ".csv")

# load data
KinMat <- read.csv(file.in.kinship)
Pheno.All <- read.csv(file.in.pheno)

# kinship matrix
K <- as.matrix(KinMat[, -1])
rownames(K) <- colnames(K) <- KinMat[, 1]

# check -> OK
all(Pheno.All$Sample.ID == rownames(K)) 

# PC
PC <- eigen(K)$vector[, 1:10]
rownames(PC) <- rownames(K)
colnames(PC) <- paste0("PC", 1:10)

# make a "full" design matrix for fixed effects
mutation.vec <- Pheno.All$Endosperm.mutation
design.X.mu <- model.matrix(~ -1 + mutation.vec)
design.X.full <- cbind(1, PC, design.X.mu)
colnames(design.X.full) <- c("intercept", colnames(PC), "sh2", "su1", "su1sh2")

# define a function to get a design matrix
myFun.GetDesignMat <- function(n.PC, use.mu, ...) {
  if ( n.PC == 0 & use.mu == "no" ) {  design.X <- design.X.full[, "intercept", drop = F] }
  if ( n.PC == 0 & use.mu == "yes" ) {  design.X <- design.X.full[, c("sh2", "su1", "su1sh2"), drop = F] }
  if ( n.PC >= 1 & use.mu == "no" ) {  design.X <- design.X.full[, c("intercept", paste0("PC", 1:n.PC)), drop = F] }
  if ( n.PC >= 1 & use.mu == "yes" ) {  design.X <- design.X.full[, c(paste0("PC", 1:n.PC), c("sh2", "su1", "su1sh2")), drop = F] }
  return(design.X)
}

# covariates
df.covariate <- expand.grid("n.PC" = 0:10, "use.mu" = c("no", "yes"))

# define a function to calculate LL
myFun.CalcLogLik <- function(ms, y, K, X) {
  n <- length(y)
  q <- length(ms$beta)
  H <- K + diag(nrow(K)) * (ms$Ve / ms$Vu)
  H.inv <- MASS::ginv(H)
  a <- (-1) * n * log(2 * pi * ms$Vu)
  b <- (-1) * determinant(H, logarithm = T)$modulus[1]
  c <- as.numeric((-1) * (1 / ms$Vu) * t(y - X %*% ms$beta) %*% H.inv %*% (y - X %*% ms$beta))
  LogLik.ML <- 0.5 * (a + b + c) # full likelihood based on REML 
  return(LogLik.ML)
}



# LOOP FOR ALL TRAITS
trait.all <- colnames(Pheno.All)[3:ncol(Pheno.All)]
res.BIC.all <- NULL
for ( i in 1:length(trait.all) ) {
  # target phenotype data
  trait <- trait.all[i]
  pheno <- Pheno.All[, trait]
  
  # remove NA
  tf <- !is.na(pheno)
  pheno.i <- pheno[tf]
  GRM.i <- K[tf, tf]
  mutation.vec.i <- mutation.vec[tf]
  PC.i <- PC[tf, ]
  
  # BIC with kinship
  res.BIC <- df.covariate
  res.BIC$trait <- trait
  res.BIC$use.kin <- "yes"
  res.BIC$LL <- res.BIC$BIC <- NA
  for ( j in 1:nrow(df.covariate) ) {
    n.PC <- df.covariate$n.PC[j]
    use.mu <- df.covariate$use.mu[j]
    design.X <- myFun.GetDesignMat(n.PC, use.mu)
    design.X.i <- design.X[tf, , drop = F]
    ms <- mixed.solve(y = pheno.i, K = GRM.i, X = design.X.i, method = "REML")
    LL <- myFun.CalcLogLik(ms, y = pheno.i, K = GRM.i, X = design.X.i)
    npar <- length(ms$beta) + 2
    BIC <- npar * log(length(pheno.i)) - 2 * LL
    res.BIC$LL[j] <- LL; res.BIC$BIC[j] <- BIC
  }
  res.BIC.with.kinship <- res.BIC
  
  # BIC without kinship
  res.BIC <- df.covariate
  res.BIC$trait <- trait
  res.BIC$use.kin <- "no"
  res.BIC$LL <- res.BIC$BIC <- NA
  for ( j in 1:nrow(df.covariate) ) {
    n.PC <- df.covariate$n.PC[j]
    use.mu <- df.covariate$use.mu[j]
    if ( n.PC == 0 & use.mu == "no" ) { lmod <- lm(pheno.i ~ 1); npar <- 2 }
    if ( n.PC == 0 & use.mu == "yes" ) { lmod <- lm(pheno.i ~ -1 + mutation.vec.i); npar <- 4 }
    if ( n.PC >= 1 & use.mu == "no" ) { lmod <- lm(pheno.i ~ 1 + PC.i[, 1:n.PC]); npar <- n.PC + 2 }
    if ( n.PC >= 1 & use.mu == "yes" ) { lmod <- lm(pheno.i ~ -1 + PC.i[, 1:n.PC] + mutation.vec.i); npar <- n.PC + 4 }
    LL <- logLik(lmod)
    BIC <- npar * log(length(pheno.i)) - 2 * LL
    res.BIC$LL[j] <- LL; res.BIC$BIC[j] <- BIC
  }
  res.BIC.no.kinship <- res.BIC
  
  # merge results 
  res.BIC <- rbind(res.BIC.with.kinship, res.BIC.no.kinship)
  
  # make a figure
  res.BIC$use.mu <- as.character(res.BIC$use.mu)
  res.BIC$use.mu[res.BIC$use.mu == "no"] <- "Without Endospearm Mutant"
  res.BIC$use.mu[res.BIC$use.mu == "yes"] <- "Use Endospearm Mutant"
  res.BIC$use.mu <- factor(res.BIC$use.mu, levels = c("Without Endospearm Mutant", "Use Endospearm Mutant"))
  res.BIC$use.kin <- as.character(res.BIC$use.kin)
  res.BIC$use.kin[res.BIC$use.kin == "no"] <- "Without Kinship"
  res.BIC$use.kin[res.BIC$use.kin == "yes"] <- "Use Kinship"
  res.BIC$use.kin <- factor(res.BIC$use.kin, levels = c("Without Kinship", "Use Kinship"))
  num.min <- which.min(res.BIC$BIC)
  res.BIC$use.kin <- as.factor(res.BIC$use.kin)
  p <- ggplot(res.BIC, aes(x = n.PC, y = BIC))
  p <- p + geom_line()
  p <- p + geom_point(size = 2)
  p <- p + geom_point(data = res.BIC[num.min, ], aes(x = n.PC, y = BIC), col = 2, size = 3)
  p <- p + facet_grid(use.mu ~ use.kin)
  p <- p + scale_x_continuous(name = "Number of PCs", breaks = 0:11)
  p <- p + ggtitle(paste0("Model Comparison: trait = ", trait))
  ggsave(filename = paste0(dir.save, "/Model_Comparison_BIC_", trait, ".png"),
         p, width = 8, height = 8)
  
  # merge
  res.BIC.all <- rbind.data.frame(res.BIC.all, res.BIC)
  print(trait)
}

# save result
write.csv(res.BIC.all, file = paste0(dir.save, "/all_BIC_result_", dat, ".csv"), row.names = F)

# save best model
vec.n.PC <- tf.vec.kin <- tf.vec.mu <- rep(NA, length(trait.all))
for ( i in 1:length(trait.all) ) {
  trait <- trait.all[i]
  res.BIC.i <- res.BIC.all[res.BIC.all$trait == trait, ]
  j <- which.min(res.BIC.i$BIC)
  res.BIC.best <- res.BIC.i[j, ]
  vec.n.PC[i] <- res.BIC.best$n.PC
  tf.vec.kin[i] <- res.BIC.best$use.kin == "Use Kinship"
  tf.vec.mu[i] <- res.BIC.best$use.mu == "Use Endospearm Mutant"
}
df.opt.model <- data.frame("trait" = trait.all,
                           "n.PC" = vec.n.PC,
                           "use.kin" = tf.vec.kin,
                           "use.mu" = tf.vec.mu)
df.opt.model$use.kin <- c("no", "yes")[df.opt.model$use.kin + 1]
df.opt.model$use.mu <- c("no", "yes")[df.opt.model$use.mu + 1]
write.csv(df.opt.model, file = paste0(dir.save, "/OptimalModel_BIC_", dat, ".csv"), row.names = F)


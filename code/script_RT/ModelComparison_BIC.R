# source
library(data.table)
library(rrBLUP)
library(ggplot2)

# mkdir
dir.save <- "RESULT/ModelComparison_BIC"
dir.create(dir.save, recursive = T)

# file I/O
file.in.PeerResid <- paste0("RESULT/4.3-Peer_UseOptFact/PeerResult_UseOptFact_BLUE_lmer_residuals.txt")
file.in.kinship <- "RAWDATA/SNPs_and_kinship/Kinship_from_GAPIT.csv"
file.in.pc <- "RAWDATA/SNPs_and_kinship/GAPIT.PCA.csv"
file.in.pheno <- "RAWDATA/Tocochromanols/toco_transformed_blups.csv"
file.in.key <- "RAWDATA/master_key.csv"

# load data
PeerResid <- fread(file.in.PeerResid, data.table = F)
KinMat <- read.csv(file.in.kinship)
Pheno.All <- read.csv(file.in.pheno)
PC <- read.csv(file.in.pc)
key.df <- read.csv(file = file.in.key)

# Peer residual matrix
PeerResidMat <- as.matrix(PeerResid[, -1])
myFun.01 <- function(VEC) { MIN <- min(VEC, na.rm = T); MAX <- max(VEC, na.rm = T); VEC.NEW <- (VEC - MIN) / (MAX - MIN); return(VEC.NEW)}
PeerResidMat.sc <- apply(X = PeerResidMat, 2, FUN = myFun.01)
rownames(PeerResidMat.sc) <- as.character(PeerResid$V1)

# TBD: remove genes without any variation
tf <- apply(is.na(PeerResidMat.sc), 2, sum) == 0
PeerResidMat.sc <- PeerResidMat.sc[, tf]

# # chr & pos
# m <- match(colnames(PeerResidMat.sc), GeneInfo$ID)
# chr <- GeneInfo$chr[m]
# pos <- GeneInfo$start[m]
# snp <- GeneInfo$ID[m] # snp = gene id

# kinship matrix
K <- as.matrix(KinMat[, -1])
rownames(K) <- colnames(K) <- KinMat[, 1]

# PC
PC <- as.matrix(PC[, -1])
rownames(PC) <- rownames(K)

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

# rename a few genotype of PC matrix (same as above)
rownames(PC)[rownames(PC) == "IL731a"] <- "Il731a"
rownames(PC)[rownames(PC) == "IL767b"] <- "Il767b"
rownames(PC)[rownames(PC) == "IL777a"] <- "Il779a"
rownames(PC)[rownames(PC) == "IaEV191"] <- "IaEv191"
rownames(PC)[rownames(PC) == "IL14H"] <- "Il14H"
rownames(PC)[rownames(PC) == "P39a"] <- "P39A"
rownames(PC)[rownames(PC) == "P39a"] <- "P39A"
rownames(PC)[rownames(PC) == "T62s"] <- "T62S"

# intersect of the genotype names among data
common.name <- intersect(rownames(PeerResidMat.sc), rownames(K))
m.peer <- match(common.name, rownames(PeerResidMat.sc))
PeerResidMat.sc.final <- PeerResidMat.sc[m.peer, ]
m.pheno <- match(common.name, Pheno.Sub$Sample.ID)
Pheno.Sub <- Pheno.Sub[m.pheno, ]
m.kinship <- match(common.name, rownames(K))
GRM <- K[m.kinship, m.kinship]
PC <- PC[m.kinship, ]

# check
all.equal(rownames(GRM), rownames(PeerResidMat.sc.final))
all.equal(rownames(GRM), Pheno.Sub$Sample.ID)

# Endosperm.mutation
mutation.vec <- Pheno.Sub$Endosperm.mutation

# remove unnecessary columns & modify column name
Pheno.Sub <- Pheno.Sub[, c(1, 6:ncol(Pheno.Sub))]
colnames(Pheno.Sub)[1] <- "Taxa"
################################################################################

# make a "full" design matrix for fixed effects
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
  b <- (-1) * log(det(H))
  c <- as.numeric((-1) * (1 / ms$Vu) * t(y - X %*% ms$beta) %*% H.inv %*% (y - X %*% ms$beta))
  LogLik.ML <- 0.5 * (a + b + c) # full likelihood based on REML 
  return(LogLik.ML)
}

# LOOP FOR ALL TRAITS
trait.all <- colnames(Pheno.Sub)[2:ncol(Pheno.Sub)]
res.BIC.all <- NULL
for ( i in 1:length(trait.all) ) {
  # target phenotype data
  trait <- trait.all[i]
  pheno <- Pheno.Sub[, trait]
  
  # remove NA
  tf <- !is.na(pheno)
  pheno.i <- pheno[tf]
  GRM.i <- GRM[tf, tf]
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

# # make figure
# df.fig.all$Model <- factor(df.fig.all$Model, levels = c(names(BIC.vec.with.K), names(BIC.vec.without.K)))
# p <- ggplot(df.fig.all, aes(x = Model, y = BIC))
# p <- p + theme_bw()
# p <- p + geom_point(size = 2)
# p <- p + facet_wrap(~ Trait, scales = "free_y")
# p <- p + theme(legend.title = element_text(size = 12,face = "bold", colour = "black"),
#                legend.text = element_text(size = 12),
#                axis.title = element_text(size = 12, face = "bold"),
#                axis.text.x = element_text(size = 8, face = "bold", angle = 90),
#                axis.text.y = element_text(size = 8, face = "bold", angle = 0),
#                legend.position = "none")
# p <- p + geom_vline(xintercept = 11.5, lty = 4, lwd = 1)
# ggsave(filename = "Model_Comparison_BIC_toco.png", p, width = 25, height = 8)


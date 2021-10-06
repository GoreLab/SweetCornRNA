# # Get GAPIT3 from GitHub
# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)

# packages
library(GAPIT3)

# mkdir
dir.save <- "RESULT/6.1-CalcKinship_GAPIT"
dir.create(dir.save)

# params
dat.all <- c("toco", "ion", "carot")
mrk.all <- c("11K", "163K")

# loop for all
for (dat in dat.all) {
	for (mrk in mrk.all) {
		# load files
		f <- paste0("RESULT/4.1-MakeBlueDatasets/PhenoData_", dat, ".csv")
		PhenoDat <- read.csv(f)
		f <- paste0("RAWDATA/master_key.csv")
		KeyDat <- read.csv(f)
		f <- paste0("RAWDATA/SNPs_and_kinship/Ion_", mrk, "_401_v4.hmp.txt")
		myG <- read.table(f, head = FALSE)
		
		# ---------------------------------------------------------------------------- #
		# Convert ID to genotype (= ion) data
		m <- match(PhenoDat$Sample.ID, KeyDat$accession_name)
		id.of.phenotyped.ones <- KeyDat$accession_name_ion[m]
		
		# convert genotype data ID (to match with phenotype data)
		id.all <- unlist(myG[1, 12:ncol(myG)], use.names = F)
		myfun <- function(x){ strsplit(x, "\\.")[[1]][1] }
		id.all.2 <- sapply(id.all, myfun, USE.NAMES = F)
		id.all.2[id.all.2 == "X2"] <- "2"
		id.all.2[id.all.2 == "X257A9"] <- "257A9"
		id.all.2[id.all.2 == "X304A"] <- "304A"
		id.all.2[id.all.2 == "X34f"] <- "34f"
		id.all.2[id.all.2 == "X471_U6"] <- "471_U6"
		id.all.2[id.all.2 == "X675A"] <- "675A"
		id.all.2[id.all.2 == "X83610b"] <- "83610b"
		id.all.2[id.all.2 == "X83612b"] <- "83612b"
		
		# match
		tf <- id.all.2 %in% id.of.phenotyped.ones
		myG.sub <- myG[, c(rep(TRUE, 11), tf)]
		
		# check
		if ( sum(tf) != nrow(PhenoDat) ) {
			print("Something Wrong!!")
		}
		
		# pseudo phenotype data
		myY <- data.frame("Taxa" = unlist(myG.sub[1, 12:ncol(myG.sub)], use.names = F),
											"Pheno" = 1:(ncol(myG.sub)-12+1))
		
		# move to tmp.dir
		tmp.dir <- paste0(dir.save, "/LogFile_", dat, "_", mrk)
		dir.create(tmp.dir, recursive = T)
		curr.dir <- getwd()
		full.path <- paste0(curr.dir, "/", tmp.dir)
		setwd(full.path)
		
		# run GAPIT
		myGAPIT <- GAPIT(Y = myY, G = myG.sub, SNP.MAF = 0.05)
		warnings() # show warnings (just in case)
		
		# move to original dir
		setwd(curr.dir)
		
		# load kinship matrix of GAPIT
		f <- paste0(tmp.dir, "/GAPIT.Kin.VanRaden.csv")
		K <- read.csv(f, header = F)
		KinMat <- as.matrix(K[, -1])
		
		# convert IDs
		myfun <- function(x){ strsplit(x, "\\.")[[1]][1] }
		id.vec <- sapply(K[, 1], myfun, USE.NAMES = F)
		id.vec[id.vec == "X2"] <- "2"
		id.vec[id.vec == "X257A9"] <- "257A9"
		id.vec[id.vec == "X304A"] <- "304A"
		id.vec[id.vec == "X34f"] <- "34f"
		id.vec[id.vec == "X471_U6"] <- "471_U6"
		id.vec[id.vec == "X675A"] <- "675A"
		id.vec[id.vec == "X83610b"] <- "83610b"
		id.vec[id.vec == "X83612b"] <- "83612b"
		rownames(KinMat) <- colnames(KinMat) <- id.vec
		
		# convert IDs
		m <- match(rownames(KinMat), KeyDat$accession_name_ion)
		id.vec.new <- KeyDat$accession_name[m]
		rownames(KinMat) <- colnames(KinMat) <- id.vec.new
		
		# sort 
		setequal(rownames(KinMat), PhenoDat$Sample.ID) # check -> OK!
		m <- match(PhenoDat$Sample.ID, rownames(KinMat))
		KinMat.sort <- KinMat[m, m]
		all(rownames(KinMat.sort) == PhenoDat$Sample.ID) # check -> OK!
		
		# save kinship matrix (or GRM)
		f <- paste0(dir.save, "/KinMat_", dat, "_", mrk, ".csv")
		write.csv(KinMat.sort, file = f)
	}
}


# package
library(data.table)
library(gdata)

# params
args <- commandArgs(trailingOnly = T)
ref <- args[1] # "B73" or "PH207"

# mkdir
dir.save <- "RESULT_NoHarvDate/4.1-MakeBlueDatasets"
dir.create(dir.save)

# ---------------------------------------------------------------------------- #
# input files
file.blue <- paste0("RESULT_NoHarvDate/3.1-BLUE/expression_BLUE_RmErr_", ref, ".csv")
file.key <- "RAWDATA/master_key.csv"

# load BLUE data
BLUE.dat <- fread(file.blue, data.table = F)
BLUE.mat.all <- t(as.matrix(BLUE.dat[, -1]))
rownames(BLUE.mat.all)[rownames(BLUE.mat.all) == "X2"] <- "2"
rownames(BLUE.mat.all)[rownames(BLUE.mat.all) == "X257A9"] <- "257A9"
rownames(BLUE.mat.all)[rownames(BLUE.mat.all) == "X304A"] <- "304A"
rownames(BLUE.mat.all)[rownames(BLUE.mat.all) == "X34f"] <- "34f"
rownames(BLUE.mat.all)[rownames(BLUE.mat.all) == "X471_U6"] <- "471_U6"
rownames(BLUE.mat.all)[rownames(BLUE.mat.all) == "X675A"] <- "675A"
rownames(BLUE.mat.all)[rownames(BLUE.mat.all) == "X83610b"] <- "83610b"
rownames(BLUE.mat.all)[rownames(BLUE.mat.all) == "X83612b"] <- "83612b"
colnames(BLUE.mat.all) <- BLUE.dat$GeneID

# load genotype data
GenDat <- fread("RAWDATA/SNPs_and_kinship/Ion_11K_401_v4.hmp.txt", data.table = F)
GenDat[1:5, 1:20]

# load ion data
IonDat <- read.xls("RAWDATA/Ionomics/SuppTableS3.BLUPs_Oct20.xlsx")
head(IonDat)
dim(IonDat)

# load toco data
TocoDat <- read.csv("RAWDATA/Tocochromanols/toco_transformed_blups.csv")
head(TocoDat)
dim(TocoDat)

# load carotenoid data
CaroDat.Raw <- read.csv("RAWDATA/Carotenoids/carotenoid_transformed_BLUPs-converted.csv",
												skip = 1)
CaroDat <- CaroDat.Raw[1:(nrow(CaroDat.Raw)-1), ] # last row is not needed (comment row in excel file)
head(CaroDat)
dim(CaroDat)

# load key file
key.dat <- read.csv(file.key)
head(key.dat)


# ---------------------------------------------------------------------------- #
# remove check lines
tf <- rownames(BLUE.mat.all) %in% c("CHECK1", "CHECK2", "CHECK3", "CHECK4")
BLUE.mat <- BLUE.mat.all[!tf, ]

# ---------------------------------------------------------------------------- #
# check overlap: key file can cover all samples
setdiff(rownames(BLUE.mat), key.dat$accession_name) # ok
setdiff(key.dat$accession_name, rownames(BLUE.mat)) # ok

# subset key file
m <- match(rownames(BLUE.mat), key.dat$accession_name)
key.dat <- key.dat[m, ] # match samples
all(key.dat$accession_name == rownames(BLUE.mat)) # check -> ok
head(key.dat)


# ---------------------------------------------------------------------------- #
# id in genotype data 
id.gen <- colnames(GenDat)[12:ncol(GenDat)]
myfun <- function(x){strsplit(x, "\\.")[[1]][1]}
id.gen <- sapply(id.gen, myfun, USE.NAMES = F)
id.gen[id.gen == "X2"] <- "2"
id.gen[id.gen == "X257A9"] <- "257A9"
id.gen[id.gen == "X304A"] <- "304A"
id.gen[id.gen == "X34f"] <- "34f"
id.gen[id.gen == "X675A"] <- "675A"
id.gen[id.gen == "X83610b"] <- "83610b"
id.gen[id.gen == "X83612b"] <- "83612b"
id.gen[id.gen == "X471_U6"] <- "471_U6"
id.gen[id.gen == "X81_1"] <- "81_1"
setequal(id.gen, IonDat$Sample.ID) # ok!


# ---------------------------------------------------------------------------- #
# IDs: from expression data to each of the three
id.expr.with.check <- rownames(BLUE.mat)
id.expr <- setdiff(id.expr.with.check, c("CHECK1", "CHECK2", "CHECK3", "CHECK4"))
id.expr.ion <- key.dat$accession_name_ion[match(id.expr, key.dat$accession_name)]
id.expr.toco <- key.dat$accession_name_toco[match(id.expr, key.dat$accession_name)]
id.expr.caro <- key.dat$accession_name_carot[match(id.expr, key.dat$accession_name)]

# IDs of phenotype data
id.ion <- IonDat$Sample.ID
id.toco <- TocoDat$Sample.ID
id.caro <- CaroDat$Sample.ID

# count numbers
length(id.expr) # 355 accessions (OK)
length(id.expr.ion[!is.na(id.expr.ion)]) # 351 accessions ovrlapped between ion & expr (OK)
length(id.expr.toco[!is.na(id.expr.toco)]) # 354 accessions ovrlapped between toco & expr (OK)
length(id.expr.caro[!is.na(id.expr.caro)]) # 285 accessions ovrlapped between caro & expr (OK)

# count numbers
length(id.ion) # 401 (OK)
length(id.toco) # 384 (OK)
length(id.caro) # 308 (OK)


# ---------------------------------------------------------------------------- #
# samples for ion data
tf <- !is.na(key.dat$accession_name_ion); table(tf) # 4 genotypes not included
id.expr.retain.ion <- key.dat$accession_name[tf]
length(id.expr.retain.ion) # number of genotypes to be retained = 351

# samples for toco data
tf01 <- !is.na(key.dat$accession_name_toco); table(tf01) # 1 genotype not included
tf02 <- !is.na(key.dat$accession_name_ion); table(tf02) # 4 genotypes not included
tf <- tf01 & tf02
id.expr.retain.toco <- key.dat$accession_name[tf]
length(id.expr.retain.toco) # number of genotypes to be retained = 351

# samples for carot data
tf01 <- !is.na(key.dat$accession_name_carot); table(tf01) # 70 genotype not included
tf02 <- !is.na(key.dat$accession_name_ion); table(tf02) # 4 genotypes not included
tf <- tf01 & tf02
id.expr.retain.carot <- key.dat$accession_name[tf]
length(id.expr.retain.carot) # number of genotypes to be retained = 283


# ---------------------------------------------------------------------------- #
# make data of ion
m.ion <- match(id.expr.retain.ion, key.dat$accession_name)
id.ion.retain <- key.dat$accession_name_ion[m.ion]
m <- match(id.ion.retain, IonDat$Sample.ID)
IonDat.m <- IonDat[m, ]
IonDat.m$Sample.ID.expr <- id.expr.retain.ion
IonDat.save <- IonDat.m[, c("Sample.ID.expr", "Endosperm.mutation",
														"Boron", "Cadmium", "Calcium",
														"Copper", "Iron", "Magnesium",
														"Manganese", "Molybdenum", "Nickel",
														"Phosphorus", "Potassium", "Rubidium",
														"Strontium", "Sulfur", "Zinc")]
colnames(IonDat.save)[1] <- "Sample.ID"
if ( ref == "B73" ) {
	fwrite(IonDat.save, file = paste0(dir.save, "/PhenoData_ion.csv"))
} # same for both B73 and PH207
BLUE.ion.save.all <- data.frame("Sample.ID" = rownames(BLUE.mat),
																BLUE.mat)
m <- match(IonDat.save$Sample.ID, BLUE.ion.save.all$Sample.ID)
BLUE.ion.save <- BLUE.ion.save.all[m, ]
fwrite(BLUE.ion.save, file = paste0(dir.save, "/BLUE_matrix_for_ion_", ref, ".csv"))


# ---------------------------------------------------------------------------- #
# make data of toco
m.toco <- match(id.expr.retain.toco, key.dat$accession_name)
id.toco.retain <- key.dat$accession_name_toco[m.toco]
m <- match(id.toco.retain, TocoDat$Sample.ID)
TocoDat.m <- TocoDat[m, ]
TocoDat.m$Sample.ID.expr <- id.expr.retain.toco
TocoDat.save <- TocoDat.m[, c("Sample.ID.expr", "Endosperm.mutation",
															"aT", "aT3", "dT", "dT3", "gT", "gT3",
															"Total.T", "Total.T3", "Total.T3_plus_T",
															"ratio_aT_gT", "ratio_dT_aT", "ratio_dt_gT",
															"dT_over_sum_gT_aT", "gt_over_sum_gT_aT",
															"ratio_aT3_gT3", "ratio_dT3_aT3", "ratio_dT3_gT3",
															"dT3_over_sum_gT3_aT3", "gT3_over_sum_gT3_aT3",
															"ratio_TotalT_TotalT3")]
colnames(TocoDat.save)[1] <- "Sample.ID"
if ( ref == "B73" ) {
	fwrite(TocoDat.save, file = paste0(dir.save, "/PhenoData_toco.csv"))
} # same for both B73 and PH207
BLUE.toco.save.all <- data.frame("Sample.ID" = rownames(BLUE.mat),
																 BLUE.mat)
m <- match(TocoDat.save$Sample.ID, BLUE.toco.save.all$Sample.ID)
BLUE.toco.save <- BLUE.toco.save.all[m, ]
fwrite(BLUE.toco.save, file = paste0(dir.save, "/BLUE_matrix_for_toco_", ref, ".csv"))


# ---------------------------------------------------------------------------- #
# make data of carot
m.carot <- match(id.expr.retain.carot, key.dat$accession_name)
id.carot.retain <- key.dat$accession_name_carot[m.carot]
m <- match(id.carot.retain, CaroDat$Sample.ID)
CaroDat.m <- CaroDat[m, ]
CaroDat.m$Sample.ID.expr <- id.expr.retain.carot
colnames(CaroDat.m) <- c("Sample.ID", "Additional.ID.name", "GBS.ID",
												 "Endosperm.mutation", "Antheraxanthin", 
												 "beta.Carotene", "beta.Cryptoxanthin", "Lutein",
												 "Violaxanthin", "Zeaxanthin", "Zeinoxanthin", 
												 "Other.carotenes", "alpha.Xanthophylls",
												 "beta.Xanthophylls", "Total.xanthophylls",
												 "Total.carotenes", "Total.carotenoids",
												 "beta.Carotene_over_beta.cryptoxanthin",
												 "beta.Carotene_over_sum_beta.cryptoxanthin_and_zeaxanthin",
												 "beta.Cryptoxanthin_over_zeaxanthin",
												 "Zeinoxanthin_over_lutein", 
												 "beta.Xanthophylls_over_alpha.xanthophylls",
												 "Total.carotenes_over_total.xanthophylls",
												 "X", "Sample.ID.expr")
CaroDat.save <- CaroDat.m[, c("Sample.ID.expr",
															"Endosperm.mutation", "Antheraxanthin", 
															"beta.Carotene", "beta.Cryptoxanthin", "Lutein",
															"Violaxanthin", "Zeaxanthin", "Zeinoxanthin", 
															"Other.carotenes", "alpha.Xanthophylls",
															"beta.Xanthophylls", "Total.xanthophylls",
															"Total.carotenes", "Total.carotenoids",
															"beta.Carotene_over_beta.cryptoxanthin",
															"beta.Carotene_over_sum_beta.cryptoxanthin_and_zeaxanthin",
															"beta.Cryptoxanthin_over_zeaxanthin",
															"Zeinoxanthin_over_lutein", 
															"beta.Xanthophylls_over_alpha.xanthophylls",
															"Total.carotenes_over_total.xanthophylls")]
colnames(CaroDat.save)[1] <- "Sample.ID"
if ( ref == "B73" ) {
	fwrite(CaroDat.save, file = paste0(dir.save, "/PhenoData_carot.csv"))
} # same for both B73 and PH207
BLUE.caro.save.all <- data.frame("Sample.ID" = rownames(BLUE.mat),
																 BLUE.mat)
m <- match(CaroDat.save$Sample.ID, BLUE.caro.save.all$Sample.ID)
BLUE.caro.save <- BLUE.caro.save.all[m, ]
fwrite(BLUE.caro.save, file = paste0(dir.save, "/BLUE_matrix_for_carot_", ref, ".csv"))


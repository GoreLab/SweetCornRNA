# Combined test

# library
library(data.table)
library(metap)
library(dplyr)

# mkdir
dir.create("RESULT/5.2-CombinedTest_with_all")

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"
trait <- args[3] # trait
chr.twas <- as.numeric(args[4]) # chr

# ###############################################
# ver <- "v1" # "v1" or "v1.1"
# ref <- "B73" # "B73" or "PH207"
# trait <- "d.T"
# chr.twas <- 1
# ###############################################

# filename I/O
if ( ref == "B73" ) { file.gff <- paste0("RAWDATA/Annotation/Zea_mays.B73_RefGen_v4.59_anno_merged.csv") }
if ( ref == "PH207" ) { file.gff <- paste0("RAWDATA/Annotation/Zm-PH207-REFERENCE_NS-UIUC_UMN-1.0.v03eksc_anno_merged.csv") }
file.gwas.all <- "RAWDATA/GWAS/vita_v6_locked_hmp_LD0.1_K_GAPIT_MLM_all_top0.1_24_traits.csv"
file.pheno.all <- "RAWDATA/Phenotype/ames_vita_blue_trans_v6version_locked_24_traits.txt"

# load phenotype file (to get names of traits)
Pheno.All <- read.delim(file.pheno.all)
trait.all <- colnames(Pheno.All)[2:ncol(Pheno.All)]

# load GFF & GWAS result
GeneInfo <- fread(file.gff, data.table = F)
gwas.res.all <- fread(file.gwas.all, data.table = F)

# print
print(paste0("Run FC test for ", trait))

# filename I/O
file.twas <- paste0("RESULT/4.1-TWAS/TwasResult_", ver, "_", ref, "_", trait, ".csv")
file.save <- paste0("RESULT/5.2-CombinedTest_with_all/Result_FCtest_", ver, "_", ref, "_", trait, "_chr", chr.twas, ".csv")

# load TWAS result
twas.res <- read.csv(file.twas, stringsAsFactors = F)
twas.res <- twas.res[twas.res$Chr %in% 1:10, ] # retain chr = 1:10

# GWAS result of the target trait
gwas.res <- gwas.res.all[gwas.res.all$trait == trait, ]
gc();gc()

# ==== NEW ==== subset of GWAS result
o <- order(gwas.res$P.value)
gwas.res.sub <- gwas.res[o[1:round(length(o) * 0.01, 1)], ]
rm(gwas.res); gc();gc()

# 1. remake twas data (assign P = 1 for non-avairable genes)
dat.twas.all <- data.frame("ID" = GeneInfo$ID,
                           "CHR" = GeneInfo$chr,
                           "START" = GeneInfo$start,
                           "END" = GeneInfo$end,
                           "GROUP.NEW" = GeneInfo$group,
                           "START.NEW" = GeneInfo$new.start,
                           "END.NEW" = GeneInfo$new.end,
                           "pvalue" = 1,
                           stringsAsFactors = F)
m <- match(twas.res$Gene, GeneInfo$ID)
dat.twas.all$pvalue[m] <- 10 ^ -(twas.res$neg.log.P) # p-values in TWAS

# ==== NEW ==== TWAS for a chromosome
dat.twas.all <- dat.twas.all[dat.twas.all$CHR == chr.twas, ]
dim(dat.twas.all)

# ==== NEW ==== 2. MERGE TWAS + GWAS for all pairs
dat.twas.to.merge <- data.frame("gene" = dat.twas.all$ID,
                                "pvalue.twas" = dat.twas.all$pvalue)
dat.gwas.to.merge <- data.frame("snp" = gwas.res.sub$SNP,
                                "pvalue.gwas" = gwas.res.sub$P.value)
dat.merged <- merge(dat.twas.to.merge, dat.gwas.to.merge)

# 3. FC.TEST
pval.mat <- cbind(dat.merged$pvalue.twas, dat.merged$pvalue.gwas)
combined.p <- apply(pval.mat, 1, FUN = function(x){sumlog(x)$p})
dat.merged$FCT.pval <- combined.p

# save result
fwrite(x = dat.merged, file = file.save)


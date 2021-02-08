# source
library(data.table)
library(lme4)
library(ggplot2)
library(ggExtra)
library(e1071)

# --------------------------------------------------------------------------- #
# ----- 1. Load data & take average
# --------------------------------------------------------------------------- #
# file I/O
dir.in.rlog <- "RAWDATA/Seetcorn_TagSeq"
file.in.rlog <- "htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.txt"
dir.in.key <- "RAWDATA"
file.in.key <- "master_key.csv"
dir.in.expinfo <- "RAWDATA/Metadata"
file.in.expinfo <- "rnaseq_trial_2019_upload.csv"

# mkdir
dir.save <- "RESULT/99-FurtherDiagnosis_Peer"
dir.create(dir.save, recursive = T)

# objects
i <- 1
name.tech.ctrl <- c("Control_NA", "fill_c1", "fill_c2", "fill_c3", "fill_c4")

# load rlog data
rlog.df <- fread(file = paste0(dir.in.rlog, "/", file.in.rlog), data.table = F) # a few min

# rlog data matrix
gene.info.df <- rlog.df[, 1:5]
rlog.mat <- as.matrix(rlog.df[, 6:ncol(rlog.df)])
rownames(rlog.mat) <- gene.info.df$gene_id

# add sample name
h <- readLines(con = paste0(dir.in.rlog, "/", file.in.rlog), n = 1) # get header
name.header <- strsplit(h, "\t")[[1]]
colnames(rlog.mat) <- name.header[6:ncol(rlog.df)]

# load master data
key.df <- read.csv(file = paste0(dir.in.key, "/", file.in.key))

# check: list of names are identical? -- yes!
setequal(key.df$rlog_ids, colnames(rlog.mat))

# load experiment info data
expinfo.df <- read.csv(file = paste0(dir.in.expinfo, "/", file.in.expinfo))

# check: info is available for all samples?
setdiff(key.df$sample_name_metadata, expinfo.df$plot_name) # OK -- "Pos" is positive control

# add rlog value of the i-th gene to the key file
df.model <- data.frame("rlog.id" = colnames(rlog.mat),
                       "Expr.rlog" = rlog.mat[i, ],
                       row.names = NULL)
df.model <- merge.data.frame(x = df.model, y = key.df, by.x = "rlog.id", by.y = "rlog_ids")

# remove technical control in RNAseq
df.model <- df.model[!(df.model$accession_name %in% name.tech.ctrl), ]

# merge metadata
length(expinfo.df$plot_name) == length(unique(expinfo.df$plot_name)) # uniqueness: OK
length(df.model$sample_name_metadata) == length(unique(df.model$sample_name_metadata)) # uniqueness: OK
setdiff(df.model$sample_name_metadata, expinfo.df$plot_name) # check inclusion: OK
df.model <- merge.data.frame(x = df.model, y = expinfo.df,
                             by.x = "sample_name_metadata", by.y = "plot_name",
                             all.x = TRUE, all.y = FALSE)

# rename columns for simplicity
colnames(df.model)[colnames(df.model) == "accession_name.y"] <- "accession_name_meta"
colnames(df.model)[colnames(df.model) == "accession_name.x"] <- "accession_name"

# Harvest.Date
days.harvest <- as.integer(as.Date(df.model$tissue_harvest_date) - min(as.Date(df.model$tissue_harvest_date)))
days.harvest.sc <- scale(days.harvest)

# days.rna.extract should be two-level factor
rna.extract <- rep("spring", nrow(df.model))
rna.extract["2020-05-01" < as.Date(df.model$nucleic_acid_extraction_date)] <- "summer"
rna.extract <- as.factor(rna.extract)

# final data
df.fit <- data.frame("y" = df.model$Expr.rlog,
                     "Genotype" = df.model$accession_name,
                     "Range" = as.factor(df.model$range_number),
                     "Block" = as.factor(df.model$block_number),
                     "Column" = as.factor(df.model$col_number),
                     "Row" = as.factor(df.model$row_number),
                     "Plate" = as.factor(df.model$plate_number),
                     "Harvest.Date" = days.harvest.sc,
                     "RnaExtract" = as.factor(rna.extract),
                     "ID" = df.model$rlog.id)

# take average
rlog.mat.blue.samples <- rlog.mat[, df.fit$ID]
gen <- levels(as.factor(df.fit$Genotype))
avg.mat <- matrix(NA, nr = nrow(rlog.mat.blue.samples), nc = length(gen))
colnames(avg.mat) <- gen
rownames(avg.mat) <- rownames(rlog.mat.blue.samples)
for ( i in 1:nrow(rlog.mat.blue.samples) ) {
   y <- rlog.mat.blue.samples[i, ]
   avg <- tapply(y, df.fit$Genotype, mean)
   m <- match(names(avg), colnames(avg.mat))
   avg.mat[i, m] <- avg
}


# --------------------------------------------------------------------------- #
# ----- 2. Load PEER result
# --------------------------------------------------------------------------- #
# PEER result
peer.raw <- fread("RESULT/4.3-Peer_UseOptFact/PeerResult_UseOptFact_BLUE_lmer_residuals.txt", data.table = F)

# df -> matrix
sample.name <- peer.raw$V1
peer.mat <- t(as.matrix(peer.raw[, -1]))
colnames(peer.mat) <- sample.name

# remove checks
setdiff(colnames(peer.mat), colnames(avg.mat)) # OK
setdiff(colnames(avg.mat), colnames(peer.mat)) # OK
avg.mat <- avg.mat[, colnames(peer.mat)]

# remove genes without BLUE
setdiff(rownames(peer.mat), rownames(avg.mat)) # OK
setdiff(rownames(avg.mat), rownames(peer.mat)) # 1160 genes
avg.mat <- avg.mat[rownames(peer.mat), ] 


# --------------------------------------------------------------------------- #
# ----- 3. Compare
# --------------------------------------------------------------------------- #
# correlation
cor.vec <- rep(NA, times = nrow(avg.mat))
names(cor.vec) <- rownames(avg.mat)
for ( i in 1:nrow(avg.mat) ) {
   cor.vec[i] <- cor(avg.mat[i, ], peer.mat[i, ])
}

# histogram
hist(cor.vec, breaks = 30, xlim = c(0, 1),
     xlab = "Pearson's correlation: cor(PEER, AVG)",
     main = paste0("Histogram of ", length(cor.vec), " genes"))

# 
df <- data.frame("Raw" = avg.mat["Zm00001d035395", ],
                 "PEER.resid" = peer.mat["Zm00001d035395", ])
p <- ggplot(df, aes(Raw, PEER.resid)) + geom_point() + theme_classic()
p <- p + ggtitle("Zm00001d035395")
p <- ggMarginal(p, type = "histogram")
ggsave(paste0(dir.save, "/Zm00001d035395.png"), p, width = 5, height = 5)

# 
df <- data.frame("Raw" = avg.mat["Zm00001d019155", ],
                 "PEER.resid" = peer.mat["Zm00001d019155", ])
p <- ggplot(df, aes(Raw, PEER.resid)) + geom_point() + theme_classic()
p <- p + ggtitle("Zm00001d019155")
p <- ggMarginal(p, type = "histogram")
ggsave(paste0(dir.save, "/Zm00001d019155.png"), p, width = 5, height = 5)


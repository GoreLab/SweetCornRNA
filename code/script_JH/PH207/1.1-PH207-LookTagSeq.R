# look at the tagseq data

# library
library(ggplot2)
library(gplots)
library(e1071)

# file I/O
dir.in <- "data" # /Sweetcorn_TagSeq_PH207"
file.in <- "htseq_count_matrix_sweetcorn_PH207_RLOG_all_info_v1.txt"

# load data
rlog.df <- read.delim(file = paste0(dir.in, "/", file.in))

# data
gene.info.df <- rlog.df[, 1:6]
rlog.mat <- as.matrix(rlog.df[, 7:ncol(rlog.df)])
rownames(rlog.mat) <- gene.info.df$gene_id

# data size
dim(rlog.mat) # 20065 genes * 433 samples (samples include checks)

# get sample info
sample.name.vec <- colnames(rlog.mat)
plate <- sapply(X = sample.name.vec, FUN = function(x){strsplit(x, "_")[[1]][3]}, USE.NAMES = F)
well <- sapply(X = sample.name.vec, FUN = function(x){strsplit(x, "_")[[1]][4]}, USE.NAMES = F)
id <- sapply(X = sample.name.vec, FUN = function(x){strsplit(x, "_")[[1]][5]}, USE.NAMES = F)
genotype <- sapply(X = sample.name.vec, FUN = function(x){strsplit(x, "_")[[1]][6]}, USE.NAMES = F)
sample.info.df <- data.frame("Sample.Name" = sample.name.vec,
                             "Plate" = plate,
                             "Well" = well,
                             "ID" = id,
                             "Genotype" = genotype)

# number of samples in each plate
table(sample.info.df$Plate) # 72 ~ 91 samples
sum(sample.info.df$ID == "Pos") # 9 positive controls
sum(sample.info.df$Genotype == "Control") # 9 positive controls (double check)
sum(sample.info.df$Genotype == "CHECK1") # 5 check1
sum(sample.info.df$Genotype == "CHECK2") # 16 check2
sum(sample.info.df$Genotype == "CHECK3") # 14 check3
sum(sample.info.df$Genotype == "CHECK4") # 5 check4

# missing value?
sum(is.na(rlog.mat)) # 20741 NA
rlog.vec.zero.count <- apply(rlog.mat == 0, 2, sum)
sum(is.na(rlog.vec.zero.count)) # all 20741 NA are for one sample
which(is.na(rlog.vec.zero.count)) # this sample is the one... the last sample.
  # any error in upload/download process? (truncation happened?)

# make a histogram
rlog.vec.non.zero.ratio <- 1 - rlog.vec.zero.count / nrow(rlog.mat)
hist(rlog.vec.non.zero.ratio,
     xlab = "Proportion of non-zero rlog2 value",
     main = "Histogram of 433 samples")

# make a boxplot for plate effects
rlog.vec.zero.count <- apply(rlog.mat == 0, 2, sum)
df.fig <- data.frame("Number.of.zero.count" = rlog.vec.zero.count,
                     "Plate" = sample.info.df$Plate)
df.fig <- df.fig[sample.info.df$ID != "Pos", ] # remove positive control
p <- ggplot(df.fig, aes(x = Plate, y = Number.of.zero.count))
p <- p + geom_violin()
p <- p + geom_boxplot(width = 0.1)
p <- p + ggtitle("# of zero count for different trait")
p

# correlation among samples
# tf.NA <- colnames(rlog.mat) %in% "SC_RNA_05_B10_19A0135_We09425"
# rlog.mat.rm.NA <- rlog.mat[, !tf.NA]
CorMat.Sample <- cor(rlog.mat)#.rm.NA)
vec.min.cor <- apply(CorMat.Sample, 1, min)
vec.lowest.cor <- apply(CorMat.Sample, 1, which.min)

# make a heatmap
heatmap.2(CorMat.Sample, trace = 'none')

# lowest correlation
table(vec.lowest.cor) # PH207 found 4 samples: 142, 165, 385, 407 ## NOTE 407 and 165 were also id'd in B73
colnames(rlog.mat)[c(142, 165, 385, 407)] # these are the three samples

# plots
plot(x = rlog.mat[, 407], y = rlog.mat[, 1], pch = 20,
     xlab = colnames(rlog.mat)[407], ylab = colnames(rlog.mat)[1])
plot(x = rlog.mat[, 385], y = rlog.mat[, 1], pch = 20,
     xlab = colnames(rlog.mat)[385], ylab = colnames(rlog.mat)[1])
plot(x = rlog.mat[, 165], y = rlog.mat[, 1], pch = 20,
     xlab = colnames(rlog.mat)[165], ylab = colnames(rlog.mat)[1])
plot(x = rlog.mat[, 142], y = rlog.mat[, 1], pch = 20,
     xlab = colnames(rlog.mat)[142], ylab = colnames(rlog.mat)[1])

which.max(rlog.mat[, 407])
hist(rlog.mat["Zm00008a023627", ], breaks = 100,
     xlab = "rlog2 of Zm00008a023627",
     main = "Histogram of 20065 genes")

which.max(rlog.mat[, 165])
hist(rlog.mat["Zm00008a022414", ], breaks = 100,
     xlab = "rlog2 of Zm00008a022414",
     main = "Histogram of 20065 genes")

which.max(rlog.mat[, 142])
hist(rlog.mat["Zm00008a027192", ], breaks = 100,
     xlab = "rlog2 of Zm00008a027192",
     main = "Histogram of 20065 genes")

# histogram of minimun cor
hist(vec.min.cor, xlab = "Minimum correlation",
     main = "Histogram of 433 samples")

# number of zero count for each gene
# tf.NA <- colnames(rlog.mat) %in% "SC_RNA_05_B10_19A0135_We09425"
# rlog.mat.rm.NA <- rlog.mat[, !tf.NA]
non.zero.count.vec.genes <- ncol(rlog.mat) - apply(rlog.mat == 0, 1, sum)
zero.count.vec.genes <- ncol(rlog.mat) - non.zero.count.vec.genes
zero.prop.vec.genes <- zero.count.vec.genes / ncol(rlog.mat)
hist(zero.prop.vec.genes, breaks = 20,
     xlab = "Proportion of zero rlog2 value",
     main = "Histogram of 20065 genes")

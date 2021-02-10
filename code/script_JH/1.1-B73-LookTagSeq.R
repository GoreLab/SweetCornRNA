# look at the tagseq data

# library
library(ggplot2)
library(gplots)
library(e1071)

# file I/O
dir.in <- "data" # /Sweetcorn_TagSeq_PH207"
file.in <- "htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.txt"

# load data
rlog.df <- read.delim(file = paste0(dir.in, "/", file.in))

# data
gene.info.df <- rlog.df[, 1:5]
rlog.mat <- as.matrix(rlog.df[, 6:ncol(rlog.df)])
rownames(rlog.mat) <- gene.info.df$gene_id

# data size
dim(rlog.mat) # 20741 genes * 433 samples (samples include checks)

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


# number of zero count for each gene
non.zero.count.vec.genes <- ncol(rlog.mat) - apply(rlog.mat == 0, 1, sum)
zero.count.vec.genes <- ncol(rlog.mat) - non.zero.count.vec.genes
zero.prop.vec.genes <- zero.count.vec.genes / ncol(rlog.mat)
hist(zero.prop.vec.genes, breaks = 20,
     xlab = "Proportion of zero rlog2 value",
     main = "B73 Histogram of 20741 genes")
abline(v = c(0.9, 0.5), col = c("red", "blue"), lty=c(2,2))

# If we apply a filter removing genes with greater than 90% of samples with a value of zero:
sum(zero.prop.vec.genes > .9) # we remove 1265 genes == 6.3%
1160/20741
# filter of 0.5
sum(zero.prop.vec.genes > .5)
1859/20741


# make a histogram
rlog.vec.non.zero.ratio <- 1 - rlog.vec.zero.count / nrow(rlog.mat)
hist(rlog.vec.non.zero.ratio,
     xlab = "Proportion of non-zero rlog2 value",
     main = "Histogram of 433 samples", col = "grey")
min(rlog.vec.non.zero.ratio)

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
CorMat.Sample <- cor(rlog.mat)
vec.min.cor <- apply(CorMat.Sample, 1, min)
vec.lowest.cor <- apply(CorMat.Sample, 1, which.min)

# make a heatmap
heatmap.2(CorMat.Sample, trace = 'none')

# lowest correlation
table(vec.lowest.cor)
colnames(rlog.mat)[c(407, 314, 165)] # these are the three samples

# plots
plot(x = rlog.mat[, 407], y = rlog.mat[, 1], pch = 20,
     xlab = colnames(rlog.mat)[407], ylab = colnames(rlog.mat)[1])
plot(x = rlog.mat[, 314], y = rlog.mat[, 1], pch = 20,
     xlab = colnames(rlog.mat)[314], ylab = colnames(rlog.mat)[1])
plot(x = rlog.mat[, 165], y = rlog.mat[, 1], pch = 20,
     xlab = colnames(rlog.mat)[165], ylab = colnames(rlog.mat)[1])
which.max(rlog.mat[, 407])
hist(rlog.mat["Zm00001d035395", ], breaks = 100,
     xlab = "rlog2 of Zm00001d035395",
     main = "Histogram of 20741 genes")
which.max(rlog.mat[,314])
hist(rlog.mat["Zm00001d019155", ], breaks = 100,
     xlab = "rlog2 of Zm00001d019155",
     main = "Histogram of 20741 genes")


# histogram of minimun cor
hist(vec.min.cor, xlab = "Minimum correlation",
     main = "Histogram of 433 samples")

# number of zero count for each gene
non.zero.count.vec.genes <- ncol(rlog.mat) - apply(rlog.mat == 0, 1, sum)
zero.count.vec.genes <- ncol(rlog.mat) - non.zero.count.vec.genes
zero.prop.vec.genes <- zero.count.vec.genes / ncol(rlog.mat)
hist(zero.prop.vec.genes, breaks = 20,
     xlab = "Proportion of zero rlog2 value",
     main = "Histogram of 20741 genes")



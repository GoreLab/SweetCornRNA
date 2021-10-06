# look at the tagseq data

# library
library(ggplot2)
library(gplots)
library(e1071)
library(data.table)

# params
ref <- "Ia453"

# mkdir
dir.save <- "RESULT/1.1-LookTagSeq_MakeData"
dir.create(dir.save, recursive = TRUE)

# load data
dir.in <- "RAWDATA/Seetcorn_TagSeq"
file.in <- "htseq_count_matrix_sweetcorn_Ia453_RLOG_all_info_v1.txt"
rlog.df <- fread(file = paste0(dir.in, "/", file.in), data.table = F)
gene.info.df <- rlog.df[, 1:6]
rlog.mat <- as.matrix(rlog.df[, 7:ncol(rlog.df)])

# data
rownames(rlog.mat) <- gene.info.df$gene_id

# data size
dim(rlog.mat) # 20216 genes * 433 samples (samples include checks)

# missing value?
sum(is.na(rlog.mat)) # no missing value!

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

# remove controls
name.tech.ctrl <- c("Control", "fill")
tf <- sample.info.df$Genotype %in% name.tech.ctrl
rlog.mat <- rlog.mat[, !tf]
sample.info.df <- sample.info.df[!tf, ]

# count numbers
tf <- sample.info.df$Genotype %in% c("CHECK1", "CHECK2", "CHECK3", "CHECK4")
table(sample.info.df$Genotype[tf])
table(table(sample.info.df$Genotype[!tf]))
length(sample.info.df$Genotype[!tf]) # 356
length(unique(sample.info.df$Genotype[!tf])) # 354

# save data
df.save <- data.frame("Sample_ID" = colnames(rlog.mat), t(rlog.mat), row.names = NULL)
fwrite(df.save, file = paste0(dir.save, "/ExpressionData_", ref, ".csv"))

# make a histogram
rlog.vec.zero.count <- apply(rlog.mat == 0, 2, sum)
rlog.vec.non.zero.ratio <- 1 - rlog.vec.zero.count / nrow(rlog.mat)
pdf(paste0(dir.save, "/Hist_sample_wise_zero_rlog_", ref, ".pdf"), width = 6, height = 4)
hist(rlog.vec.non.zero.ratio, 
     xlab = "Proportion of non-zero rlog2 value",
     main = paste0("Histogram of ", length(rlog.vec.non.zero.ratio), " samples"))
dev.off()



# make a histogram
rlog.vec.zero.count.v2 <- apply(rlog.mat == 0, 1, sum)
rlog.vec.non.zero.ratio.v2 <- 1 - rlog.vec.zero.count.v2 / ncol(rlog.mat)
pdf(paste0(dir.save, "/Hist_gene_wise_zero_rlog_", ref, ".pdf"), width = 6, height = 4)
hist(rlog.vec.non.zero.ratio.v2, 
     xlab = "Proportion of non-zero rlog2 value",
     main = paste0("Histogram of ", length(rlog.vec.non.zero.ratio.v2), " genes"))
abline(v = 0.5, lwd = 2, lty = 2, col = "red")
dev.off()

# make a boxplot for plate effects
df.fig <- data.frame("Number.of.zero.count" = rlog.vec.zero.count,
                     "Plate" = sample.info.df$Plate)
df.fig <- df.fig[sample.info.df$ID != "Pos", ] # remove positive control
p <- ggplot(df.fig, aes(x = Plate, y = Number.of.zero.count))
p <- p + geom_violin()
p <- p + geom_boxplot(width = 0.1)
p <- p + ggtitle("# of zero count of different plates")
ggsave(p, filename = paste0(dir.save, "/Boxplot_zero_count_per_plate_", ref, ".png"), width = 7, height = 5)

# correlation among samples
CorMat.Sample <- cor(rlog.mat)
vec.min.cor <- apply(CorMat.Sample, 1, min)
vec.lowest.cor <- apply(CorMat.Sample, 1, which.min)

# make a heatmap
pdf(paste0(dir.save, "/Heatmap_", ref, ".pdf"), width = 7, height = 7)
heatmap.2(CorMat.Sample, trace = 'none')
dev.off()

# # lowest correlation
# table(vec.lowest.cor)
# colnames(rlog.mat)[c(151, 355)] # these are the two samples
# 
# # plots
# plot(x = rlog.mat[, 151], y = rlog.mat[, 1], pch = 20,
#      xlab = colnames(rlog.mat)[151], ylab = colnames(rlog.mat)[1])
# 
# # plots
# plot(x = rlog.mat[, 355], y = rlog.mat[, 1], pch = 20,
#      xlab = colnames(rlog.mat)[355], ylab = colnames(rlog.mat)[1])

# histogram of minimun cor
pdf(paste0(dir.save, "/Hist_min_cor_", ref, ".pdf"), width = 7, height = 7)
hist(vec.min.cor, xlab = "Minimum correlation",
     main = paste0("Histogram of ", length(vec.min.cor), " samples"))
dev.off()



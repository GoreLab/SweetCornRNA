# summary of outlier removal

# packages
library(data.table)

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1"
ref <- args[2] # "B73" or "PH207"

# mkdir
dir.create("RESULT/3.2-OutlierRemoval_Summary")

# file I/O
file.in <- paste0("RESULT/3.1-OutlierRemoval/Bool_Matrix_RmOut_", ver, "_", ref, ".Rdata")
file.fig1 <- paste0("RESULT/3.2-OutlierRemoval_Summary/Hist_perGene_", ver, "_", ref, ".png")
file.fig2 <- paste0("RESULT/3.2-OutlierRemoval_Summary/Hist_perSample_", ver, "_", ref, ".png")
file.text <- paste0("RESULT/3.2-OutlierRemoval_Summary/Summary_of_outlier_removal_", ver, "_", ref, ".txt")

# load data
bool.matrix <- readRDS(file.in)

# 
myFun.count.outl <- function(x){ sum(x == "outlier_removed") }
num.outl.per.accession <- apply(bool.matrix, 1, FUN = myFun.count.outl)
num.outl.per.gene <- apply(bool.matrix, 2, FUN = myFun.count.outl)

# draw histogram.1
png(file.fig1)
hist(num.outl.per.accession,
     xlab = "# of genes identified as outlier",
     main = "Count the number of outlier \n for each accession")
dev.off()

# draw histogram.2
png(file.fig2)
hist(num.outl.per.gene,
     xlab = "# of samples identified as outlier",
     main = "Count the number of outlier \n for each gene")
dev.off()

# summary text
N1 <- prod(dim(bool.matrix))
n1 <- sum(num.outl.per.accession)
prcnt.1 <- round((n1 / N1) * 100, 2)
N2  <- ncol(bool.matrix)
m2 <- mean(num.outl.per.accession)
M2 <- max(num.outl.per.accession)
prcnt.m2 <- round((m2 / ncol(bool.matrix)) * 100, 2)
prcnt.M2 <- round((M2 / ncol(bool.matrix)) * 100, 2)
N3  <- nrow(bool.matrix)
m3 <- mean(num.outl.per.gene)
M3 <- max(num.outl.per.gene)
prcnt.m3 <- round((m3 / nrow(bool.matrix)) * 100, 2)
prcnt.M3 <- round((M3 / nrow(bool.matrix)) * 100, 2)
text.a <- paste0("Summary of the outlier remval for the ", ver, " ", ref, " dataset")
text.b <- paste0("Total number of outliers = ", n1, " (", prcnt.1, "%  of ", N1, " data points)")
text.c <- paste0("Average number of outlier genes per sample = ", round(m2, 2), " (", prcnt.m2, "%  of ", N2, " genes)")
text.d <- paste0("Maximum number of outlier genes per sample = ", round(M2, 2), " (", prcnt.M2, "%  of ", N2, " genes)")
text.e <- paste0("Average number of outlier samples per gene = ", round(m3, 2), " (", prcnt.m3, "%  of ", N3, " genes)")
text.f <- paste0("Maximum number of outlier samples per gene = ", round(M3, 2), " (", prcnt.M3, "%  of ", N3, " genes)")
text.save <- c(text.a, "", text.b, "", text.c, text.d, "", text.e, text.f)
write(text.save, file = file.text)




# summary of outlier removal

# packages
library(data.table)

# params
args <- commandArgs(trailingOnly = T)
ver <- args[1] # "v1" or "v1.1" (but "v1" is needed to be done)
ref <- args[2] # "B73" or "PH207"

# mkdir
dir.create("RESULT/3.3-OutlierRemoval_RemoveExtreme")

# file I/O
file.in <- paste0("RESULT/3.1-OutlierRemoval/PeerResiduals_RmOut_", ver, "_", ref, ".txt")

# load data
bool.data <- fread(file.in, data.table = F)
bool.matrix <- as.matrix(bool.data[, -1])
rownames(bool.matrix) <- bool.data$ID

# ratio of outliers
myFun.count.outl <- function(x){ sum(is.na(x)) }
num.outl.per.accession <- apply(bool.matrix, 1, FUN = myFun.count.outl)
ratio.outl.per.accession <- num.outl.per.accession / ncol(bool.matrix) # ratio

# accessions to be removed
acc.rm <- names(ratio.outl.per.accession)[which(ratio.outl.per.accession > 0.05)]
print(paste0("Remove ", acc.rm[1], " and ", acc.rm[2]))

# remove accessions (OVERWRITE!!)
acc.retain <- setdiff(bool.data$ID, acc.rm)
tf <- bool.data$ID %in% acc.retain
bool.data.new <- bool.data[tf, ]
fwrite(bool.data.new, file = file.in)

# gene-filtering

# library
library(data.table)

# params
args <- commandArgs(trailingOnly = T)
ref <- args[1] # "B73" or "PH207"
pr <- as.numeric(args[2]) # numeric value, from 0 to 1

# mkdir
dir.save <- "RESULT/1.2-GeneFiltering"
dir.create(dir.save)

# load data
file.in <- paste0("RESULT/1.1-LookTagSeq_MakeData/ExpressionData_", ref, ".csv")
rlog.df <- fread(file = file.in, data.table = F)
ExprMat <- as.matrix(rlog.df[, -1])
rownames(ExprMat) <- rlog.df[, 1]

# proportion of zeros
pr.zero.vec <- apply(ExprMat == 0, 2, sum) / nrow(ExprMat)

# visualize (1)
file.save <- paste0(dir.save, "/Fig01-Proportion_of_zero_rlog2_values_", ref, ".png")
png(file.save, width = 500, height = 400)
hist(x = pr.zero.vec, 
     xlab = paste0("% of zero samples among the ", nrow(ExprMat), " samples"),
     main = paste0("Histogram of ", length(pr.zero.vec), " genes"))
dev.off()

# visualize (2)
file.save <- paste0(dir.save, "/Fig02-Proportion_of_zero_rlog2_values_", ref, "_ver2.png")
png(file.save, width = 500, height = 400)
hist(x = pr.zero.vec[pr.zero.vec != 0], 
     xlab = paste0("% of zero samples among the ", nrow(ExprMat), " samples"),
     main = paste0("Histogram of ", sum(pr.zero.vec != 0), " genes with at least one zero-rlog2 value"))
dev.off()

# remove genes from the data
tf.rm <- pr.zero.vec > pr
ExprMat.filtered <- ExprMat[, !tf.rm]
write.csv(ExprMat.filtered, file = paste0(dir.save, "/ExpressionData_Filtered_", ref, ".csv"))

# save filtering result
df.save <- data.frame("Gene.ID" = names(pr.zero.vec),
                      "num.zero" = apply(ExprMat == 0, 2, sum),
                      "prop.zero" = pr.zero.vec,
                      "flag" = c("keep", "remove")[as.numeric(tf.rm) + 1],
                      row.names = NULL)
write.csv(df.save, file = paste0(dir.save, "/Summary_", ref, ".csv"), row.names = F)

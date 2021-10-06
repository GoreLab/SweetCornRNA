# package
library(data.table)

# params
args <- commandArgs(trailingOnly = T)
ref <- args[1] # "B73" or "PH207"

# filenames Input
filename.blue <- paste0("RESULT_NoHarvDate/3.1-BLUE/expression_BLUE_", ref, "_raw.csv")
filename.lastm <- paste0("RESULT_NoHarvDate/3.1-BLUE/expression_BLUE_", ref, "_message.csv")
filename.out <- paste0("RESULT_NoHarvDate/3.1-BLUE/expression_BLUE_RmErr_", ref, ".csv")
filename.summary <- paste0("RESULT_NoHarvDate/3.1-BLUE/Summary_RmErr_", ref, ".txt")

# load data
BLUE <- fread(filename.blue)
LastM <- fread(filename.lastm)

# BLUE
BLUE.mat <- as.matrix(BLUE[, -1])
rownames(BLUE.mat) <- unlist(BLUE[, 1])

# remove genes on which we get error message
tf <- LastM$last.m == "LogLikelihood Converged"
BLUE.RmOut <- BLUE.mat[tf, ]

# remove genes with obvious error
tf.ob <- apply(BLUE.RmOut < -100, 1, sum) == 0
BLUE.RmOut.v2 <- BLUE.RmOut[tf.ob, ]

# save final BLUE data
df.save <- data.frame("GeneID" = rownames(BLUE.RmOut.v2), BLUE.RmOut.v2)
fwrite(df.save, filename.out)


# make summary
n.sample.all <- ncol(BLUE.RmOut.v2)
n.gene.all <- nrow(LastM)
n.gene.conv <- sum(LastM$last.m == "LogLikelihood Converged")
lastm.non.conv <- LastM$last.m[LastM$last.m != "LogLikelihood Converged"]
table.lastm.non.conv <- table(lastm.non.conv)
text.err <- c()
if ( length(table.lastm.non.conv) != 0 ) {
   for ( i in 1:length(table.lastm.non.conv) ) {
      text.i <- paste0(names(table.lastm.non.conv[i]), ": ",  table.lastm.non.conv[i])
      text.err <- c(text.err, text.i)
   }
}

# write summary
text.summary <- c(paste0("Summary numbers for the BLUE calcualtion of ", ref, "-mapped dataset"), 
                  "",
                  paste0("The raw BLUE matrix has ", n.sample.all, " samples and ", n.gene.all, " genes"),
                  "",
                  paste0("Sumamry of the error message:"),
                  paste0("LogLikelihood Converged: ", n.gene.conv),
                  text.err,
                  "",
									"We also checked obviously wrong result",
									paste0("There are ", sum(!tf.ob), " additional wrong calculation results."),
									"",
                  paste0("Thus, the BLUE matrix has ", ncol(BLUE.RmOut.v2), " samples and ", nrow(BLUE.RmOut.v2),
                         " genes, after removing errors")
                  )
write(x = text.summary, file = filename.summary, sep = "/n")

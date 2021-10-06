# make scree plot

# mkdir
dir.create("RESULT_NoHarvDate/5.2-MakeScreePlot")

# for all cases
ref.all <- c("B73", "PH207", "Ia453")
dat.all <- c("toco", "ion", "carot")

# loop for all 
for(ref in ref.all) {
        for (dat in dat.all) {
                # make figure
                file.out.precision <- paste0("RESULT_NoHarvDate/5.1-Peer_Use25Fact/PeerResult_Use25Fact_", dat, "_", ref, "_precision.txt")
                file.out <- paste0("RESULT_NoHarvDate/5.2-MakeScreePlot/ScreePlot_Use25Fact_", dat, "_", ref, ".jpeg")
                precision <- read.delim(file.out.precision)
                variacne.from.second <- 1 / precision$precision[-1]
                jpeg(filename = file.out)
                plot(x = 2:25, y = variacne.from.second,
                     type = "l",
                     lwd = 2,
                     col = "red",
                     pch = 16,
                     xlab = "Factors",
                     ylab = "Variance",
                     main = paste0("Scree plot of PEER for ", ref, ", ", dat))
                points(x = 2:25, y = variacne.from.second, pch = 16)
                dev.off()
        }
}







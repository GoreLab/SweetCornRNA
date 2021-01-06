library(plotrix)

# mkdir
dir.save <- "RESULT/6.2-Prediction_toco_Makebarplot"
dir.create(dir.save)

# load phenotype data
Pheno.All <- read.csv("RAWDATA/Tocochromanols/toco_transformed_blups.csv")
id.all <- Pheno.All$Sample.ID

# traits
traits <- colnames(Pheno.All)[6:ncol(Pheno.All)]

# make a data frame to make barplot
df.fig <- NULL
for ( trait in traits ) {
  # input file
  folder <- "RESULT/6.1-Prediction_toco"
  file.mGRM <- paste0("PredictionResult_", trait, "_mGRM_BLUE_lmer.csv")
  file.tGRM <- paste0("PredictionResult_", trait, "_tGRM_BLUE_lmer.csv")
  file.mGRM.tGRM <- paste0("PredictionResult_", trait, "_mGRM.tGRM_BLUE_lmer.csv")
  
  # load
  dat.mGRM <- read.csv(file = paste0(folder, "/", file.mGRM))
  dat.tGRM <- read.csv(file = paste0(folder, "/", file.tGRM))
  dat.mGRM.tGRM <- read.csv(file = paste0(folder, "/", file.mGRM.tGRM))
  
  # match datasets
  id.pred <- dat.mGRM$accession_name_toco
  m <- match(id.pred, id.all)
  obs <- Pheno.All[[trait]][m]
  
  # calculate correlation
  cor.mGRM <- as.numeric(cor(x = obs, y = dat.mGRM[, -1], use = "p"))
  cor.tGRM <- as.numeric(cor(x = obs, y = dat.tGRM[, -1], use = "p"))
  cor.mGRM.tGRM <- as.numeric(cor(x = obs, y = dat.mGRM.tGRM[, -1], use = "p"))
  
  # make a data.frame
  df.tmp <- data.frame("Trait" = trait,
                       "Method" = c("G-only", "T-only", "G+T"),
                       "AVG" = c(mean(cor.mGRM), mean(cor.tGRM), mean(cor.mGRM.tGRM)),
                       "SE" = c(std.error(cor.mGRM), std.error(cor.tGRM), std.error(cor.mGRM.tGRM)))
  df.fig <- rbind.data.frame(df.fig, df.tmp)
  
  # show histogram of observed values
  png(file = paste0(dir.save, "/Hist_", trait, ".png"))
  hist(obs, main = paste0("Histogram of ", trait), xlab = "BLUPs", breaks = 20)
  dev.off()
}


# plot(x = obs, y = dat.tGRM[, 7], pch = 20)
# abline(h = mean(obs, na.rm = T), v = mean(dat.tGRM[, 7], na.rm = T), lty = 2)
# cor(x = obs, y = dat.tGRM[, 7], use = "p")
# 
# plot(x = obs, y = dat.mGRM[, 7], pch = 20)
# abline(h = mean(obs, na.rm = T), v = mean(dat.mGRM[, 7], na.rm = T), lty = 2)
# cor(x = obs, y = dat.mGRM[, 7], use = "p")
# 
# plot(x = obs, y = dat.mGRM.tGRM[, 7], pch = 20)
# abline(h = mean(obs, na.rm = T), v = mean(dat.mGRM.tGRM[, 7], na.rm = T), lty = 2)
# cor(x = obs, y = dat.mGRM.tGRM[, 7], use = "p")
# 
# t.test(cor.mGRM, cor.mGRM.tGRM)
# mean(cor.mGRM.tGRM / cor.mGRM)
# 
# tmp <- read.csv("RESULT/6.1-Prediction_toco/CrossValidationFold.csv")
# plot(x = obs, y = dat.tGRM[, 7], pch = 16, col = tmp[, 7])
# legend("topright", legend = c("Fold", 1:5), col = 1:5, pch = 16)
# abline(h = mean(obs, na.rm = T), v = mean(dat.tGRM[, 7], na.rm = T), lty = 2)
# 
# for (k in 1:5) {
#   tf.k <- tmp[, 7] == k
#   r.k <- cor(x = obs[tf.k], y = dat.tGRM[, 7][tf.k], use = "p")
#   print(r.k)
#   
#   plot(x = obs, y = dat.tGRM[, 7], pch = 16, col = "gray")
#   points(x = obs[tf.k], y = dat.tGRM[, 7][tf.k], pch = 16, col = 2)
#   abline(h = mean(obs, na.rm = T), v = mean(dat.tGRM[, 7], na.rm = T), lty = 2)
#   abline(h = mean(obs[tf.k], na.rm = T), v = mean(dat.tGRM[, 7][tf.k], na.rm = T), lty = 2, col = 2)
# }
# mean(dat.tGRM[, 7][tf.k])


# rename (short name)
df.fig$Trait[df.fig$Trait == "Total.T"] <- "SumT"
df.fig$Trait[df.fig$Trait == "Total.T3"] <- "SumT3"
df.fig$Trait[df.fig$Trait == "Total.T3_plus_T"] <- "SumTT3"
df.fig$Trait[df.fig$Trait == "ratio_aT_gT"] <- "aT/gT"
df.fig$Trait[df.fig$Trait == "ratio_dT_aT"] <- "dT/aT"
df.fig$Trait[df.fig$Trait == "ratio_dt_gT"] <- "dT/gT"
df.fig$Trait[df.fig$Trait == "ratio_aT3_gT3"] <- "aT3/gT3"
df.fig$Trait[df.fig$Trait == "ratio_dT3_aT3"] <- "dT3/aT3"
df.fig$Trait[df.fig$Trait == "ratio_dT3_gT3"] <- "dT3/gT3"
df.fig$Trait[df.fig$Trait == "ratio_TotalT_TotalT3"] <- "sumT/sumT3"
df.fig$Trait[df.fig$Trait == "dT_over_sum_gT_aT"] <- "dT/gTaT"
df.fig$Trait[df.fig$Trait == "gt_over_sum_gT_aT"] <- "gT/gTaT"
df.fig$Trait[df.fig$Trait == "dT3_over_sum_gT3_aT3"] <- "dT3/gT3aT3"
df.fig$Trait[df.fig$Trait == "gT3_over_sum_gT3_aT3"] <- "gT3/gT3aT3"
df.fig$Trait <- factor(df.fig$Trait, levels = unique(df.fig$Trait))

# make barplot (1)
pdf(paste0(dir.save, "/PredictionAccuracy.pdf"), width = 25, height = 10)
df.fig$Method <- factor(df.fig$Method, levels = c("T-only", "G-only", "G+T"))
Means <- xtabs(AVG ~ Method + Trait, data = df.fig)
bp <- barplot(height = Means, beside = TRUE,
              ylim = c(-0.2, 1), font.axis = 1, font.lab = 2,
              main = "",
              ylab = "Accuracy (Pearson's correlation)",
              xlab = "Traits",
              cex.lab = 10 / 8,
              border = "black",
              las = 1,
              legend.text = c("G-only", "T-only", "G+T"),
              args.legend = list(x = "topright", text.font = 3))
arrows(x0 = -2.1, y0 = 0, x1 = nlevels(df.fig$Trait) * (1 + nlevels(df.fig$Method)) + 2.1, length = 0)

# Function that puts error bars on barplot columns
error.bars <- function(x, y, se) {
  x <- as.vector(x)   # centers of bars on OX axis
  y <- as.vector(y)   # heights of bars (the means) on OY axis
  se <- as.vector(se) # standard errors
  for (i in 1:length(x))
    if (y[i] > 0 ) {
      arrows(x[i], y[i], x[i], y[i]+se[i], code=3, angle=90, length=0.1/2.54)
    } else {
      arrows(x[i], y[i]-se[i], x[i], y[i], code=3, angle=90, length=0.1/2.54)
    }
}
# call the function above to place error bars
error.bars(x = bp, y = Means, se = df.fig$SE)
dev.off()


# package
library(plotrix)
library(ggplot2)

# object

# mkdir
dir.save <- "RESULT/99.4-LookAtThePrediction"
dir.create(dir.save)

# load phenotype data
Pheno.All <- read.csv("RAWDATA/Tocochromanols/toco_transformed_blups.csv")
id.all <- Pheno.All$Sample.ID

# traits
traits <- colnames(Pheno.All)[6:ncol(Pheno.All)]
models <- c("mGRM", "tGRM", "mGRM.tGRM")

# define functions
myFun.CalcCor <- function(x, pred.vec, mu.vec) {
   r.vec <- rep(NA, nlevels(mu.vec))
   for ( i in 1:nlevels(mu.vec) ) {
      mu <- levels(mu.vec)[i]
      tf <- mu.vec == mu
      r <- cor(x[tf], pred.vec[tf], use = "p")
      r.vec[i] <- r
   }
   names(r.vec) <- levels(mu.vec)
   return(r.vec)
}
myFun.MakePlot <- function(x, y, xlab, ylab, main, ...) {
   lim <- range(c(x, y), na.rm = T)
   plot(x, y, xlim = lim, ylim = lim, pch = 20,
        xlab = xlab, ylab = ylab, main = main, col = mu.vec)
   abline(0, 1, lty = 2)
   legend("bottomright", legend = levels(mu.vec), pch = 20, col = 1:3)
}





# make a data frame to make barplot
df.cor.summary.all <- NULL
count <- 0
for ( trait in traits ) {
   for ( model in models ) {
      # count
      count <- count + 1
      
      # input file
      folder <- "RESULT/6.3-Prediction_toco_UseFull_UseMu"
      filename <- paste0("PredictionResult_", trait, "_", model, ".csv")
      
      # load
      dat <- read.csv(file = paste0(folder, "/", filename))
      
      # match samples
      id.pred <- dat$accession_name_toco
      m <- match(id.pred, id.all)
      obs <- Pheno.All[[trait]][m]
      mu.vec <- factor(Pheno.All$Endosperm.mutation[m], 
                       levels = c("su1", "sh2", "su1sh2"))
      
      # correlation within mutant group
      CorMat.mu <- apply(dat[, -1], 2, myFun.CalcCor, pred.vec = obs, mu.vec = mu.vec)
      
      # result
      m <- apply(CorMat.mu, 1, mean)
      s <- apply(CorMat.mu, 1, std.error)
      df.cor.summary <- data.frame("Trait" = trait,
                                   "Model" = model,
                                   "mutant.type" = names(m),
                                   "mean.cor" = m,
                                   "std.err.cor" = s,
                                   row.names = NULL)
      df.cor.summary.all <- rbind(df.cor.summary.all, df.cor.summary)
      
      # make figure for each
      if ( model == "mGRM" ) { model.print <- "G-only model" }
      if ( model == "tGRM" ) { model.print <- "T-only model" }
      if ( model == "mGRM.tGRM" ) { model.print <- "G+T model" }
      fig.name <- paste0(dir.save, 
                         "/Fig", formatC(count, width = 3, flag = "0"), "-",
                         trait, "_", model, ".png")
      png(fig.name, width = 500, height = 500)
      myFun.MakePlot(x = obs, 
                     y = dat[, 2],
                     xlab = "observed value",
                     ylab = "predicted value",
                     main = paste0("Prediction Result: ", 
                                   model.print, " for ", trait))
      dev.off()
   }
}   

# figure
tf <- df.cor.summary.all$Trait %in% c("aT", "aT3", "dT", "dT3", "gT", "gT3")
df.fig <- df.cor.summary.all[tf, ]
df.fig$Model[df.fig$Model == "mGRM"] <- "G-only"
df.fig$Model[df.fig$Model == "tGRM"] <- "T-only"
df.fig$Model[df.fig$Model == "mGRM.tGRM"] <- "G+T"
df.fig$Model <- factor(df.fig$Model, c("T-only", "G-only", "G+T"))
df.fig$err.bar.min <- df.fig$err.bar.max <- NA
for ( j in 1:nrow(df.fig) ) {
   if ( df.fig$mean.cor[j] < 0 ) {
      df.fig$err.bar.min[j] <- df.fig$mean.cor[j] - df.fig$std.err.cor[j]
      df.fig$err.bar.max[j] <- df.fig$mean.cor[j]
   } else {
      df.fig$err.bar.min[j] <- df.fig$mean.cor[j]
      df.fig$err.bar.max[j] <- df.fig$mean.cor[j] + df.fig$std.err.cor[j]
   }
}
p <- ggplot(df.fig, 
            aes(x = Model, y = mean.cor, fill = mutant.type))
p <- p + geom_bar(stat = "identity", position = "dodge", col = "black")
p <- p + geom_errorbar(aes(ymin = err.bar.min, 
                           ymax = err.bar.max), 
                       position = position_dodge(0.9), width = .3)
p <- p + facet_wrap(~ Trait, ncol = 6)
p <- p + xlab("Correlatrion coefficient")
ggsave(filename = paste0(dir.save, "/Correlation.png"), p, width = 12, height = 5)

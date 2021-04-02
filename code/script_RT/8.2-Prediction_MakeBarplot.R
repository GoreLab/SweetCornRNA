library(plotrix)
library(ggplot2)

# mkdir
dir.save <- "RESULT/8.2-Prediction_Makebarplot"
dir.create(dir.save, recursive = T)

# args
dat <- "toco"
mu <- "UseMu"

# load phenotype data
f <- paste0("RESULT/4.1-MakeBlueDatasets/PhenoData_", dat, ".csv")
PhenoData <- read.csv(f)

# traits
trait.all <- colnames(PhenoData)[3:ncol(PhenoData)]
trait <- traits[1]

# models
model.all <- c("mGRM", "tGRM.B73", "tGRM.PH207", "tGRM.both", 
               "mGRM.tGRM.B73", "mGRM.tGRM.PH207", "mGRM.tGRM.both")

# for all pairs...
df <- NULL
for ( j in 1:length(trait.all) ) {
  for ( k in 1:length(model.all) ) {
    # trait & model
    trait <- trait.all[j]
    model <- model.all[k]
    
    # load predicted values
    f <- paste0("RESULT/8.1-Prediction/PredictionResult_", 
                dat, "_", trait, "_", mu, "_", model, ".csv")
    PredData <- read.csv(f)
    
    # cor(obs, pred)
    cor.vec <- cor(x = PredData[, 2:ncol(PredData)], 
                   y = PhenoData[[trait]],
                   use = "pair")[, 1]
    df.jk <- data.frame("trait" = trait,
                        "model" = model, 
                        "rep" = names(cor.vec), 
                        "accuracy" = cor.vec,
                        row.names = NULL)
    
    # merge
    df <- rbind(df, df.jk)
  }
}
df$model <- factor(df$model, levels = model.all)

# add group
df$group <- NA
df$group[df$model == "mGRM"] <- "G-Only"
df$group[df$model %in% c("tGRM.B73", "tGRM.PH207", "tGRM.both")] <- "T-Only"
df$group[df$model %in% c("mGRM.tGRM.B73", "mGRM.tGRM.PH207", "mGRM.tGRM.both")] <- "G+T"

# make figure
for ( j in 1:length(trait.all) ) {
  trait <- trait.all[j]
  f <- paste0(dir.save, "/Boxplot_", dat, "_" , trait, "_", mu, ".png")
  df.tr <- df[df$trait == trait, ]
  p <- ggplot(df.tr, aes(x = model, y = accuracy, fill = group))
  p <- p + geom_boxplot()
  p <- p + ggtitle(trait)
  ggsave(p, filename = f, width = 10, height = 6)
}






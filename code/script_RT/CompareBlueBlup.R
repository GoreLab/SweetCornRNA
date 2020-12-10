# Example to interpret the BLUE/BLUP results

# library
library(ggplot2)

# load results
dat.BLUE <- read.csv("RESULT/2.1-BLUE_lmer/BLUE_005332.csv")
dat.BLUP <- read.csv("RESULT/2.2-BLUP_lmer/BLUP_005332.csv")

# e.g., BLUE vs BLUP
name.exp.gen <- dat.BLUP$Model.Term[14:371]

m1 <- match(name.exp.gen, dat.BLUP$Model.Term)
blup <- dat.BLUP$Est.effect[m1] + dat.BLUP$Est.effect[dat.BLUP$Model.Term == "Intercept"]

m2 <- match(name.exp.gen, dat.BLUE$Model.Term)
blue <- dat.BLUE$Est.effect[m2]

lim <- range(c(blue, blup))
plot(x = blue, y = blup, 
     xlab = "BLUE", ylab = "BLUP",
     main = "BLUE vs BLUP",
     xlim = lim, ylim = lim, pch = 20)
abline(0, 1, lty = 2)

# Estimated variance components in BLUP
name.varcomp <- c("var.geno", "var.range", "var.block", "var.column", "var.row",
                  "var.plate", "var.extract", "var.resid")
rename.varcomp <- c("Genotype", "Range", "Block", "Col", "Row", "Plate", "Extract.Date", "Resid")
m <- match(name.varcomp, dat.BLUP$Model.Term)
df.fig <- dat.BLUP[m, ]
df.fig$Model.Term <- factor(rename.varcomp, levels = rename.varcomp)

p <- ggplot(df.fig, aes(x = Model.Term, y = Est.effect))
p <- p + geom_bar(stat = "identity")
p <- p + ylab("Estimated Variance")
p <- p + ggtitle("Estimated Variacne in BLUP model")
p

# e.g., visualize plate effect in BLUE model
name.plate <- c("plate1", "plate2", "plate3", "plate4", "plate5")
m <- match(name.plate, dat.BLUE$Model.Term)
df.fig <- dat.BLUE[m, ]

p <- ggplot(df.fig, aes(x = Model.Term, y = Est.effect))
p <- p + geom_bar(stat = "identity")
p <- p + ylab("Estimated Effect")
p <- p + ggtitle("Plate effect in BLUE model")
p









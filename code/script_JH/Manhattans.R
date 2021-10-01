library(qqman)
library(grid)
library(gridGraphics)
library(ggplot2)
library(ggpubr)

res.twas <- read.csv("./output/NoHarvDate/TWAS_all_genes.csv") # no harvest date B73 all traits all genes


# remove ctg/Mt/Pt
res.twas.fig <- res.twas %>%
  filter(chr %in% 1:10) %>%
  mutate(chr = as.integer(chr),
         P = 10 ^ (-neg.log.P),
         trait = as.character(trait)) %>%
  rename(SNP = RefGen_v4.Gene.ID,
         CHR = chr,
         BP = start) %>%
  dplyr::select(SNP, CHR, BP, trait, neg.log.P, P) %>%
  drop_na()

# make figures
for(i in 1:length(unique(res.twas.fig$trait))){
  trait.i <- unique(res.twas.fig$trait)[i]
  res.twas.fig.i <- res.twas.fig %>%
    filter(trait == trait.i)

  # get p-adjust.fdr = 0.05 line
  p.i <- sort(res.twas.fig.i$P)
  p.fdr.i <- p.adjust(p.i, method = "fdr")
  n.signif.i <- sum(p.fdr.i < 0.05)
  if (n.signif.i == 0){
    cutoff.i <- NULL
    n.log.c.i <- FALSE
    } else if (n.signif.i >= 1){
      cutoff.i <- (p.i[n.signif.i] + p.i[n.signif.i + 1]) / 2
      n.log.c.i <- -log10(cutoff.i)
    }

  # make plot
  plotname.i <- paste0(trait.i, ".plot")
  plot.i <- res.twas.fig.i %>%
  manhattan(.,
            chr = "CHR", bp = "BP", snp = "SNP", p = "P",
            main = paste0("TWAS results for ", trait.i),
            suggestiveline = FALSE,
            genomewideline = n.log.c.i,
            annotatePval = cutoff.i)
  print(paste0(trait.i, " has ", n.signif.i, " genes above the FDR threshold"))
  assign(plotname.i, plot.i)

  # transform plot to ggplot object so can save as png
  plot.record.i <- recordPlot()
  grid.plot.record.i <- grid.grabExpr(grid.echo(plot.record.i))
  ggsave(grid.plot.record.i, filename = paste0("./output/NoHarvDate/Manhattans/manhattan_", trait.i, ".png"),
         device = "png", units = "in", width = 12, height = 9)
}




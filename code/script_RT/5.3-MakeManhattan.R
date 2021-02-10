# source


# params
library(qqman)

# mkdir
dir.save <- "RESULT/5.3-MakeManhattan"
dir.create(dir.save, recursive = T)

# object
traits <- c("aT", "aT3", "dT", "dT3", "gT", "gT3", "Total.T", "Total.T3", "Total.T3_plus_T",
            "dT_over_sum_gT_aT", "dT3_over_sum_gT3_aT3", 
            "gt_over_sum_gT_aT", "gT3_over_sum_gT3_aT3", 
            "ratio_aT_gT", "ratio_aT3_gT3",
            "ratio_dT_aT", "ratio_dt_gT",
            "ratio_dT3_aT3", "ratio_dT3_gT3",
            "ratio_TotalT_TotalT3")

# loop for all traits
fig.file <- paste0(dir.save, "/TWAS_manhattan_all.pdf")
pdf(fig.file, width = 8, height = 4)
for ( trait in traits ) {
  # load data
  file <- paste0("RESULT/5.2-TWAS_toco_ver2/TwasResult_", trait, ".csv")
  res.twas <- read.csv(file)
  
  # remove ctg/Mt/Pt
  tf <- res.twas$Chr %in% 1:10
  res.twas.fig <- res.twas[tf, ]
  res.twas.fig$Chr <- as.integer(res.twas.fig$Chr)
  
  # Make the Manhattan plot on the gwasResults dataset
  res.twas.fig$P <- 10 ^ (-res.twas.fig$neg.log.P)
  colnames(res.twas.fig) <- c("SNP", "CHR", "BP", "neg.log.P", "P")
  
  # get p-adjust.fdr = 0.05 line
  p <- sort(res.twas.fig$P)
  p.fdr <- p.adjust(p, method = "fdr")
  n.signif <- sum(p.fdr < 0.05)
  if ( n.signif == 0 ) { cutoff <- NULL; n.log.c <- FALSE }
  if ( n.signif >= 1 ) { cutoff <- (p[n.signif] + p[n.signif+1]) / 2; n.log.c <- -log10(cutoff) }

  # make figure
  manhattan(res.twas.fig,
            chr = "CHR", bp = "BP", snp = "SNP", p = "P",
            main = paste0("TWAS result ", trait),
            suggestiveline = FALSE,
            genomewideline = n.log.c,
            annotatePval = cutoff)
  print(n.signif)
}
dev.off()


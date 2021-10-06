# 14.6-Kernel_Manhattan
# https://www.r-graph-gallery.com/101_Manhattan_plot.html
library(tidyverse)
library(qqman)
library(grid)
library(gridGraphics)
library(ggplot2)
library(ggpubr)
library(ggtext)

# endosperm mutation results
res.twas <- read.csv("./RESULT/9.7-TWAS_mutant_type/TwasResult_toco_B73_endosperm.mutation.csv") %>%  # no harvest date B73 all traits all genes
  mutate(Chr = as.character(Chr)) %>%
  filter(Chr %in% as.character(1:10)) %>%
  mutate(P = 10 ^ (-neg.log.P))

# calculate FDR cutoff
p.i <- sort(res.twas$P)
p.fdr.i <- p.adjust(p.i, method = "fdr")
n.signif.i <- sum(p.fdr.i < 0.05)
if (n.signif.i == 0){
  cutoff.i <- NULL
  n.log.c.i <- FALSE
} else if (n.signif.i >= 1){
  cutoff.i <- (p.i[n.signif.i] + p.i[n.signif.i + 1]) / 2
  n.log.c.i <- -log10(cutoff.i)
}


twas.prepped <- res.twas %>%
  mutate(Chr = factor(Chr, levels = 1:10)) %>%
  # Compute chromosome size
  group_by(Chr) %>%
  summarise(chr_len = max(Pos)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(res.twas, ., by = c("Chr" = "Chr")) %>%
  mutate(Chr = factor(Chr, levels = as.character(1:10))) %>%
  # Add a cumulative position of each SNP
  arrange(Chr, Pos) %>%
  mutate(Pos = as.numeric(Pos),
         tot = as.numeric(tot)) %>%
  mutate(BPcum=Pos+tot)


axisdf <- twas.prepped %>%
  group_by(Chr) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

sc_genes <- twas.prepped %>% filter(Gene == "Zm00001d044129" | Gene == "Zm00001d049753") %>%
  cbind(nickname = c("sh2", "su1")) %>% mutate(label = paste0(nickname, " (", Gene, ")"))


# plot
kernel_manhattan <- ggplot(twas.prepped, aes(x=BPcum, y=neg.log.P)) +

  # Show all points
  geom_point( aes(color=Chr), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +

  # FDR
  geom_hline(yintercept = n.log.c.i, col = "firebrick", linetype = "dashed") +
  geom_text(data = data.frame(x = 1, y = n.log.c.i, label = "0.05 FDR"),
            aes(x = x, y = y, label = label, vjust = -1, hjust = "inward"), col = "firebrick") +

  # Label sh2 and su1
  geom_point(data = sc_genes, aes(x = BPcum, y = neg.log.P), col = "black", size = 1)+
  geom_text(data = sc_genes[1,],
            aes(x = BPcum, y = neg.log.P, label = nickname, hjust = 1.5, vjust = -0.1),
            fontface = "italic") +
  geom_text(data = sc_genes[2,],
            aes(x = BPcum, y = neg.log.P, label = nickname, hjust = -0.5, vjust = -0.1),
            fontface = "italic") +
  geom_vline(xintercept = sc_genes$BPcum[1], linetype = "dashed") +
  geom_vline(xintercept = sc_genes$BPcum[2], linetype = "dashed") +

  # custom X axis:
  scale_x_continuous(label = axisdf$Chr, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(twas.prepped$neg.log.P) + 2)) +     # remove space between plot area and x axis

  # titles
  labs(x = "Chromosome", y = "-log<sub>10</sub>*P*") +

  # Customize the theme:
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown()
  )
kernel_manhattan

ggsave(kernel_manhattan, filename = "./RESULT/14.6-Kernel_Manhattan/Figure1_manhattan_kernel_type.png",
       device = "png",
       units = "in", width = 6, height = 4
       )


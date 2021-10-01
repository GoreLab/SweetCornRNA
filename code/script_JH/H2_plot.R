
library(tidyverse)

carot.h2 <- read.csv("./data/SC_carotenoid_H2s.csv")
toco.h2 <- read.csv("./data/SC_tocochromanol_H2s.csv")

h2s <- rbind(carot.h2, toco.h2)

hist(h2s$HeritabilityEstimate)

h2s %>%
  ggplot(aes(x = HeritabilityEstimate)) +
  geom_histogram(bins = 15) + theme_bw() +
  labs(x = "Broad sense heritability estimate", y = "Count",
       title = "Broad sense heritability estimates for 20 tocochromanol and 19 carotenoid traits")

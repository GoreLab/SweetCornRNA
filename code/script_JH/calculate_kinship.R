# calculate kinship matrix with GAPIT
# 12/15/2020
# Jenna Hershberger
# jmh579@cornell.edu

source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

hmp_input <- read.table("./data/Ion_11K_401_v4.hmp.txt", head = FALSE)
setwd("./ouput/GAPIT")

myGAPIT <- GAPIT(
  G = hmp_input,
  kinship.algorithm = "VanRaden",
  file.output = T
  )

# 9/15/2021

Ia453_0.05FDR <- read.csv("./output/NoHarvDate/noharv_Ia453_overallFDR_0.05.csv")
B73_0.05FDR <- read.csv("./output/NoHarvDate/noharv_B73_overallFDR_0.05.csv")

intersect(Ia453_0.05FDR$B73v4_id, B73_0.05FDR$RefGen_v4.Gene.ID)

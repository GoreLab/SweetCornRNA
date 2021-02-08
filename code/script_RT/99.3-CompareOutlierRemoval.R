# 
library(data.table)
Res1 <- fread("RESULT/4.4-OutlierRemoval/StudentizedRes_UseOptFact_BLUE_lmer_alpha_5.txt",
              data.table = F)
Res2 <- fread("RESULT/4.4-OutlierRemoval_asreml/StudentizedRes_UseOptFact_BLUE_lmer_alpha_5.txt",
              data.table = F)
Res1[1:5, 1:5]
Res2[1:5, 1:5]
all(which(is.na(Res1)) == which(is.na(Res2)))



Res1 > 

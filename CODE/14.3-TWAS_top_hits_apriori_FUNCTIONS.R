# Functions for compiling TWAS results and matching with annotations, a priori candidate lists, and NAM JL QTL intervals
# 5/11/2021

library(tidyverse)
library(magrittr)
library(fs)


read_csv_adapted_B73 <- function(path){
  file.df <- read.csv(file = path)
  pat <- "(?<=B73_)(.*)(?=.csv)"
  trait.name <- str_extract(path, pat)
  file.df$trait <- trait.name
  return(file.df)
}

process_infile_B73 <- function(input.df, exclude.vec, apriori.df, annotations.df){
  output.df <- input.df %>%
    rename(RefGen_v4.Gene.ID = Gene) %>%
    filter(!trait %in% exclude.vec) %>%
    arrange(trait, -neg.log.P) %>%
    group_by(trait) %>%
    mutate(overall.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr"),
           rank = rank(-neg.log.P),
           rank_percent = rank/18765*100) %>%
    full_join(apriori.df) %>%
    filter(!trait %in% exclude.vec) %>%
    dplyr::select(RefGen_v4.Gene.ID:rank_percent, Apriori, Diepenbrock_21, Wu_21) %>%
    left_join(annotations.df) %>%
    dplyr::select(RefGen_v4.Gene.ID, neg.log.P:rank_percent, chr:logic_name, Apriori, Diepenbrock_21,
                  Wu_21, QTL_ID:NAM.JL.QTL.traits, Marker.at.the.Peak:Overlap) %>%
                # Wu_21, QTL_ID:NAM.JL.QTL.traits, A.Priori.Gene.in.Support.Interval, Overlap)  %>%
    distinct()
}

match_JL_to_TWAS <- function(SC_TWAS, NAM_JL){
  SC_TWAS <- as_tibble(SC_TWAS)
  for(i in 1:nrow(SC_TWAS)){
    cat(i, ", ")
    keep.NAM_JL <- NAM_JL %>%
      filter(`Chr` == unlist(SC_TWAS[i,"chr"]),
             # start gene <= end interval and end gene >= start interval
             (unlist(SC_TWAS[i,"start"]) <= `Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)`) & # switch to CSI
               (unlist(SC_TWAS[i,"end"]) >= `Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)`)) %>%
      mutate(Overlap = ifelse((`Support Interval (α = 0.01), Left Bound (bp) (RefGen v4)` <= unlist(SC_TWAS[i,"start"])) &
                                (`Support Interval (α = 0.01), Right Bound (bp) (RefGen v4)` >= unlist(SC_TWAS[i,"end"])),
                              "Complete", "Partial"))

    if(nrow(keep.NAM_JL) < 1){
      keep.NAM_JL <- t(matrix(data = rep(NA, ncol(NAM_JL)+1)))
      colnames(keep.NAM_JL) <- c(colnames(NAM_JL), "Overlap")
    }
    new.i <- cbind(SC_TWAS[i,], as_tibble(keep.NAM_JL))

    if(i == 1){
      new.SC_TWAS <- new.i
    } else{
      new.SC_TWAS <- rbind(new.SC_TWAS, new.i)
    }
  }

  return(new.SC_TWAS)
}


########## Ia453 ##############

read_csv_adapted_Ia453 <- function(path){
  file.df <- read.csv(file = path) %>%
    dplyr::select(Gene, neg.log.P)
  pat <- "(?<=Ia453_)(.*)(?=.csv)"
  trait.name <- str_extract(path, pat)
  file.df$trait <- trait.name
  return(file.df)
}

process_infile_Ia453 <- function(input.df, exclude.vec, annotations.df){
  output.df <- input.df %>%
    rename(Ia453_id = Gene) %>%
    filter(!trait %in% exclude.vec) %>%
    arrange(trait, -neg.log.P) %>%
    group_by(trait) %>%
    mutate(overall.fdr = p.adjust((10 ^ (-neg.log.P)), method = "fdr"),
           rank = rank(-neg.log.P),
           rank_percent = rank/18477*100) %>%
    left_join(annotations.df, by = "Ia453_id") %>%
    filter(!trait %in% exclude.vec) %>%
    dplyr::select(Ia453_id, B73v4_id, neg.log.P:rank_percent, Ia453_chr:logic_name, Apriori, Diepenbrock_21,
                  Wu_21, QTL_ID:NAM.JL.QTL.traits, Marker.at.the.Peak:Overlap) %>%
    distinct()
}


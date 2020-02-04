#!/usr/bin/env Rscript

# GDD functions
# Jenna Hershberger
# jmh579@cornell.edu
# 02/04/2020

GDDtotals <- function(start_date_list, end_date_list, temp.df){
  GDD.list <- c()
  for(i in 1:length(start_date_list)){
    # date only > (and not =) because Paul sends temperature updates for the previous day
    # example// temperature on July 29th in the input file was *recorded* on the 28th but
    # *reported* on the 29th. Any ears pollinated on the 29th did not experience those temps.
    total.GDD <- temp.df %>% dplyr::filter(date > start_date_list[i]) %>%
      dplyr::filter(date <= end_date_list[i]) %>%
      summarize(total.GDD = sum(GDD, na.rm = T)) %>% unlist()
    GDD.list[i] <- total.GDD
  }
  return(GDD.list)
}

calculateGDD <- function(max.temp = 86, min.temp = 50){
  if(is.na(max.temp)|is.na(min.temp)){
    return(NA)
  }
  if(max.temp > 86){
    max.temp <- 86
  }
  if(min.temp > 86){
    min.temp <- 86
  }
  if(min.temp < 50){
    min.temp <- 50
  }
  if(max.temp < 50){
    max.temp <- 50
  }
  GDD <- ((max.temp + min.temp) / 2) - 50
  return(GDD)
}

# vectorized GDD function
lcalculateGDD <- function(max.temp.list, min.temp.list){
  if(length(max.temp.list) != length(min.temp.list)){
    stop("Must have same number of minimum and maximum temperatures!")
  }
  GDD.list <- c()
  for(i in(1:length(max.temp.list))){
    GDD.list[i] <- calculateGDD(max.temp.list[i], min.temp.list[i])
  }
  return(GDD.list)
}

calcDAP <- function(start_date_list, end_date_list){
  DAP.list <- c()
  for(i in 1:length(start_date_list)){
    DAP.list[i] <- end_date_list[i] - start_date_list[i]
  }
  return(DAP.list)
}


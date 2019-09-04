#!/usr/bin/env Rscript
 
# Fresh harvest @ 400 GDD script
# Jenna Hershberger
# jmh579@cornell.edu
# 08/06/2019

# TODO remove doubles from output
intercross = F

# required input:
# args[1] = working directory ex// "~/Desktop/SweetCorn/2019"
# args[2] = cross file path ex// "~/Desktop/crosses_jenna_2019-07-12_10:34:03.csv"
# args[3] = output path for running list of harvested plots ex// "~/Desktop/SweetCorn/2019/running_list_harvested_2019-07-15.csv"

# example:
# Rscript GDD_app.R "~/Desktop/crosses_jenna_2019-07-29_12/34/03.csv" "~/Desktop/SweetCorn/2019" NULL
library(tidyverse)
library(googledrive)
library(readxl)
library(lubridate)

# for debugging
#args = c("~/Desktop/crosses_jenna_2019-07-12_10:34:03.csv", "~/Desktop/SweetCorn/2019", "~/Desktop/SweetCorn/2019/running_list_harvested_2019-07-15.csv")  

args = commandArgs(trailingOnly=TRUE)

# test if there are three arguments: if not, return an error
if (length(args)<1) {
  stop("An input argument is required:\n1) Working directory\nwith optional second and third arguments:\n
       2) Cross file path\n
       3) Path for running list of harvested plots.\nType 'NA' for the second argument
       to leave it blank if you just want to enter the first and third.", call.=FALSE)
}

setwd(args[1]) #"~/Desktop/SweetCorn/2019"
working.directory <- getwd()

cross.path <- args[2] #"~/Desktop/crosses_jenna_2019-07-12_10:34:03.csv"

harvested.path <- args[3] #"~/Desktop/SweetCorn/2019/running_list_harvested_2019-08-21.csv"

#### functions ####
## Calculate GDD
# GDD = The Daily Average Temp (°F) = (Daily Max Temp °F + Daily Min Temp °F) / 2
# Daily Corn GDD (°F) = Daily Average Temperature °F - 50 °F
# Constraints on maximum and minimum temperatures are used to eliminate the effect of low or high 
#   temperatures that prevent or retard growth. For corn these constraints are:
#   - If the daily max and/or min Temp < 50 °F (10 °C), it's set equal to 50 °F (10 °C).
#   - If the daily max and/or min Temperature > 86 °F (30 °C), it's set equal to 86 °F (30 °C)

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

# Sum GDDs (vectorized)
sumGDDs <- function(pol_date.list, temp.df){
  GDD.list <- c()
  for(i in 1:length(pol_date.list)){
    # date only > (and not =) because Paul sends temperature updates for the previous day
    # example// temperature on July 29th in the input file was *recorded* on the 28th but 
    # *reported* on the 29th. Any ears pollinated on the 29th did not experience those temps.
    total.GDD <- temp.df %>% filter(date > pol_date.list[i]) %>% 
      summarize(total.GDD = sum(GDD, na.rm = T)) %>% unlist()
    GDD.list[i] <- total.GDD
  }
  return(GDD.list)
}


#### input data and run calculations ####
# Get cross list 
if(intercross){
  cross.input <- read.csv(cross.path)
  cross <- cross.input %>% mutate(pol_date = as.Date(timestamp, origin = "1899-12-30")) %>% 
    rename(plot_id = female) %>% dplyr::select(pol_date, plot_id, cross_count)
} else{
  drive_download(file = "crosses_2019", overwrite = T) # saves to working directory
  cross.input <- read_excel("crosses_2019.xlsx")
  cross <- cross.input %>% 
    mutate(pol_date = as.Date(pollination_date, format = "%m/%d/%Y")) %>% 
    dplyr::select(pol_date, plot_id, cross_count) %>% drop_na(plot_id) %>% arrange(pol_date, plot_id)
}

first.pollinations <- cross[match(unique(cross$plot_id), cross$plot_id),]



# Get list of previously harvested plots
if(is.null(harvested.path)|is.na(harvested.path)) {
  harvested.input <- NA
  } 

if(!is.null(harvested.path)&!is.na(harvested.path)){
  harvested.input <- read.csv(harvested.path)
}


# Get temperature file and calculate GDDs
drive_download(file = "Temperature_2019.xlsx", overwrite = T) # saves to working directory
temp.input.df <- read_excel("Temperature_2019.xlsx")
temp.df <- transform(temp.input.df, max.temp = as.numeric(max.temp),
                     min.temp = as.numeric(min.temp))
temp.df <- temp.input.df %>% mutate(GDD = lcalculateGDD(max.temp, min.temp)) %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

# calculate totals
running.totals <- first.pollinations %>% #TODO changed this
  mutate(total.GDD = sumGDDs(pol_date, temp.df)) 
harvest.today <- running.totals %>% mutate(harvest.date = today()) %>% 
  filter(total.GDD >= 400) 
if(!is.na(harvested.input)){ # TODO fix me
  harvest.today <- harvest.today %>% filter(!plot_id %in% as.character(harvested.input$plot_id))
}
  

# generate new list of harvested plots including those in harvest.today
if(is.na(harvested.input)){
  new.harvested <- harvest.today
} else{
  new.harvested <- rbind(harvested.input, harvest.today)
}
cat("Yesterday's weather:\n")
temp.df.2 <- temp.df %>% dplyr::select(-other)
print(temp.df[nrow(na.omit(temp.df.2)),])
cat("\n")

harvested.total <- cross.input %>% dplyr::select(plot_id) %>% distinct() %>% nrow()
cat(paste0( "\n", harvested.total, " plots have been pollinated so far.\n"))

#### save to .csv files ####
if(nrow(harvest.today) > 0){
  write.csv(harvest.today, paste0(working.directory, "/harvest_", today(), ".csv"), row.names = F)
  write.csv(new.harvested, paste0(working.directory, "/running_list_harvested_", today(), ".csv"), row.names = F)
  cat(paste0(nrow(harvest.today), " plot(s) to harvest today.\n
             Running list of harvested plots has been updated.\n\n"))
} else{
  cat("No plots to harvest today.\nRunning list of harvested plots has not changed.\n")
}
write.csv(running.totals, paste0(working.directory, "/GDD_totals_", today(), ".csv"), row.names = F)
cat("GDD totals have been updated:\n\n")
print(running.totals)

# est.date.next.pol <- 400 - # TODO finish this
# cat("\nEstimated date of next pollination: ")
# print(mean(temp.df$GDD, na.rm = T)) 

# save to google drive
#drive_upload(media = harvest.today.path, path = "~/Work/SweetCornPollinations2019/")

########## Google form version ########## 
# responses1 <- read_excel("~/Downloads/Untitled form (Responses).xlsx")
# responses2 <- responses1 %>% mutate(Plot = NA) %>% dplyr::select(-Pollinated)
# for(i in 1:nrow(responses1)){
#   plots.list <- unlist(strsplit(as.character(responses1[i, "Pollinated"]), split = "\n"))
#   for(j in 1:length(plots.list)){
#     responses2[(i+j-1),"Timestamp"] <- date(responses1$Timestamp[i]) # this changes dttm to a random number...
#     responses2[(i+j-1),"Plot"] <- plots.list[j]
#   }
# }
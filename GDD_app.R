#!/usr/bin/env Rscript
 
# Fresh harvest @ 400 GDD script
# Jenna Hershberger
# jmh579@cornell.edu
# 09/09/2019

# required input:
# args[1] = working directory ex// "~/Desktop/SweetCorn/2019"
# args[2] = Last pollination date harvested ex// "2019-08-10"

# example:
# Rscript GDD_app.R "~/Desktop/SweetCorn/2019" "2019-08-10"
library(tidyverse)
library(googledrive)
library(readxl)
library(lubridate)

# for debugging
#args = c("~/Desktop/SweetCorn/2019", "2019-08-10")  

args = commandArgs(trailingOnly=TRUE)

# test if there are three arguments: if not, return an error
if (length(args)<1) {
  stop("An input argument is required:\n1) Working directory\nwith optional second argument:\n
       2) Last pollination date harvested (example: '2019-08-10')\n", call.=FALSE)
}

setwd(args[1]) #"~/Desktop/SweetCorn/2019"
working.directory <- getwd()

last.harvested.date <- args[2] #"2019-08-10"


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
drive_download(file = "crosses_2019", overwrite = T) # saves to working directory
cross <- read_excel("crosses_2019.xlsx", na = "NA") %>% 
  mutate(pollination_date = as.Date(pollination_date, format = "%m/%d/%Y")) %>% 
  mutate(harvested = as.Date(harvested, format = "%m/%d/%Y"))
# Filter to get list of previously harvested plots
first.pollinations <- cross[match(unique(cross$plot_id), cross$plot_id),]
first.pol.dates <- first.pollinations %>% dplyr::select(pollination_date, harvested) %>% distinct() # TODO this leaves 2 rows if there are any non-harvested plots on that day
harvested <- cross %>% filter((harvested <= today()))

# Get temperature file and calculate GDDs
drive_download(file = "Temperature_2019.xlsx", overwrite = T) # saves to working directory
temp.input.df <- read_excel("Temperature_2019.xlsx")
temp.df <- transform(temp.input.df, max.temp = as.numeric(max.temp),
                     min.temp = as.numeric(min.temp))
temp.df <- temp.input.df %>% mutate(GDD = lcalculateGDD(max.temp, min.temp)) %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

# calculate GDD totals for each pollination date and write to .csv
running.totals <- first.pol.dates %>% 
  mutate(total.GDD = sumGDDs(pollination_date, temp.df)) %>% 
  filter(pollination_date >= last.harvested.date | pollination_date < last.harvested.date & !is.na(harvested))
write.csv(running.totals, paste0("~/Desktop/SweetCorn/2019/running.totals.", today(), ".csv"), row.names = F)

# print yesterday's weather and running totals of GDDs for each pollination date
cat("Yesterday's weather:\n")
temp.df.2 <- temp.df %>% dplyr::select(-other)
print(temp.df[nrow(na.omit(temp.df.2)),])
cat("\nGDD running totals:")
print(as.data.frame(running.totals))
cat("\n")

# print count of those harvested
pollinated.total <- cross %>% dplyr::select(plot_id) %>% distinct() %>% nrow()
cat(paste0( "\n", nrow(harvested), " plots have been harvested so far \nout of the ", pollinated.total, " plots that have been pollinated.\n"))


#### create and write harvest checklists ####
# Create full list of plots with each harvest event as a column
df2 <- as.data.frame(matrix(nrow = 462, ncol = 5, data = NA))
colnames(df2) <- c("plot_id", "pol_1", "pol_2", "pol_3", "pol_4")
df2[,1] <- c(paste0("19A000", 1:9), paste0("19A00", 10:99), paste0("19A0", 100:462))

cross$plot_number <- as.numeric(unlist(str_split(cross$plot_id, "A"))[c(FALSE, TRUE)])
for(i in 1:nrow(cross)){
  if(is.na(df2[cross$plot_number[i],2])){
    df2[cross$plot_number[i],2] <- cross$pollination_date[i]
  } else if(is.na(df2[cross$plot_number[i],3])){
    df2[cross$plot_number[i],3] <- cross$pollination_date[i]
  } else if(is.na(df2[cross$plot_number[i],4])){
    df2[cross$plot_number[i],4] <- cross$pollination_date[i]
  } else if(is.na(df2[cross$plot_number[i],5])){
    df2[cross$plot_number[i],5] <- cross$pollination_date[i]
  }
}
df2[,2] <- as.Date(df2[,2], origin = "1970-1-1")
df2[,3] <- as.Date(df2[,3], origin = "1970-1-1")
df2[,4] <- as.Date(df2[,4], origin = "1970-1-1")
df2[,5] <- as.Date(df2[,5], origin = "1970-1-1")

# Don't need to update pol.date.by.plot unless new pollinations are made
# write.csv(df2, "~/Desktop/SweetCorn/2019/pol.date.by.plot.20190903.csv", row.names = F, na = "")

# harvests by plot
# remove plots with no pollinations and those that have been harvested already
harvests.by.plot <- df2 %>% filter(!plot_id %in% harvested$plot_id) %>% filter(!is.na(pol_1))
write.csv(harvests.by.plot, paste0("~/Desktop/SweetCorn/2019/harvests.by.plot.", today(), ".csv"), row.names = F, na ="")

# harvests by date
# crosses sheet with harvested plots removed
h.plots <- harvested$plot_id
harvests.by.date <- cross %>% filter(!plot_id %in% h.plots) %>% filter(pollination_date > last.harvested.date) %>% 
  arrange(pollination_date, plot_number)
write.csv(harvests.by.date, paste0("~/Desktop/SweetCorn/2019/harvests.by.date.", today(), ".csv"), row.names = F, na ="")


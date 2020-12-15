

# Plot GDD vs DAP for 2014, 2015, 2019

crosses2019 <- read_excel("~/Desktop/SweetCorn/2019/crosses_2019.xlsx", na = "NA") %>% 
  mutate(pollination_date = as.Date(pollination_date, format = "%m/%d/%Y")) %>% 
  mutate(harvested = as.Date(harvested, format = "%m/%d/%Y"))
sofar2019 <- crosses2019 %>% filter(!is.na(harvested)) %>% 
  dplyr::select(pollination_date, harvested) %>% distinct()

temp.input.df <- read_excel("Temperature_2019.xlsx")
temp.df <- transform(temp.input.df, max.temp = as.numeric(max.temp),
                     min.temp = as.numeric(min.temp))
temp.df <- temp.input.df %>% mutate(GDD = lcalculateGDD(max.temp, min.temp)) %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"))


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

sofar2019$GDD_at_harvest <- GDDtotals(start_date_list = sofar2019$pollination_date, 
                                      end_date_list = sofar2019$harvested, temp.df = temp.df)
sofar2019$DAP <- calcDAP(start_date_list = sofar2019$pollination_date, 
                         end_date_list = sofar2019$harvested)


dates.2015 <- read_excel("~/Desktop/matt_2014_2015_GDD/Harvest.xlsx", sheet = "Sheet2")
for(i in 1:nrow(dates.2015)){
  year(dates.2015$`First Pollination`[i]) <- 2015
  year(dates.2015$`First Harvest Date`[i]) <- 2015
}
weather.2015 <- read.csv("~/Desktop/matt_2014_2015_GDD/2015 Gore Musgrave Weather by day.csv", header = T) %>% 
  mutate(Date = as.Date(Date, format= "%m/%d/%Y")) %>% dplyr::mutate(max.temp=as.numeric(High)) %>% 
  dplyr::mutate(min.temp=as.numeric(Low)) %>% mutate(GDD = lcalculateGDD(max.temp, min.temp)) %>% 
  dplyr::rename(date = Date)
for(i in 1:nrow(weather.2015)){
  year(weather.2015$date[i]) <- 2015
}

dates.2015$GDD_at_harvest <- GDDtotals(start_date_list = dates.2015$`First Pollination`, 
                                      end_date_list = dates.2015$`First Harvest Date`, 
                                      temp.df = weather.2015)
dates.2015$DAP <- calcDAP(start_date_list = dates.2015$`First Pollination`, 
                         end_date_list = dates.2015$`First Harvest Date`)


# no pollination dates for 2014?
GDD.2014 <- read.csv("~/Desktop/matt_2014_2015_GDD/Pollination_2014.csv")
dates.2014 <- read_excel("~/Desktop/matt_2014_2015_GDD/Harvests_20180201.xlsx", sheet = "Sheet2") 
colnames(dates.2014)[1] <- "Row.ID"
colnames(dates.2014)[2] <- "harvest.date"
dates.2014 <- dates.2014 %>% 
  mutate(harvested = as.Date(harvest.date, format = "%m-%d"))
for(i in 1:nrow(dates.2014)){
  year(dates.2014$harvested[i]) <- 2014
}


tomerge.2015 <- dates.2015 %>% dplyr::rename(pollination_date = `First Pollination`) %>% 
  dplyr::rename(harvested = `First Harvest Date`) %>% 
  dplyr::select(pollination_date, harvested, DAP, GDD_at_harvest) %>% mutate(year = 2015)
tomerge.2019 <- sofar2019 %>% 
  mutate(year = 2019) %>% 
  dplyr::select(pollination_date, harvested, DAP, GDD_at_harvest, year)
tomerge.2014 <- dates.2014 %>% 
  mutate(year = 2014) %>% mutate(pollination_date = NA) %>% 
  mutate(DAP = NA) %>% 
  full_join(GDD.2014, by = "Row.ID") %>% 
  dplyr::rename(GDD_at_harvest = First.Harvest.GDD) %>% 
  dplyr::select(pollination_date, harvested, DAP, GDD_at_harvest, year) %>% 
  filter(!is.na(harvested)) %>% filter(harvested > "2014-07-01")
# some dates in 2014 are listed as 01-01 in the input file. filtered them out.

tomerge.2014$GDD_at_harvest <- as.numeric(as.character(tomerge.2014$GDD_at_harvest))
tomerge.2015$GDD_at_harvest <- as.numeric(tomerge.2015$GDD_at_harvest)
tomerge.2019$GDD_at_harvest <- as.numeric(tomerge.2019$GDD_at_harvest)

full_dataset <- rbind(tomerge.2014, tomerge.2015, tomerge.2019)
full_dataset$year <- as.factor(full_dataset$year)

ggplot(full_dataset, aes(x = harvested, y = GDD_at_harvest, color = year)) + geom_point()
hist(tomerge.2015$harvested, breaks = 50)
hist(tomerge.2015$GDD_at_harvest, breaks = 50)
hist(tomerge.2014$harvested, breaks = 50)
min(tomerge.2014$harvested)


full_same_year <- full_dataset %>% mutate(harvested = as.Date(harvested)) %>% distinct()
for(i in 1:nrow(full_same_year)){
  year(full_same_year$harvested[i]) <- 2020
}
gdd.date.plot <- ggplot(full_same_year, aes(x = harvested, y = GDD_at_harvest, color = year)) + geom_point() +
  labs(x = "Harvest date", y = "GDD at harvest", title = "GDD at harvest over time", 
       subtitle = "year set to same for all to see overlap\neach point is a unique pollination date/harvest date combination")
dap.date.plot <- ggplot(full_same_year, aes(x = harvested, y = DAP, color = year)) + geom_point() +
  labs(x = "Harvest date", y = "DAP at harvest", title = "DAP at harvest over time", 
       subtitle = "year set to same for all to see overlap\neach point is a unique pollination date/harvest date combination\nno DAP, harvest, or pollination dates given for 2014")

ggsave(gdd.date.plot, filename = "~/Desktop/gdd_by_date.png",  bg = "transparent", height=5, width=9)
ggsave(dap.date.plot, filename = "~/Desktop/dap_by_date.png",  bg = "transparent", height=5, width=9)

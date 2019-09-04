# Pollination date master spreadsheet script for 2019 sweet cornRNA-seq
# Jenna Hershberger
# jmh579@cornell.edu

# 8/26/2019
# Updated 9/3/2019

crosses823 <- read_excel("~/Desktop/SweetCorn/2019/crosses_2019.xlsx", na = "NA") %>% 
  mutate(pollination_date = as.Date(pollination_date, format = "%m/%d/%Y")) %>% 
  mutate(harvested = as.Date(harvested, format = "%m/%d/%Y"))

df2 <- as.data.frame(matrix(nrow = 462, ncol = 5, data = NA))
colnames(df2) <- c("plot_id", "pol_1", "pol_2", "pol_3", "pol_4")

df2[,1] <- c(paste0("19A000", 1:9), paste0("19A00", 10:99), paste0("19A0", 100:462))

crosses823$plot_number <- as.numeric(unlist(str_split(crosses823$plot_id, "A"))[c(FALSE, TRUE)])
for(i in 1:nrow(crosses823)){
  if(is.na(df2[crosses823$plot_number[i],2])){
    df2[crosses823$plot_number[i],2] <- crosses823$pollination_date[i]
  } else if(is.na(df2[crosses823$plot_number[i],3])){
    df2[crosses823$plot_number[i],3] <- crosses823$pollination_date[i]
  } else if(is.na(df2[crosses823$plot_number[i],4])){
    df2[crosses823$plot_number[i],4] <- crosses823$pollination_date[i]
  } else if(is.na(df2[crosses823$plot_number[i],5])){
    df2[crosses823$plot_number[i],5] <- crosses823$pollination_date[i]
  }
}
df2[,2] <- as.Date(df2[,2], origin = "1970-1-1")
df2[,3] <- as.Date(df2[,3], origin = "1970-1-1")
df2[,4] <- as.Date(df2[,4], origin = "1970-1-1")
df2[,5] <- as.Date(df2[,5], origin = "1970-1-1")

# Don't need to update pol.date.by.plot unless new pollinations are made
#write.csv(df2, "~/Desktop/SweetCorn/2019/pol.date.by.plot.20190903.csv", row.names = F, na = "")

harvested <- crosses823 %>% filter((harvested <= today()))
# remove plots with no pollinations
df3 <- df2 %>% filter(!plot_id %in% harvested$plot_id) %>% filter(!is.na(pol_1))

write.csv(df3, "~/Desktop/SweetCorn/2019/harvests.by.plot.20190903.csv", row.names = F, na ="")

h.plots <- harvested$plot_id
df4 <- crosses823 %>% filter(!plot_id %in% h.plots) %>% filter(pollination_date > "2019-08-10") %>% 
  arrange(pollination_date, plot_number)
write.csv(df4, "~/Desktop/SweetCorn/2019/harvests.by.date.20190903.csv", row.names = F, na ="")
  
###### other tries ####
# filter out plots with stand count of zero
stand <- read.csv("~/Desktop/db_uploads/standcounts_upload_rnaseq_20190713.csv") %>% 
  tidyr::separate(`pstand_ct_plntplot.CO_322.0000891`, into = c("pstand_ct", "tsmp"), sep = ",") %>% 
  mutate(pstand_ct = as.numeric(pstand_ct)) %>% mutate(tsmp = as.Date(tsmp)) %>% 
  dplyr::rename(plot_id = observationunit_name) 
# remove duplicated stand count for plot 231. Pick one at random for now and go back to field and count later.
stand <- stand[-220,]

try2 <- df2 %>% mutate(plot2 = plot_id) %>% 
  tidyr::separate(plot2, into = c("nineteen", "plot_number"), sep = "A") %>% 
  dplyr::select(-nineteen) %>% mutate(plot_number = as.numeric(plot_number)) %>% 
  full_join(stand) %>% 
  full_join(rnaseq_coordinates) %>% 
  ggplot(aes(x=row_number, y=col_number,
             #fill = pstand_ct)) + geom_tile()
             fill = pol_1)) + geom_tile()


#df3 <- df2 %>% gather(pol_1:pol_4, key = "pol_event", value = "pol_date") %>% tidyr::spread(pol_date, pol_event)
#
# gdd.totals <- read.csv("~/Desktop/SweetCorn/2019/GDD_totals_2019-08-23.csv")
# gdd.by.plot <- gdd.totals %>% tidyr::spread(key = pol_date, value = cross_count)
# for(col.i in 1:ncol(gdd.by.plot)){
#   if(col.i>2){
#     for(row.j in 1:nrow(gdd.by.plot)) {
#       if(!is.na(gdd.by.plot[row.j, col.i])){
#         gdd.by.plot[row.j, col.i] <- colnames(gdd.by.plot)[col.i]
#       }
#     }
#   }
# }
# 
# 
# write.csv(gdd.by.plot, "~/Desktop/gdd.by.plot.20190823.csv", row.names = F)
# 
# 
# 
# for(row.i in 1:nrow(gdd.by.plot)){
#   if(sum(is.na(gdd.by.plot[row.i,3:20]))==1){
#     gdd.by.plot[row.i,3]) <- gdd.by.plot[row.i,which(!is.na(gdd.by.plot[row.i,3:20]))+2]
#   }
#   if(sum(is.na(gdd.by.plot[row.i,3:20]))==2){
#     
#   }
#   if(sum(is.na(gdd.by.plot[row.i,3:20]))==3){
#     
#   }
#   if(sum(is.na(gdd.by.plot[row.i,3:20]))==4){
#     
#   }
# }

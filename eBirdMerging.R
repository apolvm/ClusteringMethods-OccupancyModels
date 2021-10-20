#eBird WETA and Covariates files merged code
#August 11, 2021
##libraries needed
library(dplyr)  
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(lubridate)
library(unmarked)
library(Metrics)

#Loading data
covData <- rbind( read.csv("~/RStudio/DREU2021/REU files/data/2017_UPDATED_COVS_df.csv"))
wetaData <- rbind( read.csv("~/RStudio/DREU2021/REU files/data/ebd_weta_breeding_or_zf copy.csv"))

#Filter the WETA file to only get the 2017 data
  #originally the data had 79847 observations and 19 variables
wetaData <- wetaData%>%filter(year == 2017)
  #after the filtering the data was reduced to 13010 observations and 19 variables

#check the head of both dataframes
names(wetaData)
names(covData)

#both dataframes have the following variables in common:
# "checklist_id" <- maybe
# "observer_id" <- NO
# "locality_id"  <- maybe
# "latitude" <- maybe
# longitude <- maybe
# observation_date <- no
# "day_of_year" <- no
# "time_observations_started" <- no
# "duration_minutes" <- no
# "effort_distance_km" <- no
# "number_observers" <- no

##Merging both dataset into one with all data

#dataframe size: 13010 . after: 9170 same size as covData
checklistData <- merge(wetaData, covData, by= c("checklist_id"), all.x = T)
checklistData <- checklistData%>%na.omit()
colnames(checklistData)

##Visualizing checkList data

#plotting which hour of the day is the most popular to do surveys
ggplot(checklistData) + geom_histogram(mapping = aes( x= time_observations_started.x), stat = "bin",
                                   position = "stack", binwidth = 0.05,
                                   na.rm = T 
) + labs( title= "Time of the day of surveys", x= "hour", y= "")

#type of birds recorded: only 1
ggplot(checklistData)+ geom_bar( aes(x= scientific_name), stat = "count", na.rm = T) 
+ labs( title = "Type of Birds", x= "", y = "")

#detection: yes/ no
ggplot(checklistData)+ geom_bar( aes(x= species_observed), stat = "count", na.rm = T) + 
  labs( title = "Species Observed?", x= "", y = "")
#stationary or traveling surveys
ggplot(checklistData)+ geom_bar( aes(x= protocol_type), stat = "count", na.rm = T) + 
  labs( title = "Type of Protocol", x= "", y = "")
#how long are the surveys
ggplot(checklistData%>% filter(duration_minutes.x < 75)) + geom_histogram(mapping = aes( x= duration_minutes.x), stat = "bin",
                                       position = "stack", binwidth = 1.0,
                                       na.rm = T 
) + labs( title= "How long lasted the surveys", x= "minutes", y= "")

#Using the as.Date() function, a new column Weekday was created. The as.Date() function take numeric data and represents it into calendar dates.
checklistData$weekday <- weekdays(as.Date(checklistData$observation_date.x, "%m/%d/%y"))

with(diamonds, barplot(rev(sort(table(checklistData$weekday) )[1:7]), width = 5, main= "Busiest Day of the Week", xlab= "Day of the Week")) 

ggplot(checklistData, aes(x= weekday, fill= protocol_type)) + geom_bar( stat = "count" , na.rm= T, 
                                                              position = position_dodge(width = 0.25)) + labs(
                                                                title= "Weekday and Protocol Type",
                                                                x= "Week day", y= ""
                                                              )
#separting observation date
checklistData <- checklistData%>% separate(observation_date.x, c("month", "day", "year"), "/")

#plotting which month had the most surveys recorded
ggplot(checklistData) + geom_bar(aes(x= checklistData$month),stat = "count", na.rm= T)+ 
  labs( title= "Which month had the most surveys?", x= "month")


#trying to check how many unique locality IDs there are
ggplot(checklistData)+ geom_bar( aes(x= checklist_id), stat = "count", na.rm = T) + 
  labs( title = "Locality IDs", x= "", y = "")

length(unique(checklistData$checklist_id))

length(unique(checklistData$locality_id))

#Mapping location
qmplot(longitude.x,latitude.x, data = checklistData, maptype = "toner-lite", color = I("brown"))+ 
  labs( title= "Location of surveys")

##Making "Site" column

#*spare code*
latLongDF<- cbind(checklistData$latitude.x, checklistData$longitude.x)
latLongDF <- data.frame(latLongDF)
#numbering sites
p<- unique(latLongDF)
p$site <- 1:nrow(p)
colnames(p) <- c("latitude.x", "longitude.x", "site")
#******************

#similar latitude-longitude are grouped into one site
checklistData<- merge(checklistData, p, by=c("latitude.x", "longitude.x"))
checklistData<- checklistData %>% arrange(site)
checklistData <- checklistData%>%select(site, everything())

##Splitting data into training and testing sets by unique sites

sp<- sample(c(rep(0, 0.8 * nrow(checklistData)),
              rep(1, 0.2 * nrow(checklistData))))
  
data_train<- checklistData[sp == 0, ] ##size: 7336
data_test <- checklistData[sp == 1, ] ##size: 1834

##-------Training data

data_train<- data_train[!duplicated(data_train$site), ] ##size: 4285
data_train <- data.frame(data_train)

data_train$species_observed<- as.integer(data_train$species_observed)
names(data_train) <- str_replace(names(data_train),'y.0','y.2')
names(data_train) <- str_replace(names(data_train),'sy.1n_y','syn_y')
names(data_train) <- str_replace(names(data_train),'date.1.2','date.2')
data_train<- data_train%>%relocate(y.1, .after= site)
data_train<- data_train%>%relocate(date.1, .after= aspect_mean_300)
data_train<- data_train%>%relocate(y.1, .after= site)
data_train<- data_train%>%relocate(date.2, .after= date.1)

data_train$y.1 <- NA
data_train$date.2 <-NA

head(data_train)
drops <- c("observer_id.x","sampling_event_identifier", "scientific_name", "observation_count",
           "state_code", "locality.1_id.x", "protocol_ty.1pe", "all_species_reported",
           "observation_date.1.x", "y.1ear", "day.1_of_year.x", "time_observations_started.x",
           "duration_minutes.x", "effort_distance_km.x", "number_observers.x",
           "time_observations_started.y.1", "duration_minutes.y.1",
           "effort_distance_km.y.1", "number_observers.y.1" ,
           "observation_date.1.y.1", "formatted_date.1", "latitude.y.1",
           "longitude.y.1", "locality.1_id.y", "observer_id.y.1", "as_date.1", "occupied sy.1n_y", "occupied_prob")

d1<- c("X", "occupied","sy.1n_y" )
data_train<- data_train[ , !(names(data_train) %in% d1)]

#dataframe with duplicated sites
dupSites <- checklistData[duplicated(checklistData$site), ]
dupSites<- dupSites%>%group_by(site)  %>% 
  summarise_all(funs(trimws(paste(., collapse = ''))))

##Fit occupancy models into train data to make predictions
write.csv(data_train,"~/RStudio/DREU2021/eBird-trainData.csv", row.names = FALSE)
##convert to UMF
Y1<- csvToUMF("~/RStudio/DREU2021/eBird-trainData.csv",long = FALSE, type = "unmarkedFrameOccu")

names(head(data_train))

##fit occu
obsCovs(Y1) <- scale(obsCovs(Y1))
fm1 <- occu(~1~1, Y1)
fm2 <- occu(~ date~ elevation_stdDev_150 + fall_nbr_TCA_mean_75+ fall_nbr_B4_stdDev_150+
              spring_nbr_B7_stdDev_300+ aspect_mean_300, Y1)
fm2

occuPred <- predict(fm2, type = "state", newdata= data_test, na.rm= T, inf.rm= T)


levelplot(Predicted ~ data_test$date.1 + data_test$aspect_mean_300,
          data = occuPred,
          col.regions = rev(terrain.colors(100)),
          at = seq(0,1,length.out=101))

write.csv(checklistData,"~/RStudio/DREU2021/WETAandCOVARIATESdata.csv", row.names = FALSE)

plot(testing_data$longitude.x, testing_data$latitude.x, col= occuPred, cex= 2)
qplot(testing_data$longitude.x,testing_data$latitude.x, data = occuPred, maptype = "toner-lite", 
       color = Predicted)




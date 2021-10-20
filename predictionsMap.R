##Training data and making predictions
library(dplyr)  
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(lubridate)
library(unmarked)
library(Metrics)
library(raster)
library(patchwork)
library(maps)
library(rasterVis)
library(ggthemes)

data_train <- rbind( read.csv("~/RStudio/DREU2021/REU files/data/WETAandCOVARIATESdata_v2.csv"))
#arranging sites column
data_train<- data_train %>% arrange(site)
data_train <- data_train%>%select(site, everything())

##Splitting data into training and testing sets by unique sites
sp<- sample(c(rep(0, 0.8 * nrow(data_train)),
              rep(1, 0.2 * nrow(data_train))))

training<- data_train[sp == 0, ] ##size: 7336
testing <- data_train[sp == 1, ] ##size: 1834

#cleaning and arranging training set
training<- training[!duplicated(training$site), ] ##size: 4292
training <- data.frame(training)

training$species_observed <- as.integer(training$species_observed)
names(training) <- str_replace(names(training),'species_observed','y.1')
training$y.2 <- NA
names(training) <- str_replace(names(training),'day_of_year.y','date.1')
training$date.2<- NA

names(training)
drop<- c("latitude.x", "longitude.x", "checklist_id", "observer_id.x","sampling_event_identifier",
         "scientific_name","observation_count",  "state_code","locality_id.x" ,"protocol_type",
         "all_species_reported","observation_date.x","year", "day_of_year.x","time_observations_started.x",
         "duration_minutes.x", "effort_distance_km.x","number_observers.x","X","time_observations_started.y",
         "duration_minutes.y","effort_distance_km.y","number_observers.y","observation_date.y","formatted_date",
         "latitude.y", "longitude.y", "locality_id.y","observer_id.y","as_date", "occupied","y.1_syn",
         "occupied_prob")
training<- training[ , !(names(training) %in% drop)]

training<- training%>%relocate(y.2, .after= y.1)
training<- training%>%relocate(date.1, .after= aspect_mean_1200)
names(training)

#fitting occupancy models
write.csv(training,"~/RStudio/DREU2021/eBird-trainingData.csv", row.names = FALSE)
#convert to UMF
Y1<- csvToUMF("~/RStudio/DREU2021/eBird-trainingData.csv",long = FALSE, type = "unmarkedFrameOccu")

#fit occu
obsCovs(Y1) <- scale(obsCovs(Y1))
fm1 <- occu(~1~1, Y1)
fm2 <- occu(~ date~ TCA_mean_75 + TCB_mean_75+ TCW_std_150+
              TCW_std_1200+ aspect_mean_1200, Y1)
fm2

##Making predictions
test <- rbind( read.csv("~/RStudio/DREU2021/REU files/data/OR_features.csv"))
occuPred <- predict(fm2, type = "state", newdata= test, na.rm= T, inf.rm= T)
qplot(test$longitude.x,test$latitude.x, data = occuPred, maptype = "toner-lite", 
      color = Predicted)

states<- map_data("state")
OR<- subset(states, region%in%c("oregon"))
geom_polygon(data = OR, aes(x = test$longitude.x, y = test$latitude.x), fill=NA, color = "black", size = 1)

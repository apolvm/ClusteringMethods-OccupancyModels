##Demo 1.4.2
library(unmarked)
library(Metrics)
library(caTools)
library(ggplot2)
library(ggpubr)
#########################F U N C T I O N S#####################################
subsettingData <- function(myFrame, nSites, nVisits, yDatesDF){
  nSites2 = nSites/2
  nVisits2 = nVisits*2
  #-sites-
  sitesDF <- 1:nSites2
  sitesDF <- data.frame(sitesDF)
  #-ys-
  #Copy the Y values * 2 from myFrame into yDF
  yDF <- data.frame(matrix(NA, nrow=nSites2, ncol = nVisits2))
  for(i in 1: nVisits){
    for(j in 1: nSites2){
      yDF[j, i] = yDatesDF[j, i]
    }
  }
  i=1
  k= nVisits
  while(i <= nVisits){
    for(j in 1: nSites2){
      yDF[j, k+ 1] = yDatesDF[j, i]
    }
    k= k+1
    i = i+1
  }
  
  #-dates-
  #Copy the dates values * 2 from myFrame into datesDF
  dateDF <- data.frame(matrix(NA, nrow=nSites2, ncol = nVisits2))
  index2 = nVisits +1
  k = 1
  for(i in index2: nVisits2){
    for(j in 1: nSites2){
      dateDF[j, k] = yDatesDF[j, i]
    }
    k= k+1
  }
  
  index = nVisits +1
  for(i in index: nVisits2){
    for(j in 1: nSites2){
      dateDF[j, i] = yDatesDF[j, i]
    }
  }
  
  #elev-
  elevNum <- which( colnames(myFrame) == "elev")
  elevDF <- matrix(data= NA, nrow = nSites2, ncol=1)
  i= 1
  j= i
  k= 1
  while(i <= nSites){
    meanElev <- mean(myFrame[i,elevNum], myFrame[j+1,elevNum])
    elevDF[k, ] <- meanElev
    i = i+2
    j = i
    k= k+1
  }
  elevDF <- data.frame(elevDF)
  
  #temp-
  tempNum <- which( colnames(myFrame) == "temp")
  tempDF <- matrix(data= NA, nrow = nSites2, ncol=1)
  i= 1
  j= i
  k= 1
  while(i <= nSites){
    meanTemp <- mean(myFrame[i,tempNum], myFrame[j+1,tempNum])
    tempDF[k, ] <- meanTemp
    i = i+2
    j = i
    k= k+1
  }
  tempDF <- data.frame(tempDF)
  
  cols <- c("site")
  myFrame2 <- data.frame(matrix(NA, nrow= nSites2, ncol = nVisits2))
  
  ##Creating dataframe columns
  for (i in 1:nVisits2) {
    cols[i+1] <- addString("y.",i)
  }
  
  lengthCols <- length(cols)
  cols[lengthCols +1] <- "elev"
  cols[lengthCols +2] <- "temp"
  lengthCols <- length(cols)
  j <- 1
  k<- lengthCols + 1
  
  while (j <= nVisits2) {
    cols[k] <- addString("date.",j)
    j <- j + 1
    k <- k + 1
  }
  
  ##setting columns into dataframe
  # myFrame2 <- data.frame(matrix(nrow= nSites2, ncol = length(cols)))
  
  myFrame2 <- cbind(sitesDF, yDF, elevDF, tempDF, dateDF)
  colnames(myFrame2)<- cols
  
  return(myFrame2)
}

`addString` <- function(e1, e2) {
  if (is.character(e1) | is.character(e2)) {
    paste0(e1, e2)
  } else {
    base::`+`(e1, e2)
  }
}

#Function to plot the occu rmse and det rmse based on
# row=10 col=4 matrix of rmseValues from runExperiment function
`plotExperiment` <- function(detRmseValues, occuRmseValues) {
  sites <- c("original dataset", "wrong dataset")
  yMax <- 0
  
  tempDetDF <- data.frame(colMeans(detRmseValues))
  tempOccuDF <- data.frame(colMeans(occuRmseValues))
  
  if(max(tempDetDF) > max(tempOccuDF)) {
    yMax <- max(tempDetDF)
  } else {
    yMax <- max(tempOccuDF)
  }
  
  df1 <- transform(tempDetDF, mean=rowMeans(tempDetDF), sd=apply(detRmseValues,2, sd))
  
  detPlot <- ggplot(df1, aes(sites,mean)) +
    ylim(0,yMax) +
    geom_point() +
    geom_bar(position=position_dodge(), stat="identity", colour='black') +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))
  
  df2 <- transform(tempOccuDF, mean=rowMeans(tempOccuDF), sd=apply(occuRmseValues,2, sd))
  
  occuPlot <- ggplot(df2, aes(sites,mean)) +
    ylim(0,yMax) +
    geom_point() +
    geom_bar(position=position_dodge(), stat="identity", colour='red') +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))
  
  
  
  
  ggarrange(detPlot, occuPlot, labels=c("detPlot","occuPlot"),ncol=2,nrow=1)
}

calculateRMSE <- function(occuModel,beta_int_occ, beta_temp, beta_elev,beta_int_det, beta_date ){
  occuUMF <- unname(occuModel@estimates[1]@estimates[1])
  elevUMF <- unname(occuModel@estimates[1]@estimates[2])
  tempUMF <- unname(occuModel@estimates[1]@estimates[3])
  
  detcUMF <- unname(occuModel@estimates[2]@estimates[1])
  dateUMF <- unname(occuModel@estimates[2]@estimates[2])
  
  #original estimates
  actualOccu <- c(beta_int_occ, beta_temp, beta_elev)
  actualDect <- c(beta_int_det, beta_date)
  
  ##how close are the estimates?
  
  #RMSE for occupancy
  predictedOccu <- c(occuUMF, elevUMF, tempUMF)
  occuRMSE <- rmse(actualOccu, predictedOccu)
  
  ##RMSE for detection
  predictedDect <- c(detcUMF, dateUMF)
  dectRMSE <- rmse(actualDect, predictedDect)
  
  results <- c(occuRMSE, dectRMSE)
  return(results)
}

###############################################################################
nSites= 200 #will keep changing
nVisits= 2 #will keep changing

#matrices to save the RMSE of each OM
occuEstimates <- matrix(data= NA, nrow = 10, ncol=1)
dectEstimates <- matrix(data= NA, nrow = 10, ncol=1)
occuEstimates2 <- matrix(data= NA, nrow = 10, ncol=1)
dectEstimates2 <- matrix(data= NA, nrow = 10, ncol=1)

##coefficients
beta_int_occ = 0.3
beta_date = 1.25
beta_elev = -1.6
beta_temp = -0.85
beta_int_det = 0.7

logistic = function(x) { exp(x)/(1+exp(x)) }

for(i in 1: 10){
  ##covariates
  dates = matrix(nrow=nSites,ncol=nVisits)
  for (s in 1:nSites) {
    dates[s,] = sort(round(runif(nVisits,min=1,max=365)))
  } 
  dates= scale(dates)
  elev= rnorm(nSites)
  temp= rnorm(nSites)
  
  ##probabilities
  psi = logistic(beta_int_occ + (beta_elev*elev) + (beta_temp * temp) )
  p = logistic(beta_int_det + (beta_date * dates))
  
  ##data
  Z = rbinom(n=nSites,size=1,prob=psi)
  
  ## Y ~ Bern(Zp)
  Y = matrix(NA,nrow=nSites,ncol=nVisits)
  for (s in 1:nSites){
    Y[s,] = rbinom(n=nVisits,size=1,prob=Z[s]*p[s,])
  }
  
  ##dataframe
  myFrame <- data.frame(
    "site" = matrix(1:nSites,nrow=nSites,ncol=1),
    "y" = Y,
    "elev" = elev,
    "temp" = temp,
    "date" = dates
  )
  yDatesDF <- data.frame(
    "y" = Y,
    "date" = dates
  )
  
  myFrame2<- subsettingData(myFrame, nSites, nVisits, yDatesDF)
  
  write.csv(myFrame,"/Users/anapatriciaolvera/RStudio/DREU2021/demo_wrong_dataset_experiment.csv", row.names = FALSE)
  ##convert to UMF
  Y1<- csvToUMF("/Users/anapatriciaolvera/RStudio/DREU2021/demo_wrong_dataset_experiment.csv",long = FALSE, type = "unmarkedFrameOccu")
  
  ##fit occu
  obsCovs(Y1) <- scale(obsCovs(Y1))
  fm1 <- occu(~1~1, Y1)
  fm2 <- occu(~ date~  temp + elev, Y1)
  fm2
  
  ##calculating RMSE and storing values into vector
  occu_results <- calculateRMSE(fm2, beta_int_occ, beta_temp, beta_elev, beta_int_det, beta_date )
  
  rmse_results <- c(occu_results[1], occu_results[2])
  
  ##storing values into matrices
  occuEstimates[i,] <- rmse_results[1]
  dectEstimates[i, ] <- rmse_results[2]
  
  write.csv(myFrame2,"/Users/anapatriciaolvera/RStudio/DREU2021/demo_wrong_dataset_experiment_W.csv", row.names = FALSE)
  ##convert to UMF
  Y2<- csvToUMF("/Users/anapatriciaolvera/RStudio/DREU2021/demo_wrong_dataset_experiment_W.csv",long = FALSE, type = "unmarkedFrameOccu")
  
  ##fit occu
  obsCovs(Y2) <- scale(obsCovs(Y2))
  fm1_1 <- occu(~1~1, Y2)
  fm2_1 <- occu(~ date~  temp + elev, Y2)
  fm2_1
  
  ##calculating RMSE and storing values into vector
  occu_results <- calculateRMSE(fm2_1, beta_int_occ, beta_temp, beta_elev, beta_int_det, beta_date )
  
  rmse_results <- c(occu_results[1], occu_results[2])
  
  ##storing values into matrices
  occuEstimates2[i,] <- rmse_results[1]
  dectEstimates2[i, ] <- rmse_results[2]
  
}

occuEstimatesResults <- cbind(occuEstimates, occuEstimates2)
detecEstimatesResults <- cbind(dectEstimates, dectEstimates2)

colnames(occuEstimatesResults) <- c( "original dataset", "wrong dataset")
colnames(detecEstimatesResults) <- c( "original dataset", "wrong dataset")

plotExperiment(occuEstimatesResults, detecEstimatesResults)


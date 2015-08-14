#functions for calculating common hydrologic signatures from time-series 
#of discharge (Q) and rainfall (P)


#Calculates Runoff Coefficient
runCoeff <- function(P,Q){
  
  #function returns the runoff coefficient for the supplied
  #rainfall (P) and discharge (Q)
  #both in mm 
  RC <- sum(Q, na.rm = TRUE)/sum(P, na.rm = TRUE)
  
  return(RC)
}

#Calculates the Variance Ratio of discharge to precipitation
varRatio <- function(P,Q){
  
  #function returns the ratio of precipitation variance to flow variance
  #a measure of how damped the system response is
  #var(Q)/var(P)
  VR <- var(Q, na.rm = TRUE)/var(P, na.rm = TRUE)
  
  return(VR)
}

#calculates the rising and failling limb density of the hydrograph
limbDensity <- function(Q){
  
  #returns Rising Limb Density (RLD) 
  #returns falling limb Density (FLD)
  #both calculated as #peaks divided by number of steps rising and falling, respectively.
 
  #multiply output by input timestep to give appropriate units in time.
  
  peaks <- 0 #peaks gives the number of rising limbs
  rising <- 0
  falling <- 0
  check <- 0
  end <- length(Q)-1
  for(i in 2:end){
    if(Q[i] > Q[i-1]){
      rising <- rising + 1
      if(Q[i]> Q[i+1]){
        peaks <- peaks + 1
      }
    }
    if(Q[i] < Q[i-1]){
      falling <- falling + 1
    }
  }
  
  RLD <- 0
  FLD <- 0
  if(rising > 0) RLD <- peaks/rising
  if(falling > 0) FLD <- peaks/falling
    
  return(cbind(RLD,FLD))
}

#calculates the slope of the Flow duration Curve
slopeFDC <- function(Q, ql, qu){
  
  #slope of the flow duration curve, between supplied lower (ql) and upper (qu) quantiles (fractions)
  #see equation 3 in Sawicz et al 2011 catchment classification empirical analysis of hydrologic similarity
  percentiles <- c(ql*100,qu*100)
  
  #inverse of percentiles for quantile function
  temp <- as.numeric(quantile(Q, c(1-ql, 1-qu)))
  
  #points(percentiles, temp, col = "green", pch = 19)
  
  SFDC <- (log(temp[1]) - log(temp[2]))/((qu - ql)*100)
  
  remove(temp)
  remove(percentiles)
  
    
  return(SFDC)
}

#function returns the flow duration curve for a given discharge time series (Q)
flowDurationCurve <- function(Q){
  
  percentiles <- seq(0,1, by = 0.005)
  dist <- quantile(Q,percentiles)
  p <- sort(percentiles, decreasing = TRUE)
  temp <- cbind(p*100, dist)
    
  return(temp)
}

#function calculates the ratio of Q90 to Q50 as a metric of low flow variability
Q90Q50 <- function(Q){
  
  fdc <-flowDurationCurve(Q)
  ratio <- fdc[20,2]/fdc[101,2]
  
  return(ratio)
}

slopePeakDistribution <- function(Q, ql, qu){
  
  #slope of the peak distribution curve, between supplied lower (ql) and upper (qu) quantiles
  #see section 4.2.3 in Sawicz et al 2011 catchment classification empirical analysis of hydrologic similarity
  
  len <- length(Q)
  before <- sign(Q[1:len-1]-Q[2:len])[2:(len-1)]
  after <- sign(Q[2:len]-Q[1:len-1])[1:(len-2)]
  before <- replace(before, before==0, -1)
  after <- replace(after, after==0, -1)
  comb <- (before + after)/2
  peak <- ifelse(comb > 0, 1, 0)
  peak <- c(0,peak,0)
  peak <- data.frame(Q,peak)
  peak <- peak[peak$peak == 1,1]

  temp <- as.numeric(quantile(peak, c(1-ql, 1-qu)))
  
  #points(percentiles, temp, col = "green", pch = 19)
  
  peakSlope <- (log(temp[1]) - log(temp[2]))/((qu - ql)*100)
  
  remove(peak)
  remove(temp)
  remove(comb)
  remove(after)
  remove(before)
  
    
  return(peakSlope)
}

highPulseCount <- function(Q){
  
  #returns the average number of high flow events per year
  #which is calculated as the number of events above three times median daily flow
  
  fdc <-flowDurationCurve(Q)
  thresh <- fdc[101,2]*3 #threshold 3 times the median daily flow
  eventCount <- 0
  #dayCount <- 0
  for(i in 2:length(Q)){
    if(Q[i] > thresh){
      if(Q[i-1] < thresh){
        eventCount <- eventCount + 1
      }
      #dayCount <- dayCount + 1
    }
  }
  
  return(eventCount/(length(Q)/365.25))
}

#calculates the baseflow index as the ratio of baseflow to discharge
baseFlowIndex <- function(Q){
  
  baseflow <- calcBaseLH(Q, 0.925)
    
  return(sum(baseflow)/sum(Q))  
}


#function calculates the base flow using a low pass filter, controlled by the alpha parameter
calcBaseLH <- function(discharge, a){
  
  len <- length(discharge)
  discharge <- as.numeric(discharge)
  baseflow <- vector(length = len)
  baseflow[1] <- 0.5*discharge[1]
  
  # any data missing are set to mean discharge for the time-period
  discharge[is.nan(discharge)] <- mean(discharge, na.rm = TRUE)
  
  
  for(i in 2:len){
    c1 <- a*baseflow[i-1]
    c2 <-(1-a)*0.5*(discharge[i]+discharge[i-1])
    baseflow[i] <- c1+c2
    if (baseflow[i] > discharge[i]){
      baseflow[i] <- discharge[i]
    }
  }
  
  return(baseflow)
}

#function calculates all signatures. used should be a binary vector, the length of all signtures,
#specifying whic hshould be calculated. Used if, in optimisation, signatures are not required to
#speed up calculation.

calcAllSigs <- function(Q, P, used){
  
  temp <- vector(length = 9)
  
  #runoff coefficient
  if(used[1] == 1) temp[1] <- runCoeff(P, Q)  
  
  #variance ratio
  if(used[2] == 1) temp[2] <- varRatio(P, Q)
  
  #Rising Limb Density
  if(used[3] == 1 | used[4] == 1) temp[3:4] <- limbDensity(Q)
  
  #Falling Limb Density
  #temp[4] <- limbDensity(Q)[2]
  
  #Slope Flow Duration Curve
  if(used[5] == 1) temp[5] <- slopeFDC(Q, 0.33, 0.66)
  
  #slopepeakDistribution
  if(used[6] == 1) temp[6] <- slopePeakDistribution(Q, 0.1, 0.5)  
  
  #baseflow index
  if(used[7] == 1) temp[7] <- baseFlowIndex(Q)
  
  #high pulse count
  if(used[8] == 1) temp[8] <- highPulseCount(Q)
  
  #Q90Q50
  if(used[9] == 1) temp[9] <- Q90Q50(Q)
  
  return(temp)  
}






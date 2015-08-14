#code the analysie time-varying parameter sensitivity in a series of moving windows
#output is from the hymod model.

library(ggplot2)

#source("E:/Dropbox/R_scripts/models/SA_functions.R")
source("C:/Users/Christopher/Dropbox/R_scripts/models/SAFunctions.R")

##### call in window simulation file 

    #mDir <- "E:/Dropbox/Switch-On_not_shared/HJ_Andrews/Simulations/hymod/WS01/sim1/"
    mDir <- "C:/Users/Christopher/Dropbox/Switch-On_not_shared/HJ_Andrews/Simulations/hymod/WS01/sim1/"
    #read in the behavioural model samples
    setwd(mDir)
    store <- read.table("behavioural_models.txt", header = TRUE)
    store <- data.frame(store)

##### call in window signatures file 

    #should try and call as much of this in as possible to make the script generic

    #read in the windows, and constraints
    setwd("C:/Users/Christopher/Dropbox/Switch-On_not_shared/HJ_Andrews/Data/WS1/")
    sigMat <- read.table("window_signatures.txt", header = TRUE)
    nSig <- length(sigMat[1,]) - 2
    nPar <- length(store[1,])-nSig-1
    nWin <- as.numeric(length(sigMat[,1]))
    parnames <- colnames(store)[(nSig+1):(nPar+nSig)]
    mmy <- (sigMat[,1]+sigMat[,2]+1)/2
    used <- c(1,1,0,0,0,0)

    pmin <- c(1, 1, 0, -3, -3, 0, 0, 0)
    pmax <- c(1000, 7, 20, 3, 3, 1, 7, 2000)
    pRange <- cbind(pmin,pmax)

##### specify thresholds to use when producing graphs, SA etc.

    #within the retained samples. 
    #the window size/distance from each signature is reduced to evaluate its effect on performance
    nThresh <- 6
    maxbound <- 0.3 #this is the maximum bound used in sampling as a percentage of the signature value
    thresh <- vector(length = nThresh)
    for(i in 1:nThresh) thresh[i] = maxbound - (maxbound/nThresh)*(i-1)


#two loops for summary information and plotting:
#1. loop by window, producing information for each window
#2. loop by threshold, producing information for each threshold, across windows

##### produce RSA data for each window
     
      #hydrologic signatures used in the actual calibration - see window_optimisation.R
      used <- c(1,1,0,0,1,1,0,0,0)


###### loop through each window producing RSA sensitivity output: this is sensitivity output in each window.
for(i in 1:nWin){
   
    #####define windows in each signature for sensitivity
        
        minWin <- matrix(ncol = sum(used), nrow = nThresh)
        maxWin <- matrix(ncol = sum(used), nrow = nThresh)
        
        for(j in 1:nThresh){
          #find behavioural models for the threshold across signatures used
          count <- 1
          for(k in 1:nSig){
            if(used[k] == 1){
              minWin[j, count] <- sigMat[i,(k+2)]*(1-thresh[j])
              maxWin[j, count] <- sigMat[i,(k+2)]*(1+thresh[j])
              count <- count + 1
            }
          }
        }
    
        #get data for that window
        temp <- store[store$window == mmy[i],]
    
        #remove unwanted signatures from the matrix
        #sums <- which(used == 0)
        #temp <- temp[,-sums] 
        temp <- temp[,-length(temp[1,])]
        sigUsed <- sum(used)
    
        #calculate CDFs for each threshold
        parThreshCDFs <- parCDFsens(temp, minWin, maxWin, sigUsed, pRange)
        
        #produce CDF plots for each threshold
        parPlots <- plotparCDFs(parThreshCDFs)
        
        #calculate distance from uniform prior distribution for each threshold
        pardiffs <- diffparCDFs(parThreshCDFs)
        
        #plot distances from uniform prior to compare sensitivity
        cdfdiffplot <- plotparDiffs(pardiffs)
        
        #produce par-coord plot as function of thresholds
        
            #create matrix of parameters and a column of associated thresholds
            #temp will only contain signatures used, as with minWin and maxWin
            resStore <- matrix(ncol = (nThresh + 1), nrow = resCDF)
            Threshold <- vector(length = length(temp[,1]))
            nSig <- 2
            
            for(j in 1:nThresh){
                tempThresh <- vector(length = length(temp[,1]))
                for(k in 1:nSig){
                  tt <- which(temp[,k] >= minWin[j,k] & temp[,k] <= maxWin[j,k]) 
                  tempThresh[tt] <- tempThresh[tt] + 1
                }
                indexes <- which(tempThresh == nSig)
                Threshold[indexes] <- j
            }
            ppData <- cbind(temp,Threshold)
            ppData <- ppData[,-c(1:nSig)]
            ppData <- ppData[ppData$Threshold > 0,]
               
            #produce the parallel coordinate plot for data and thresholds
            pcplot <- parCordPlot(ppData)
        
        #parameter interactions as a function of threshold too for that window  
        parInter <- parInteractions(temp, minWin, maxWin, sigUsed)
        
        #plot parameter interactions as a function of threshold applied
        parInterPlot <- plotParInteraction(parInter, 0.4)
                
}
       
#code produces evaluation of model sensitivity across different windows

#stores the number of samples stored in each window, based on the threshold
retSamp <- matrix(ncol = nWin, nrow = nThresh)

for(j in 1:nThresh){
  
  #return to top directory
  setwd(mDir)
  
  
  #create new directory to store the results of the new threshold
  dir <- toString(thresh[j])
  dir <- paste("win_dist_",gsub("[.]","point", dir), sep = "") #square brackets required to not a literal dot
  
  #### toggle this code on and off to either create the directories, or just re-populate the existing directories
  dir.create(dir)

  setwd(paste(mDir,dir, sep = ""))
  
  #####select a subset of the parameter sets that conform to new constraint
  
    #first group is the largest threshold with all samples
    if(j == 1){
      subset <- store
    }
    
    #select subgroups for each window
    if(j > 1){
        #loop through each window selecting values
        for(i in 1:nWin){
          temp <- store[store$window == mmy[i],]
          temp <- temp[abs(temp$RC-sigMat$RC[i]) < thresh[j], ] 
          temp <- temp[abs(temp$VR-sigMat$VarRatio[i]) < thresh[j], ] 
          plot(temp[,1], temp[,2], ylim = c(0,0.3), xlim = c(0,0.7))
          if(i == 1) subset <- temp
          if(i > 1) subset <- rbind(subset,temp)
        }         
    }
  
    #count and store retained samples in each window
    for(i in 1:nWin) retSamp[j,i] <- sum(subset$window == mmy[i])
  
  #create boxplots for each parameter
  for(i in (nSig+1):(nPar+nSig)){
    
    temp <- data.frame(subset[,i],subset[,length(store[1,])])
    colnames(temp) <- c(parnames[i-nSig],"Window")
    
    p <- ggplot(data = temp, aes(y = temp[,1], x = temp[,2], group = temp[,2])) + 
      geom_boxplot(fill = "white", colour = "deepskyblue4", outlier.colour="deepskyblue3") + 
      xlab("Window centre") +
      ylab(parnames[i-nSig]) + 
      #scale_x_discrete(breaks = round(seq(1963, 1989, by = 1),1)) +
      scale_x_continuous(breaks = round(seq(min(temp[,2]), max(temp[,2]), by = 3),1)) +
      #theme_bw()+
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
            axis.title.x = element_text(size = 19), axis.title.y = element_text(size = 19),
            legend.position = c(0.93,0.7),
            legend.title = element_text(size = 16, face= "bold"),
            legend.text = element_text(size = 12))
     print(p)
    
     ggsave(filename = paste("hymod__window_sensitivity_",parnames[i-nSig],".png", sep = ""), dpi = 200)
        
     print(i)
  } # end of loop through parameters















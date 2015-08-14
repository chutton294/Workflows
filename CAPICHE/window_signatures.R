library(ggplot2)
library(reshape)

#call input data
setwd("C:/Users/Christopher/Dropbox/Switch-On_not_shared/HJ_Andrews/Data/WS1")
#setwd("E:/Dropbox/Switch-On_not_shared/HJ_Andrews/Data/WS1")
data <- read.table("HJA_W1.txt", header = TRUE)

start <- 0
end <- length(data[,2])
rain <- data[start:end,3]
evap <- data[start:end,2]
temp <- data[start:end,4]

#source("E:/Dropbox/R_scripts/models/signatures.R")
source("C:/Users/Christopher/Dropbox/code_tools/Signatures/R/signatures.R")

day <- as.numeric(paste(substr(data[,1], 1, 2)))
month <- as.numeric(paste(substr(data[,1],4, 5)))
year <- as.numeric(paste(substr(data[,1],7, 10)))

aa <- cbind(day, month, year, data[,2:5])

  ##### set parameters for calculating the runoff signatures

    #burn in period (years)
    BIP <- 3
  
    #window_size (years)
    WS <- 2
    
    #signatures - see signatures.R file
    nSig <- 9
  
  ##### derive matrix of periods and runoff signatures for each period of time 
  
    #start year
    sy <- min(aa$year)
    ey <- max(aa$year)
  
    #window starts and ends
    ssy <- seq(from = sy+BIP, to = ey, by = WS)
    eey <- seq(from = sy+BIP+WS-1, to = ey, by = WS)
    ssy <- ssy[1:length(eey)]
    
    #window middle
    mmy <- (ssy+eey+1)/2
       
    #number of windows
    nWin <- length(ssy)
  
  ##### calculate signature values for each window
      
      sigMat <- matrix(nrow = nWin , ncol = 2 + nSig)
      sigMat[,1] <- ssy
      sigMat[,2] <- eey

      hydrology <- matrix(nrow = nWin, ncol = 2)

      for(i in 1:nWin){
        
          datWin <- aa[aa$year >= ssy[i] & aa$year <= eey[i], ]
          
          hydrology[i, 1] <- sum(datWin$P, na.rm = TRUE)
          hydrology[i, 2] <- sum(datWin$Q, na.rm = TRUE)
          
          sigMat[i,3:11] <- calcAllSigs(datWin$Q,datWin$P, used = rep(1,9))
        
      } # end of loop through windows

   ##### plot signatures over time
      colnames(sigMat) <- c("start", "end", "RC", "VarRatio", "RLD", "FLD", "SFDC","SPDC","BFI","HPC","Q90Q50")
      colnames(hydrology) <- c("totalP","totalQ")
      
   #ggplot of signatures over time



#for(i in 1:nSig){
        #ind <- i + 2
       # plot(mmy,sigMat[,ind], type = "l", ylab = colnames(sigMat)[ind])
      #}
    
   ##### store the signatures to file, to be called in by subsequent runs.
   write.table(sigMat, "window_signatures.txt")



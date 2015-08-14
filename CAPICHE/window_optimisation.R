library("Rcpp")
library(ggplot2)
library(reshape)
library(DEoptim)

#HYMOD model code
#setwd("E:/Dropbox/R_scripts/models/")
setwd("C:/Users/Christopher/Dropbox/code_tools/hymod/R/")
sourceCpp("hymodmain.cpp")

#signature code
source("C:/Users/Christopher/Dropbox/code_tools/Signatures/R/signatures.R")

source("c:/Users/Christopher/Dropbox/R_scripts/models/targets_for_optimisation.R")

#call input data
setwd("C:/Users/Christopher/Dropbox/Switch-On_not_shared/HJ_Andrews/Data/WS1")
data <- read.table("HJA_W1.txt", header = TRUE)


start <- 0
end <- length(data[,2])
rain <- data[start:end,3]
evap <- data[start:end,2]
temp <- data[start:end,4]

##### single run of the model to check code link is working - compiled correctly
  param <- c(3,    #ks (1/days) slow flow time constant
              1.1,    #Kq (1/days) quick flow time constant 
              8,    #DDF (mm/degC/day) degree day factor
              0,    #T2 (degC) threshold for snow melt
              0,    #T1 (degC) threshold for snow
              0.3,  #alpha (-)   reservoir split parameter
              2,  #beta (-) 0-7 distribution of soil moisture stores
              2000  #cMax (mm) 0-2000 maximum soil moisture storage
              )

  a <- Rhymod(rain, evap, temp, param)
  a <- data.frame(a)
  #windows()
  plot(data$Q[0:1000], type = "l", col = "blue")
  lines(a[1], col = "red")

  plot(data$Q[start:end],a[,1])
  abline(0,1)

##### prepare data for modelling

    day <- as.numeric(paste(substr(data[,1], 1, 2)))
    month <- as.numeric(paste(substr(data[,1],4, 5)))
    year <- as.numeric(paste(substr(data[,1],7, 10)))
    
    aa <- cbind(day, month, year, data[,2:5])

##### call in window and signature data 
    
    sigMat <- read.table("window_signatures.txt", header = TRUE)

    #total number of windows
    nWin <- as.numeric(length(sigMat[,1]))
  
    #burn in period (years) for each window during simulation
    BIP <- 1

    #signatures
    nSig <- 9

    #centre of each window
    mmy <- (sigMat[,1]+sigMat[,2]+1)/2

    #of all of the nSigs a binary vector indicating those used in calibration (1) and not(0)
    used <- c(1,1,0,0,1,1,0,0,0)
    sigUsed <- sum(used)

    #target function needs to change based on the signatures actually used.
    
    #threshold number fo signatures that need to be satisfied to store the value
    sigThresh <- 4


##### define constraints to be used within the model

    #bounds for each signatures as percentage of signature value
    bound <- 0.2
    conLower <- matrix(ncol = nSig, nrow = nWin)
    conUpper <- matrix(ncol = nSig, nrow = nWin)

    for(i in 1:nSig){
        plot(mmy, sigMat[,(i+2)], type = "l", ylim = c(0.5*min(sigMat[,(i+2)]),1.5*max(sigMat[,(i+2)])), ylab = (colnames(sigMat)[i+2]))
        conLower[,i] <- sigMat[,(i+2)]*(1-bound)
        conUpper[,i] <- sigMat[,(i+2)]*(1+bound)
        lines(mmy, conUpper[,i], col = "blue")
        lines(mmy, conLower[,i], col = "blue")
    }

    
##### loop through all windows and find a set of optimised, behavioural models for each window.

    ####define parameter ranges for the hymod model
    # can change these to "fix parameters for further simulations
    pmin <- c(1, 1, 0, -3, -3, 0, 0, 0)
    pmax <- c(1000, 7, 20, 3, 3, 1, 7, 2000)
    #pmin <- c(639, 1.19, 0, -3, -3, 0.85, 0, 0)
    #pmax <- c(639, 1.19, 20, 3, 3, 0.85, 7, 2000)
    pnames <- c("Ks","Kq","DDF","Tb","Tth","alpha","beta","Cmax")

    for(i in 1:nWin){
      
      #define constraints for the window
      min <- c(conLower[i,])
      max <- c(conUpper[i,])
      
      #get data for the window -BIP is burn in period in years
      datWin <- aa[aa$year >= (sigMat$start[i]-BIP) & aa$year <= sigMat$end[i], ]
      burnPeriod <- as.numeric(sum(datWin$year < sigMat$start[i]))
      
      #initialise targets in targets_for_optimisation - based on those signatures used.
      tt <-  data.frame(cbind(min,max,used))  
      tt <- tt[tt$used == 1,]
      #nTargets - number of targets set across each signatures - min = 2
      ttar <- setTargets(sigUsed,tt$min,tt$max,3)
      tar <- matrix(ncol = (nSig+1), nrow = length(ttar[,1]))
      tar[,nSig+1] <- ttar[,length(ttar[1,])]
      cc <- 1
      for(k in 1:nSig){
          if(used[k] == 1){
            tar[,k] <- ttar[,cc]
            cc <- cc+ 1
          }   
      }
      
      precip <- datWin$P
      temp <- datWin$T
      potevap <- datWin$PET
      runoff <- datWin$Q
         
      nPar <- length(pmin)
      nPop <- nPar*10
      len <- length(tar[,1])
      store <- matrix(nrow = 200000, ncol = (nPar+nSig+1)) 
      count <- 0
      
      ##### loop through all of the targets
      for(j in 1:len){
        
        targets <- as.numeric(c(tar[j,1:nSig]))
        
        #initial loop for model
        if(j == 1){
          
          
          #run model
          junk <- DEoptim(fn=hymod_so_multiobj, lower=pmin, upper=pmax,
                        control= DEoptim.control(VTR = 0, NP=nPop, itermax=500, reltol=0.01, steptol=100, trace=10, parallelType=0),
                        precip=precip, temp=temp, potevap=potevap, runoff=runoff, min = min, max = max, targets = targets, burnPeriod = burnPeriod, sigThresh = sigThresh, used = used)
        }
        
        
        #main loop, using previously found targets in the initialisation
        if(j > 1){
                            
          ##### initialise new search using previous seaches closest to the target
          #### alternatively use randomly generated initial values from the store
          #ind <- length(store[1,])
          #for(k in 1:count){
               #store[k, ind] <-  calDist(store[k,1:nSig],targets) 
          #}
          #sort store from smallest to largest based on distance from the current target.
          #store <- store[order(store[,ind]),]
          inpop <- matrix(nrow = nPop, ncol = nPar)
          inpop[1:nPop,1:nPar] <- as.numeric(store[1:nPop,(nSig+1):(nPar+nSig)])
          
          #assign parameter sets to init_pop based on random allocation
          inds <- runif(nPop,1,count)
          #inpop[1:nPop,1:nPar] <- as.numeric(store[inds,(nSig+1):(nPar+nSig)])
          
          #inscor- matrix(nrow = nPop, ncol = 2)
          #inscore[1:nPop,1:2] <- as.numeric(store[1:nPop,1:2])
 
          #if population is greater than stored samples, generate the rest randomly
          if(nPop > count){
            start <- count 
            end <- nPop
            for(j in start:end){
              inpop[j,] <- runif(nPar,pmin, pmax)
            }
          }
          
          junk <- DEoptim(fn=hymod_so_multiobj, lower=pmin, upper=pmax,
                          control= DEoptim.control(initialpop=inpop, VTR = 0, NP=nPop, itermax=50, reltol=0.01, steptol=100, trace=10, parallelType=0),
                          precip=precip, temp=temp, potevap=potevap, runoff=Q, min = min, max = max, targets = targets, burnPeriod = burnPeriod, sigThresh = sigThresh, used = used)
          #junk <- DEoptim(fn=hymod_so_multiobj, lower=pmin, upper=pmax,
                          #control= DEoptim.control(VTR = 0, NP=nPop, itermax=100, reltol=0.01, steptol=100, trace=10, parallelType=0),
                          #precip=precip, temp=temp, potevap=potevap, runoff=Q, min = min, max = max, targets = targets, burnPeriod = burnPeriod, sigThresh = sigThresh, used = used)
          
          remove(inpop)
          
        }   
      
        print(i)
        print(j)
      } #end of loop through all targets

      #cut_samples from store
      store <- store[1:count, ]
      window <- rep(mmy[i],count)
      store <- cbind(store,window)
      
      #add code in here that if the store is nearly full then we remove dupliactes from the store and randomise and then write over
      
      
      
      if(i == 1) master <- store
      if(i > 1) master <- rbind(master,store)
          
      
    } # end of loop through all windows
    

#set column names for master
master <- master[,-(nSig+nPar+1)]
colnames(master)[(nSig+1):(nSig+nPar)] <- pnames
colnames(master)[1:nSig] <- colnames(sigMat)[-c(1:2)]
colnames(master)[length(master[1,])] <- "window"


#remove duplacted rows from master - results from optimisation algorithm
master <- unique(master)
length(master[,1])

##### write the output from all of the simulations to file
setwd("E:/Dropbox/Switch-On_not_shared/HJ_Andrews/Simulations/hymod/WS01/sim3/")
write.table(master, "behavioural_models.txt")
    
   
    


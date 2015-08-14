setTargets <- function(nObj, min, max, nTarget){

##### function description
  
    #initialises the targets for each optimisation run 
    #nObj dimensions of target space - number of signatures
    #min vector of mimimum values for the signatures
    #max vector of maximum signatures of the signatures
    #nTarget - the number of targets required in each dimension/nObj
    
    # define number of optima to be searched for in the window space
    nOptima <- nTarget^nObj
    
    
##### define the master storage array to store all of the behavioural results #######
  
    #store <- matrix(nrow = 10000, ncol = (nPar+nObj+1)) #extra 1 for distance from points
    #count <- 1
    #total_function_evaluations <- 0
  
    ##### set up targets, a priori at edges of grid too.
        tar <- matrix(nrow = nOptima, ncol = nObj)
        cc <- 1
        
        #first, create the individual, marginal targets to be combined
        marTar <- matrix(nrow = nTarget, ncol = nObj)

        for(i in 1:nObj){
          for(j in 1:nTarget){
            marTar[j,i] <- (((j-1)/(nTarget-1))*(max[i]-min[i])) + min[i]
          }
        }
        
        #combine the individual targets together
        temp <- split(marTar, rep(1:ncol(marTar), each = nrow(marTar)))
        temp <- expand.grid(temp)

        ######start targets in the centre and work outwards...based on a normaliseed distance from centre
        centre <- (max + min)/2
        dist <- vector(length = nOptima)
        zeros <- rep(0,nObj)
        for(i in 1:nOptima){
          normDist <- abs(temp[i,]-centre)/(max-min)
          dist[i] <- calDist(zeros,normDist)
        }
        tar <- cbind(temp,dist)
        tar <- tar[order(tar[,(nObj+1)]),]
        tar <- data.frame(tar)
        
        #p <- ggplot(tar) + 
          #geom_point(data=tar,aes(x = tar[,1], y = tar[,2], colour = tar[,3])) +
          #scale_colour_gradientn(colours=rainbow(5), name = "distance from centre") +
          #geom_path(data=a, aes(x = a[,1], y = a[,2])) +
          #coord_fixed() +
          #theme(text = element_text(size=20))
        #print(p)
    
    return(tar)
}

hymod_single_objective <- function (param, precip, temp, potevap, runoff, min, max, targets, burnPeriod){
  
  simu <- Rhymod(precip, potevap, temp, param)
  
  #BIP is the burn in period to remove from the start of each run.
  if(burnPeriod > 0){
    precip <- precip[-c(1:burnPeriod)]
    simu <- simu[-c(1:burnPeriod)]
  }
  
  VR <- varRatio(precip, simu)
  RC <- runCoeff(precip, simu)
  
  if(VR < 0 | VR > 1){ lh <- 0}
  if(RC < 0 | RC > 1){ lh <- 0}
  
  lh <- calDist(cbind(RC,VR), targets)
  
  #set lh to 0 if within a tolerable distance of the target
  if(lh < 0.005) lh <- 0
  
  #identify models that fall within the high-performing space - e.g. satisfy both constraints
  if(VR < max[2] & VR > min[2]){
    if(RC < max[1] & RC > min[1]){
      
      count <<- count + 1
      #store the RC and varRatio of the problem
      store[count, 1] <<- RC
      store[count, 2] <<- VR
      #store the parameter set of the model
      temp <- length(param) + 2
      store[count, 3:temp] <<- param
      
    }
  }
  
  
  if(is.nan(lh)) lh <- 1
  
  #print(lh)
  
  #total_function_evaluations <<- total_function_evaluations + 1
  
  return(lh)
}

TUWmodel_single_objective <- function (param, precip, temp, potevap, runoff, min, max, targets, burnPeriod, area){
  
  simu <- as.numeric(TUWmodel(prec=as.numeric(precip), airt=as.numeric(temp), ep=as.numeric(potevap), area=area, 
                             param=param)$q)
  simu[is.na(simu)] <- 0
  
  #BIP is the burn in period to remove from the start of each run.
  if(burnPeriod > 0){
    precip <- precip[-c(1:burnPeriod)]
    simu <- simu[-c(1:burnPeriod)]
  }
  
  VR <- varRatio(precip, simu)
  RC <- runCoeff(precip, simu)
  
  if(VR < 0 | VR > 1){ lh <- 0}
  if(RC < 0 | RC > 1){ lh <- 0}
  
  lh <- calDist(cbind(RC,VR), targets)
  
  #set lh to 0 if within a tolerable distance of the target
  if(lh < 0.005) lh <- 0
  
  #identify models that fall within the high-performing space - e.g. satisfy both constraints
  if(VR < max[2] & VR > min[2]){
    if(RC < max[1] & RC > min[1]){
      
      count <<- count + 1
      #store the RC and varRatio of the problem
      store[count, 1] <<- RC
      store[count, 2] <<- VR
      #store the parameter set of the model
      temp <- length(param) + 2
      store[count, 3:temp] <<- param
      
    }
  }
  
  
  if(is.nan(lh)) lh <- 1
  
  #print(lh)
  
  #total_function_evaluations <<- total_function_evaluations + 1
  
  return(lh)
}

hbv_single_objective <- function (param, precip, temp, potevap, runoff, min, max, targets){
  
  simu <- Rhbv(precip, potevap, temp, param)
  
  VR <- varRatio(precip, simu)
  RC <- runCoeff(precip, simu)
  
  if(VR < 0 | VR > 1){ lh <- 0}
  if(RC < 0 | RC > 1){ lh <- 0}
  
  lh <- calDist(cbind(RC,VR), targets)
  
  #set lh to 0 if within a tolerable distance of the target
  if(lh < 0.01) lh <- 0
  
  #identify models that fall within the high-performing space - e.g. satisfy both constraints
  if(VR < max[2] & VR > min[2]){
    if(RC < max[1] & RC > min[1]){
      
      count <<- count + 1
      #store the RC and varRatio of the problem
      store[count, 1] <<- RC
      store[count, 2] <<- VR
      #store the parameter set of the model
      temp <- length(param) + 2
      store[count, 3:temp] <<- param
      
    }
  }
  
  
  if(is.nan(lh)) lh <- 1
  
  #print(lh)
  
  #total_function_evaluations <<- total_function_evaluations + 1
  
  return(lh)
}

calDist <- function(a,b){
  
  #calculates distance between vectors a and b in multiple dimensions
  c <- (a-b)^2
  c <- sum(c)^0.5
  
  return(c)
}

#single objective calibraton of hymod, where single objective is derived from multiple signatures
hymod_so_multiobj <- function (param, precip, temp, potevap, runoff, min, max, targets, burnPeriod, used, sigThresh){
  
  #min is minimum bound of target
  #max is maximum bound of target
  #targets are actual target values
  #used is a binar vector saying whether or not the signature is to be used in calibration
  #sigThresh: the number of signatures satisfied in order to store the parameter set
    
  simu <- Rhymod(precip, potevap, temp, param)
  
  #BIP is the burn in period to remove from the start of each run.
  if(burnPeriod > 0){
    precip <- precip[-c(1:burnPeriod)]
    simu <- simu[-c(1:burnPeriod)]
  }
  
  #calculate modelled signatures
  simSig <- calcAllSigs(simu,precip, used)
  nSig <- length(simSig)
  simSig[is.na(simSig)] <- 0
  
  #calculate whether models are behavioural
  behavSig <- rep(0, nSig)
  
  for(i in 1:nSig){
    if(used[i] == 1){
        if (simSig[i] >= min[i] & simSig[i] <= max[i]) behavSig[i] <- 1
    }
  }
  
  #for optimisation calculate the normalised distance so all signatures are balanced when
  #calculating distance for the single objective optimisation
  normSimSig <- vector(length = nSig)
  normSimSig <- abs(simSig-targets)/targets
  
  #get vector of behavioural signatures based on those signatures actually used in calibration
  temp <- cbind(used,normSimSig)
  temp <- temp[temp[,1] == 1, ]
  temp <- temp[,2]
   
  #of those signatures used, check distance from maximum single objective.
  #distance used need to be a normalised distance from the target...
  tt <- rep(0,length(temp))
  lh <- calDist(temp, tt)
  
  #set lh to 0 if within a tolerable distance of the target for convergence
  if(lh < 0.06) lh <- 0
  
   #identify models that fall within the high-performing space - e.g. satisfy  constraints
   #behavsig multiplied by used so that only those signatures that need to be satisfied are.
  if(sum(behavSig*used) >= sigThresh){
      
      count <<- count + 1
      #store the RC and varRatio of the problem
      store[count, 1:nSig] <<- simSig
      #store the parameter set of the model
      temp <- length(param) + nSig
      store[count, (nSig+1):temp] <<- param
  }
  
  if(is.nan(lh)) lh <- 1
  
  #print(lh)
  
  #tidy up
  remove(simu)
  remove(behavSig)
  remove(normSimSig)
  remove(temp)
  remove(tt)
  
  
  return(lh)
}















    
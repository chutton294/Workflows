library(ggplot2)
library(reshape2)
library(GGally)
library(RColorBrewer)
library(grid)

#calculates the cumulative distribution function of the supplied variable.
#returns a matrix with two columns:
#first, the value of the input variable (parvals), 
#secod, the probability of that value.
parCDF <- function(parvals, min, max, resCDF){
  
  #min <- min(parvals)
  #max <- max(parvals)
  parval <- seq(min,max, by = (max-min)/(resCDF-1))
  prob <- rep(0,1000)
  if(length(parvals)> 0) prob <- ecdf(parvals)(parval)
  p <- data.frame(cbind(parval,prob))
  
  return(p)
}

#calculates parameter sensitivity for given window threshold sensitivity for all supplied 
#parameters
parCDFsens <- function(temp, minWin, maxWin, nSig, pRange){
  
  #temp: matrix, first cols nSigs, remaining cols parameters
  #minWin/maxWin: window ranges for each signature to threshold
  #nSig: number of signatures - first nSig columns in temp
  #prange: parameter prior range - nPar rows, 2 columns: min and max of prior
  
  nPar <- length(temp[1,])- nSig
  
  resCDF <- 1000
  
  results <- vector("list", nPar)
  
  nThresh <- length(minWin[,1])
  
  for(m in 1:nPar){
    resStore <- matrix(ncol = (nThresh + 1), nrow = resCDF)
    for(j in 1:nThresh){
      #find behavioural models for the threshold across signatures used
      for(k in 1:nSig){
          if(k == 1) tt <- temp[temp[,k] >= minWin[j,k] & temp[,k] <= maxWin[j,k], ] 
          if(k > 1)  tt <- tt[tt[,k] >= minWin[j,k] & tt[,k] <= maxWin[j,k], ] 
      }
      #get behavioural paramerters, calculate and store CDF
      par <- tt[ ,nSig+m]
      #print(length(pars[,1]))
      if(j == 1){ 
        t <- parCDF(par, pRange[m,1], pRange[m,2],resCDF)
        resStore[,1] <- t[,1]
        resStore[,2] <- t[,2]
      }
      if(j > 1) resStore[,j+1] <- parCDF(par, pRange[m,1], pRange[m,2],resCDF)[,2]
    }
    a <- vector(length = nThresh + 1)
    a[1] <- colnames(temp)[m+nSig]
    a[2:(nThresh+1)] <- seq(1,nThresh, by = 1)
    colnames(resStore) <- a
    results[[m]] <- resStore
  }  
  
  return(results)
}

#plots the parameter CDFs stored in resStore
plotparCDF <- function(resStore){
  
  #prior CDF
  aa <- seq(0,1, by = (1/999))
  
  colnames <- colnames(resStore)
  resStore <- data.frame(cbind(resStore,aa))
  colnames(resStore) <- colnames
  colnames(resStore)[1] <- "par"
  colnames(resStore)[length(resStore[1,])] <- "prior"
  res <- melt(resStore, id = "par")
  pp <- res[res$variable == "prior",]
  res <- res[res$variable != "prior",]
    #colnames(tf) <- c("xval","yval")
  
  colnames(res) <- c(colnames[1], "Threshold", "value")
                   
  plotList <- list(
      geom_line(data = res, aes_string(x = colnames[1],y = "value", group = "Threshold", colour = "Threshold"), size = 1),
      geom_line(data = pp, aes_string(x = "par", y = "value"), colour = "grey53", size = 1)
  )
  
  p <- ggplot() + 
    plotList +
    xlab(colnames[1]) +
    ylab("CDF") +
    scale_x_continuous(limits = c(min(resStore[,1]),max(resStore[,1])), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_colour_brewer(palette = "Spectral") +
    theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
          axis.title.x = element_text(size = 19), axis.title.y = element_text(size = 19),
          #legend.position = c(0.93,0.7),
          legend.title = element_text(size = 16, face= "bold"),
          legend.text = element_text(size = 6),
          legend.position="none",
          plot.margin = unit(c(0.01,0.01,0.2,0.1), "cm"),
          aspect.ratio=1
          ) 
  
  remove(aa)
     
  #print(p)
  return(p)
}

#calculates the difference in CDF between the prior and posterior distribution
diffCDF <- function(resStore){
  
  aa <- seq(0,1, by = (1/999))
  
  nThresh <- length(resStore[1,])-1
  maxDiff <- vector(length = nThresh)
  for(i in 2:(nThresh + 1)) maxDiff[i-1] <- max(abs(aa-resStore[,i]))
  
  return(maxDiff)
}

diffparCDFs <- function(parCDFs){

  maxDiffMat <- matrix(ncol = length(parCDFs), nrow = length(parCDFs[[1]][1,])-1)
  colnames <- vector(length = length(parCDFs))
  for(i in 1:length(parCDFs)) {
    maxDiffMat[,i] <- diffCDF(parCDFs[[i]])
    colnames[i] <- colnames(parCDFs[[i]])[1]
  }
  colnames(maxDiffMat) <- colnames
  maxDiffMat <- cbind(seq(1,6,by = 1),maxDiffMat)
  colnames(maxDiffMat)[1] <- "Threshold"
  
  return(maxDiffMat)
}

#plots all of the CDFs for the different supplied parameters
plotparCDFs <- function(parCDFs){
  
  plots <- vector("list", length(parCDFs))
  
  for(j in 1:length(parCDFs)){
    plots[[j]] <- plotparCDF(parCDFs[[j]])
  }
  
   return(plots) 
}

#plots differences in CDFs as a function of window size.
plotparDiffs <- function(pardiffs){
  
  pd <- melt(pardiffs, id.vars = c("Threshold"))
  pd <- pd[-c(1:length(pardiffs[,1])),]
  colnames(pd) <- c("Threshold", "Parameter", "value")
  nThresh <- length(pardiffs[,1])
                        
  p <- ggplot(data = pd, aes_string(x = "Threshold", y = "value", group = "Parameter", colour = "Parameter")) +
    geom_line() +
    geom_point(size = 3) +
    xlab("Threshold") +
    ylab("CDF Difference") +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_x_continuous(limits = c(1,nThresh), breaks = 1:nThresh) + 
    #scale_colour_brewer(palette = "Spectral") +
    theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
          axis.title.x = element_text(size = 19), axis.title.y = element_text(size = 19),
          #legend.position = c(0.93,0.7),
          legend.title = element_text(size = 16, face= "bold"),
          legend.text = element_text(size = 12),
          aspect.ratio=1
    ) 
    
  return(p)
}

#produces a parallel coordinates plot of model parameters, as a function of window size
parCordPlot <-function(ppData){
  
  nPar <- length(ppData[1,])-1
  ppData <- ppData[order(ppData[,(nPar+1)]),]
  ppData[,nPar+1] <- as.numeric(ppData[,nPar+1])
  
  
  # Generate basic parallel coordinate plot
  p <- ggparcoord(data = ppData,                 
                  # Which columns to use in the plot
                  columns = 1:nPar,                 
                  # Which column to use for coloring data
                  groupColumn = nPar + 1,                 
                  # Allows order of vertical bars to be modified
                  #order = "columns",
                  order = c(1:nPar),
                  # Do not show points
                  showPoints = FALSE,                
                  # Turn on alpha blending for dense plots
                  alphaLines = 0.6,                
                  # Turn off box shading range
                  shadeBox = NULL,                
                  # Will normalize each column's values to [0, 1]
                  scale = "uniminmax" # try "std" also
  )
  
  # Start with a basic theme
  p <- p + theme_minimal()
  
  # Decrease amount of margin around x, y values
  p <- p + scale_y_continuous(expand = c(0, 0.0))
  p <- p + scale_x_discrete(expand = c(0.0, 0.0))
  
  
  # Remove axis ticks and labels
  p <- p + theme(axis.ticks = element_blank())
  p <- p + theme(axis.title = element_blank())
  p <- p + theme(axis.text.y = element_blank())
  
  mypalette<-brewer.pal(5,"Spectral")
  p <- p + scale_colour_gradientn(colours = mypalette, guide = "colorbar")
  
  # Clear axis lines
  p <- p + theme(panel.grid.minor = element_blank())
  p <- p + theme(panel.grid.major.y = element_blank())
  
  # Darken vertical lines
  p <- p + theme(panel.grid.major.x = element_line(color = "#bbbbbb"))
  
  #axis label sizes
  p <- p + theme(axis.text=element_text(size=16))
  
  # Move label to bottom
  p <- p + theme(legend.position = "bottom")
  p <- p + theme(legend.title = element_text(size = 20, face = "plain"),
                legend.text = element_text(size = 16))
  
  # Figure out y-axis range after GGally scales the data
  min_y <- min(p$data$value)
  max_y <- max(p$data$value)
  pad_y <- (max_y - min_y) * 0.1
  
  # Calculate label positions for each veritcal bar
  lab_x <- rep(1:nPar, times = 2) # 2 times, 1 for min 1 for max
  lab_y <- rep(c(min_y - pad_y, max_y + pad_y), each = nPar)
  
  # Get min and max values from original dataset
  lab_z <- c(sapply(ppData[, 1:nPar], min), sapply(ppData[, 1:nPar], max))
  
  # Convert to character for use as labels
  lab_z <- as.character(lab_z)
  lab_z <- as.character(c(pmin,pmax))
  
  # Add labels to plot
  p <- p + annotate("rect", xmin = min(lab_x)-0.5, xmax = max(lab_x)+0.5, ymin =  min(lab_y)-0.05, ymax = min(lab_y)+0.05,
                    fill = "white")
  p <- p + annotate("rect", xmin = min(lab_x)-0.5, xmax = max(lab_x)+0.5, ymin = max(lab_y)-0.05, ymax = max(lab_y)+0.05,
                    fill = "white")
  p <- p + annotate("text", x = lab_x, y = lab_y, label = lab_z, 
                    size = 5, col = "black")
  
  
  # Display parallel coordinate plot
  print(p)
  return(p)
}

#spearmans rank correlation between all parameters.
corrCoeff <- function(pars){
  
  #calculates the spearmans correlation between parameters for supplied matrix: pars,
  #with ncol: parameters; nrow: samples
  #column names must be supplied in the matrix too
  #returns matrix with par1, par2, and corr as columns.
  parnames <- colnames(pars)
  
  #get correlations between pars for the window
  z <- cor(pars, use="complete.obs", method="spearman") 
  count <- 1
  colnames(z) <-  parnames
  rownames(z) <- parnames
  z[lower.tri(z,diag=TRUE)]=NA  #Prepare to drop duplicates and meaningless information
  z=as.data.frame(as.table(z))  #Turn into a 3-column table
  z=na.omit(z)  #Get rid of the junk we flagged above
  colnames(z)[3] <- c("corr")
  colnames(z)[1:2] <- c("par1", "par2")
  return(z)
}

#evaluates parameter interactions using spearmans rank correlation
parInteractions <- function(temp, minWin, maxWin, sigUsed){
  
    
  nPar <- length(temp[1,])- sigUsed
  
  resCDF <- 1000
  
  results <- vector("list", nPar)
  
  nThresh <- length(minWin[,1])
  
    resStore <- matrix(ncol = (nThresh + 1), nrow = resCDF)
    for(j in 1:nThresh){
      #find behavioural models for the threshold across signatures used
      #loop shrinks the behavioural set (tt) by looping through the signature thresholds
      for(k in 1:sigUsed){
        if(k == 1) tt <- temp[temp[,k] >= minWin[j,k] & temp[,k] <= maxWin[j,k], ] 
        if(k > 1)  tt <- tt[tt[,k] >= minWin[j,k] & tt[,k] <= maxWin[j,k], ] 
      }
      #behavioural parameters
      par <- tt[ ,-(1:sigUsed)]
      #calculate correlations
      parCor <- corrCoeff(par)
      if(j == 1) parMaster <- parCor
      if(j > 1) parMaster <- cbind(parMaster,parCor[,3]) 
    }
      
    colnames(parMaster)[3:(nThresh+2)] <- seq(1,nThresh, by = 1)
  
  return(parMaster)    
}

#produces plot of parameter interaction
plotParInteraction <- function(parIntMat, thresh){
  
  #parIntMat <- is a nthresh + 2 column matrix. first two columns are the parameters interacting
  #remaining columns are the interaction coefficients (e.g. spearmans rank, but can be other measures
  
  #thresh: fractional value (0,1] on fraction of interactions to plot
  #interactions are ranked by the maximum value over any threshold, 
  #and the top "thresh" fraction of interactions plotted.
  
  temp <- paste(parIntMat[,1],"-", parIntMat[,2], sep = "")
  master <- cbind(parIntMat,temp)
  
  colnames(master)[length(master[1,])] <- c("Interaction")
  master <- master[,-c(1:2)]
  nInt <- as.numeric(length(master[,1]))
  nthresh <- as.numeric(length(master[1,])) -1
  
  #insert code here to only allow interactions above and below a pre-specified threshold
  keep <- vector(length = length(master[,1]))
  for(i in 1:nInt) keep[i] <- max(abs(master[i ,1:nthresh]))
  master <- cbind(master, keep)
  master <- master[order(-master$keep),]
  keep <- as.integer(nInt*thresh) #number of retained samples to plot
  master <- master[1:keep,]
  master <- master[,-length(master[1,])]
  pd <- melt(master, int.var = "Interaction")
 
  #create colours for plots from spectral palette based on interactions
  getPalette = colorRampPalette(brewer.pal(5, "Spectral")) 
  a <- getPalette(keep) #used to generate nInt colours from spectral palette.
  
  p <- ggplot(data = pd, aes(x= variable, y = value, group = Interaction, col = Interaction)) + 
    geom_point(size = 4) +
    geom_line() +
    xlab("Threshold") +
    ylab("Spearman's Rank Correlation") +  
    #scale_colour_manual(values = a) +
    #theme_bw()+
    #facet_wrap(~ sign, ncol = 1, scales = "free_y") +
    theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
          #legend.position = c(0.7,0.7), #first value is x-axis, second y-axis
          legend.title = element_text(size = 16, face= "bold"),
          legend.text = element_text(size = 12),
          strip.background = element_blank(),
          strip.text.x = element_blank()
    ) 
  
  #print(p)
  
  return(p)
}

#multiplot function to produce multiple plots in grid format
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
  


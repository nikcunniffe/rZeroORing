#
# Used for numerical integration
#
library(sfsmisc)

#
# Dispersal kernel
#
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075892
#
dispKernel <- function(r,alpha,c)
{
  retVal <- c * exp(-(r/alpha)^c) / (2 * pi * alpha^2 * gamma(2/c))
  return(retVal)
}

#
# function to calculate the log likelihood given components of R0
#
calcLogLikelihood <- function(paramVect,simulationData,maxGen)
{
  R0_11 <- exp(paramVect[1])
  R0_21 <- exp(paramVect[2])
  R0_12 <- exp(paramVect[3])
  R0_22 <- exp(paramVect[4])
  LL <- 0
  for(i in 1:nrow(simulationData))
  {
    thisRun <- simulationData[i,]
    for(g in 1:maxGen)
    {
      lastIOne <- as.numeric(thisRun[2+2*(g-1)])
      lastITwo <- as.numeric(thisRun[2+2*(g-1)+1])
      thisIOne <- as.numeric(thisRun[2+2*(g)])
      thisITwo <- as.numeric(thisRun[2+2*(g)+1])
      #print(paste("it=",i-1," g=",g," last1=",lastIOne," last2=",lastITwo," this1=",thisIOne," this2=",thisITwo,sep=""))

      ###
      ### https://www.nature.com/articles/nature04153.pdf
      ###
      ### Lloyd Smith: number of offspring from a single individual is geometric with mean R0 (i.e. prob = 1/R0)
      ###
      ### So from X independent individuals, the sum of the number of infections they cause is -ve binomial 
      ###       with prob = 1/R0 and size = X
      ###
      ### If both species are still active, can get the required contribution to the likelihood by convolution
      ###
      
      
      #      
      # Account for the iOnes in the new generation
      #
      
      #
      # both pathogens were still active in last generation => must use convolution
      #
      if(lastIOne > 0 & lastITwo > 0)
      {
        thisContrib <- 0
        for(x in 0:thisIOne)
        {
          fromOne <- x
          fromTwo <- thisIOne - x
          p1 <- dnbinom(fromOne,prob=1/(1+R0_11),size=lastIOne,log = FALSE)
          p2 <- dnbinom(fromTwo,prob=1/(1+R0_12),size=lastITwo,log = FALSE)
          thisContrib <- thisContrib + p1*p2
        }
        LL <- LL + log(thisContrib)
      }
      if(lastIOne == 0) # must've all come from type 2
      {
        LL <- LL + dnbinom(thisIOne,prob=1/(1+R0_12),size=lastITwo,log = TRUE)
      }
      if(lastITwo == 0) # must've all come from type 1
      {
        LL <- LL + dnbinom(thisIOne,prob=1/(1+R0_11),size=lastIOne,log = TRUE)
      }
      
      #      
      # Account for the iTwos in the new generation
      #
      
      #
      # both pathogens were still active in last generation => must use convolution
      #
      if(lastIOne > 0 & lastITwo > 0)
      {
        thisContrib <- 0
        for(x in 0:thisITwo)
        {
          fromOne <- x
          fromTwo <- thisITwo - x
          p1 <- dnbinom(fromOne,prob=1/(1+R0_21),size=lastIOne,log = FALSE)
          p2 <- dnbinom(fromTwo,prob=1/(1+R0_22),size=lastITwo,log = FALSE)
          thisContrib <- thisContrib + p1*p2
        }
        LL <- LL + log(thisContrib)
      }
      if(lastIOne == 0) # must've all come from type 2
      {
        LL <- LL + dnbinom(thisITwo,prob=1/(1+R0_22),size=lastITwo,log = TRUE)
      }
      if(lastITwo == 0) # must've all come from type 1
      {
        LL <- LL + dnbinom(thisITwo,prob=1/(1+R0_21),size=lastIOne,log = TRUE)
      }
    }
  }
  return(LL)
}

findRZero <- function(topLevelDir,jobName,maxR0Gen,printToScreen)
{

  #
  # Set top level directory (caching where you start, so can set it back)
  #
  startWD <- getwd()
  setwd(topLevelDir)

  #
  # Read in parameters
  #
  paramFName <- paste("Outputs\\",jobName,"_param.csv",sep="")
  myParam <- read.csv(paramFName)
  
  oRingStub <- gsub(".csv", "", myParam$xyFile)
  oRingStub <- gsub("_xy", "", oRingStub)
  
  #
  # Do the calculation to find R_0 analytically from the O-ring statistics
  #
  for(i in 1:4)
  {
    if(i == 1)
    {
      s <- "11"
      fromInf <- myParam$thetaOne
      toSusc  <- myParam$rhoOne
      fromDeath <- myParam$muOne
    }
    if(i == 2)
    {
      s <- "12"
      fromInf <- myParam$thetaOne
      toSusc  <- myParam$rhoTwo
      fromDeath <- myParam$muOne
    }
    if(i == 3)
    {
      s <- "21"
      fromInf <- myParam$thetaTwo
      toSusc  <- myParam$rhoOne
      fromDeath <- myParam$muTwo
    }
    if(i == 4)
    {
      s <- "22"
      fromInf <- myParam$thetaTwo
      toSusc  <- myParam$rhoTwo
      fromDeath <- myParam$muTwo
    }
    sFile <- paste(oRingStub,"_ORing_",s,".csv",sep="")
    oRingData <- read.csv(sFile)
    
    toInt <- (fromInf*toSusc/fromDeath) * oRingData$oRing2PiR * dispKernel(oRingData$r,myParam$dispA,myParam$dispC)
    rZeroORing <- integrate.xy(oRingData$r,toInt)
    
    if(i == 1)
    {
      R0_11 <- rZeroORing
    }
    if(i == 2)
    {
      R0_12 <- rZeroORing
    }
    if(i == 3)
    {
      R0_21 <- rZeroORing
    }
    if(i == 4)
    {
      R0_22 <- rZeroORing
    }
  }
  
  A <- matrix(c(R0_11,R0_12,R0_21,R0_22),2,2,byrow=T)
  
  eVals <- eigen(A)$values
  
  rZeroAnalytic <- max(eVals)
  
  if(printToScreen)
  {
    print(paste("Calculated:", 
                sprintf("%.3f",R0_11),
                sprintf("%.3f",R0_12),
                sprintf("%.3f",R0_21),
                sprintf("%.3f",R0_22),
                sprintf("%.3f",rZeroAnalytic)))
  }
  
  #
  # Do the estimation from this simulation's output
  #
  simDataFName <- paste("Outputs\\",jobName,".csv",sep="")
  simData <- read.csv(simDataFName)
  
  # sanity check output file is in correct format
  if(dim(simData)[1] == myParam$numIts & dim(simData)[2] == (2*(myParam$maxGen+1)+1))
  {
    # if have sufficient data
    if(maxR0Gen <= myParam$maxGen)
    {
      thisFit <- optim(c(1,1,1,1), 
                       calcLogLikelihood, 
                       control = list(fnscale = -1,maxit = 20000,trace=F), 
                       simulationData = simData,
                       maxGen = maxR0Gen)
      estimatedM <- matrix(exp(thisFit$par),nrow=2,byrow=F)
      rZero <- max(eigen(estimatedM)$values)
      if(printToScreen)
      {      
        print(paste("Estimated: ", 
                  sprintf("%.3f",estimatedM[1,1]),
                  sprintf("%.3f",estimatedM[1,2]),
                  sprintf("%.3f",estimatedM[2,1]),
                  sprintf("%.3f",estimatedM[2,2]),
                  sprintf("%.3f",rZero)))
      }
    }
  }
  setwd(startWD)
  return(c(rZeroAnalytic,rZero))
}
rm(list=ls())

#
# load spatial statistics library
#
library(spatstat)

#
# load the spline library (for calculating numerical derivatives)
#
library(pspline)

#
# for numerical integration
#
library(sfsmisc)

setwd("./Inputs/")

#
# function to lay down host locations given number of clusters as well as their width
#
makePoints <- function(xMax,yMax,nPoints,clusterNum,clusterWidth,scaleFact)
{
  allX <- numeric(nPoints)
  allY <- numeric(nPoints)
  allCentreX <- numeric(clusterNum)
  allCentreY <- numeric(clusterNum)

  goodClusters <- FALSE
  while(!goodClusters)
  {
    clusterSizes <- rpois(clusterNum-1,nPoints/clusterNum)
    if(sum(clusterSizes) < nPoints)
    {
      clusterSizes <- c(clusterSizes,nPoints-sum(clusterSizes))
      goodClusters <- TRUE
    }
  }
  thisIDX <- 1
  for(i in 1:clusterNum)
  {
    allCentreX[i] <- runif(1,max=xMax*scaleFact)
    allCentreY[i] <- runif(1,max=yMax*scaleFact)
    for(j in 1:clusterSizes[i])
    {
      v <- rnorm(2,sd=clusterWidth*scaleFact)
      allX[thisIDX] <- allCentreX[i] + v[1]
      allY[thisIDX] <- allCentreY[i] + v[2]
      thisIDX <- thisIDX + 1
    }
  }
  return(list(allX,allY,clusterSizes,allCentreX,allCentreY))
}

#
# These are parameters used in building the landscapes
#
xMaxStart <- 1
yMaxStart <- 1

nSpOne <- 800
nSpTwo <- 1200

#
# How many landscapes to build
#
nLS <- 2

for(lsNum in 1:nLS)
{
  clusterNumSpOne <- 1+ceiling(9*runif(1))
  clusterNumSpTwo <- 1+ceiling(9*runif(1))
  clusterWidthSpOne <- 0.025+0.175*runif(1)
  clusterWidthSpTwo <- 0.025+0.175*runif(1)

  if(runif(1) < 0.5)
  {
    scaleFact <- runif(1,min=0.1,max=1.0)
  }
  else
  {
    scaleFact <- runif(1,min=1.0,max=5.0)
  }

  lsName <- paste("ls_",lsNum,sep="")

  sMeta <- paste(lsName,"_meta.csv",sep="")
  metaDF <- data.frame(N1=nSpOne,
                       c1=clusterNumSpOne,
                       w1=clusterWidthSpOne,
                       N2=nSpTwo,
                       c2=clusterNumSpTwo,
                       w2=clusterWidthSpTwo)
  
  write.table(metaDF,sMeta,row.names=F,quote=F,sep=",")

  t <- makePoints(xMaxStart,
                  yMaxStart,
                  nSpOne,
                  clusterNumSpOne,
                  clusterWidthSpOne,
                  scaleFact)

  sTitle <- sprintf("%d Type 1 (c1=%d,w1=%.2f); %d Type 2 (c2=%d,w2=%.2f); scale=%.4f", 
                    nSpOne, clusterNumSpOne, clusterWidthSpOne, 
                    nSpTwo, clusterNumSpTwo, clusterWidthSpTwo, 
                    scaleFact)

  globMax <- max(max(t[[1]]),max(t[[2]]))
  globMin <- min(min(t[[1]]),min(t[[2]]))

  typeOneCentreX <- t[[4]]
  typeOneCentreY <- t[[5]]

  ppSpOne <- ppp(t[[1]],t[[2]], c(globMin,globMax), c(globMin,globMax))

  t <- makePoints(xMaxStart,
                  yMaxStart,
                  nSpTwo,
                  clusterNumSpTwo,
                  clusterWidthSpTwo,
                  scaleFact)

  typeTwoCentreX <- t[[4]]
  typeTwoCentreY <- t[[5]]

  globMax <- max(max(t[[1]]),max(t[[2]]),globMax)
  globMin <- min(min(t[[1]]),min(t[[2]]),globMin)

  ppSpTwo <- ppp(t[[1]],t[[2]], c(globMin,globMax), c(globMin,globMax))

  ppSpOne <- ppSpOne %mark% factor(rep(1,npoints(ppSpOne)))
  ppSpTwo <- ppSpTwo %mark% factor(rep(2,npoints(ppSpTwo)))

  pp <- superimpose(ppSpOne,ppSpTwo)

  sFile <- paste(lsName,"_pic.png",sep="")
  png(sFile,width=8,height=8, units='in',res=900)
  par(mfrow=c(1,1))

  # plot in white to get the right bounds to the figure...
  plot(pp$x,pp$y,col="white",pch=20,asp=1,xlab="",ylab="",
       main=sTitle,lwd=2,bty="n",xaxt="n",yaxt="n")
  points(pp$x[pp$marks==2],pp$y[pp$marks==2],col="blue",
         pch=1,asp=1,xlab="",ylab="",main=sTitle,lwd=2,bty="n",xaxt="n",yaxt="n")#)
  points(pp$x[pp$marks==1],pp$y[pp$marks==1],col="red",
         pch=4,lwd=2)
  par(lwd=2)
  legend("bottom",
         c("Type 1 host", "Type 2 host"),
         pch=c(4,1),
         lwd=c(NA,NA),
         lty=c(1),
         col=c("red","blue"),
         inset=c(0,-0.1),
         xpd=TRUE,
         ncol = 2,
         bty="n",
         cex=1.5,
         pt.cex=c(1,1))

  par(xpd=TRUE)

  yPos <- min(pp$y) + 0.1
  xPosStart <- max(pp$x) - 0.2

  lines(c(xPosStart,xPosStart+0.25),c(yPos,yPos),col="black",lwd=4)
  text(xPosStart+0.125,yPos+0.05,"0.25 km",cex=1.5)

  dev.off()

  xMax <- (max(pp$x) - min(pp$x))
  yMax <- (max(pp$y) - min(pp$y))

  fv <- Kest(pp,correction="none",rmax=sqrt(xMax*xMax+yMax*yMax))

  fv_11 <- Kcross(pp, 1, 1, r=fv$r, correction="none")
  fv_12 <- Kcross(pp, 1, 2, r=fv$r, correction="none")
  fv_21 <- Kcross(pp, 2, 1, r=fv$r, correction="none")
  fv_22 <- Kcross(pp, 2, 2, r=fv$r, correction="none")

  globalMax <- 0

  # 
  # Do two passes as a simple mechanism to allow all graphs to have same y-axis maximum
  #
  for(z in 1:2)
  {
    for(i in 1:4)
    {
      if(i == 1)
      {
        fv <- fv_11
        thisDen <- ppSpOne$n/((globMax-globMin)^2)
        s <- "11"
        thisCol <- "red"
        thisReqArea <- ppSpOne$n
      }
      if(i == 2)
      {
        fv <- fv_12
        thisDen <- ppSpTwo$n/((globMax-globMin)^2)
        s <- "12"
        thisCol <- "purple"
        thisReqArea <- ppSpTwo$n
      }
      if(i == 3)
      {
        fv <- fv_21
        thisDen <- ppSpOne$n/((globMax-globMin)^2)
        s <- "21"
        thisCol <- "purple"
        thisReqArea <- ppSpOne$n
      }
      if(i == 4)
      {
        fv <- fv_22
        thisDen <- ppSpTwo$n/((globMax-globMin)^2)
        s <- "22"
        thisCol <- "blue"
        thisReqArea <- ppSpTwo$n
      }
      sFile <- paste(lsName,"_ORing_",s,".png",sep="")
      png(sFile,width=8,height=8, units='in',res=900)
      par(las=1)
      par(mar=c(4.1,6.1,4.1,2.1))


      sFile <- paste(lsName,"_ORing_",s,".csv",sep="")
      fittedSpline <- sm.spline(fv$r, fv$un)
      numDer <- predict(fittedSpline,fv$r,1)
      oRingSmooth <- thisDen*(1/(2*pi*(fv$r)))*numDer
      oRingSmooth[1] <- 0
      oRingSmooth[oRingSmooth<0] <- 0

      dKbydR <- diff(fv$un)/diff(fv$r)
      oRing <- thisDen*(1/(2*pi*fv$r[-1]))*dKbydR

      area <- integrate.xy(fv$r,2*pi*fv$r*oRingSmooth) 

      allORing <- data.frame(r=fv$r,oRing=oRingSmooth,oRing2PiR=2*pi*fv$r*oRingSmooth)

      write.table(allORing,sFile,row.names=F,quote=F,sep=",")

      sTitle <- sprintf("ls=%s type=%s (area=%.3f)", lsName, s, area)

      print(sTitle)

      fromType <- strsplit(s,"")[[1]][1]
      toType <- strsplit(s,"")[[1]][2]
      sTitle <- sprintf("Density of type %s at radius r of type %s", toType, fromType)

      toPlot <- oRingSmooth*2*pi*fv$r
      thisMax <- max(toPlot)
      if(thisMax > globalMax)
      {
        globalMax <- thisMax
      }
      if(z == 1)
      {
        # Note the first past is just used to get y-axes of all graphs to be the same (via globalMax)
      }else{
        plot(fv$r,toPlot,ty="l",ylim=c(0,globalMax),
             xlab="Radius, r (km)",ylab="",main=sTitle, 
             cex.lab=1.5,cex.axis=1.5,lwd=3,xaxs="i",yaxs="i",col=thisCol,cex.main=1.5)
        title(ylab=expression(paste("Host density, ",2* pi * r *O(r))),line=4,cex.lab=1.5)
      }
      dev.off()
    }
  }

  # write out the locations
  df <- data.frame(x=pp$x,y=pp$y,t=pp$marks)
  sFile <- paste(lsName,"_xy.csv",sep="")
  write.table(df,sFile,row.names=F,quote=F,sep=",")
}

setwd("../")
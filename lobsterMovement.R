#updated Dec 2, 2020
library(raster)
library(ggplot2)
library(reshape2)
library(geosphere)

# A function to get the distance of each lobster from the trap
# By default the trap is at (5,5)
# The function needs the x and y coordinate of the lobster
# adding trap location makes it generic
distanceToTrapCalculator<- function(xLobster, yLobster,xtrap=5, ytrap=5){
    distanceToTrap<- sqrt((xLobster - xtrap)^2 + (yLobster -ytrap)^2)
    return(distanceToTrap)
}

#when a lobster is outside the influence of bait, we need a function
#to randomly move the lobster for a fixed length = dStep 
# The function needs the current x and y coordinate of lobster + dStep (here 5)

randomMove<- function(xLobster, yLobster, dStep){
  
    randomAngle<- runif(n=1, min = 0, max=360) #selects a random angle of direction 
    xNew<- dStep * sin(randomAngle * pi / 180) + xLobster #moves dStep based on the angle
    yNew<- dStep * cos(randomAngle * pi / 180) + yLobster
    
    return( list(EASTING = xNew, NORTHING = yNew) )
}

# Another function is needed to  move the lobster for a fixed length = dStep toward the trap (+some randomness) 
#when it is within bait influence
# The function needs the current x and y coordinate of lobster + dStep + distance of the lobster from the bait 
#to determine if the condition for directionalMove is met
   
directionalMove<- function(xLobster, yLobster, dStep, distanceToTrap,
                           xtrap=5, ytrap=5,radius_of_influence=15,ZoI){

  thetaT = bearing(c(xLobster,yLobster),c(xtrap,ytrap))#calculates the angle to the trap
  
  # Calculating tethaR
    #ZoI is adjusted each time step for decay (shrinkage) of radius of influence (equation 3 in A)
    b = 1 + 0.9 * (distanceToTrap - ZoI) / radius_of_influence
    thetaR = -180:180
    P = 1/(180^b) * abs(thetaR) ^ b
    Prtheta_R = (1-P) / sum(1-P)
    theta_r = sample(thetaR,size = 1, prob=Prtheta_R)  
  
  
   tetha      <- thetaT + theta_r
  
  xNew   <- dStep * sin(tetha * pi / 180) + xLobster
  yNew   <- dStep * cos(tetha * pi / 180) + yLobster
  
  return( list(EASTING = xNew, NORTHING = yNew) )
  
}
# A function to get the coordinate at time t, update the coordinate, 
#and return the coordinate at time t + 1
updateGrid    = function(lobsterCoordinates,trapCoordinates, radius_of_influence=15, dstep = 5, currentZoI){
  
  # Takes the argument lobsterCoordinates as the x and y coordinates and updates it accordingly
  numberOfLobsters <- nrow(lobsterCoordinates)
  xNew <- vector(mode = 'numeric', length = numberOfLobsters)
  yNew <- vector(mode = 'numeric', length = numberOfLobsters)
  
  # Takes the argument trapCoordinates as the x and y coordinates
  numberOfTraps <- nrow(trapCoordinates)
  xTrap <- vector(mode = 'numeric', length = numberOfTraps)
  yTrap <- vector(mode = 'numeric', length = numberOfTraps)
  
  # Loops over each lobster and update their coordinate accordingly
  for( lobsterIndex in 1:numberOfLobsters ){
    
    xOld <- lobsterCoordinates[lobsterIndex,1]
    yOld <- lobsterCoordinates[lobsterIndex,2]
    distanceToTrap <- distanceToTrapCalculator(xLobster = xOld,yLobster = yOld,xtrap=xTrap, ytrap = yTrap)
    
    if( distanceToTrap > radius_of_influence){
      temp <- randomMove(xLobster = xOld , yLobster = yOld , dStep = dstep)
      xNew[lobsterIndex] <- temp$EASTING
      yNew[lobsterIndex] <- temp$NORTHING
    }else{
      temp <- directionalMove(xLobster = xOld , yLobster = yOld , distanceToTrap = distanceToTrap, radius_of_influence = radius_of_influence, dStep = dstep, ZoI = currentZoI)
      xNew[lobsterIndex] <- temp$EASTING
      yNew[lobsterIndex] <- temp$NORTHING
    }
  }
  
  updatedGrid <- data.frame(EASTING = xNew, NORTHING = yNew)
  return(updatedGrid)
  
}

#### initial coordinate of lobsters is simulated (written by Adam Cook)
rpoisD<-function (n, lambda,D=1) {
  if (D==1){
    rpois(n, lambda)
  }  else {
    sz = lambda^2/(D*lambda-lambda) #this puts the overdispersion in terms of lambda or the mean and gives you appropriate size for negbin
    rnbinom(n, size=sz, mu=lambda)
  }
}



#AMC 
catchability<- function(q0= 0.5,qmin=0, saturationTreshhold=5, Ct=0 ){
   r = (log(0.01) - log(q0 - qmin))/- saturationTreshhold
    qo = (q0-qmin) / exp(r*Ct) + qmin
   return(qo)
 }

#AMC Need to id if trap in path then lobster is caught


trapInPath = function(loc1, loc2, trap_loc,how_close=0.1){
  x = seq(loc1$EASTING,loc2$EASTING,length.out = 10)
  y = seq(loc1$NORTHING,loc2$NORTHING,length.out = 10)
  path = data.frame(EASTING=x, NORTHING=y)
  dB = distanceToTrapCalculator(loc1,trap_loc)
  cl_trap = trap_loc[dB[3],]
  ds = unlist(apply(path,1,unitD,x=cl_trap))
  out = list(NA)
  if(any(ds<how_close)) {
    i= min(which(ds<how_close))
    path = path[i,]
    out[[1]] = cl_trap
    out[[2]] = path
  }
  return(out)
}


#density per grid 100 grids
# each element in the vector represents the starting value in each grid
initialGrid = rpoisD(n=100,lambda=.1, D = 3)
LobsterStart = data.frame(EASTING = rep(1:10,times=10), NORTHING = rep(1:10,each=10), Lobs = initialGrid)

LobsterStart <- subset(LobsterStart,Lobs>0)
replicateCoordinates <- function(d){ rep(d[1:2], d[3]) } #replicates coordinates for grids with more than 1 lobster
tt <- unlist( apply(X = LobsterStart, MARGIN = 1, FUN = replicateCoordinates) )
tt<- matrix(tt, ncol = 2, byrow = TRUE)
colnames(tt)<- c("EASTING","NORTHING")
initialxyCoordinate  = as.data.frame(tt)
coordinatesOverTime <- list()
coordinatesOverTime[[1]] <- initialxyCoordinate 
currentZoI<- 1
s =  0.993


for(t in 2:100){
  currentZoI<- currentZoI * s
  coordinatesOverTime[[t]] <- updateGrid( lobsterCoordinates = coordinatesOverTime[[t-1]], currentZoI = currentZoI)
  par( mfrow=c(1,2) ) # create a plot with 1 row and 2 columns to show plots side by side
  plot( coordinatesOverTime[[t-1]], xlim = c(0,10), ylim = c(0,10), main = paste0('Time = ', t-1) )
  points(x = 5, y = 5, col = 'red')
  plot( coordinatesOverTime[[t]],   xlim = c(0,10), ylim = c(0,10), main = paste0('Time = ', t) )
  points(x = 5, y = 5, col = 'red')
  
}




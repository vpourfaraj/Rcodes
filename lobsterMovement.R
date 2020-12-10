#updated Dec 2, 2020
library(raster)
library(ggplot2)
library(reshape2)
library(geosphere)

# A function to get the distance of each lobster from the trap
# By default the trap is at (5,5)
# The function needs the x and y coordinate of the lobster
# adding trap location makes it generic
distanceToTrapCalculator<- function(Lobster,trap = x(5,5)){
    xLobster = Lobster[1]
    yLobster = Lobster[2]
    xtrap = trap[1]
    ytrap = trap[2]
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


#### initial coordinate of lobsters is simulated (written by Adam Cook)
rpoisD<-function (n, lambda,D=1) {
  if (D==1){
    rpois(n, lambda)
  }  else {
    sz = lambda^2/(D*lambda-lambda) #this puts the overdispersion in terms of lambda or the mean and gives you appropriate size for negbin
    rnbinom(n, size=sz, mu=lambda)
  }
}


dispersion <- function(x) {
          var(x) / mean(x)
}


#AMC need to make general if we have multiple traps, so need to id closest trap
distanceClosestTrap <- function(lob_loc, trap_loc){
      #trap_loc can be a matrix (n x 2)
      #lob_loc is an individual lobster
  ds = unlist(apply(trap_loc,1,distanceToTrapCalculator,lob_loc))
  dmin = which.min(ds)
  return(c(distance=ds[dmin], trapid =dmin))
}




#AMC 
catchability<- function(q0= 0.5,qmin=0, saturationThreshold=5, Ct=0 ){
   r = (log(0.01) - log(q0 - qmin))/- saturationThreshold
    qo = (q0-qmin) / exp(r*Ct) + qmin
   return(qo)
 }

#AMC Need to see if trap in path then lobster is caught
trapInPath = function(loc1, loc2, trap_loc,how_close=0.1){
  x = seq(loc1[1],loc2[1],length.out = 10)
  y = seq(loc1[2],loc2[2],length.out = 10)
  path = data.frame(EASTING=x, NORTHING=y)
  ds = unlist(apply(path,1,distanceToTrapCalculator,trap=trap_loc))
  if(any(ds<how_close)) {
    i= min(which(ds<how_close))
    path = c(path[i,1],path[i,2])
    trapped = 1
    } else {
     path = loc2
     trapped = 0
    }
  return(c(path, trapped))
}



# A function to get the coordinate at time t, update the coordinate, 
#and return the coordinate at time t + 1
updateGrid    = function(lobsterCoordinates,trapCoordinates, trapCatch, radius_of_influence=15, dstep = 5, currentZoI, how_close=0.1, q0=.5, qmin=0, saturationThreshold=5){
  
  # Takes the argument lobsterCoordinates as the x and y coordinates and updates it accordingly
  numberOfLobsters <- nrow(lobsterCoordinates)
  xNew <- vector(mode = 'numeric', length = numberOfLobsters)
  yNew <- vector(mode = 'numeric', length = numberOfLobsters)
  trappedLobster<- vector(mode = 'numeric', length = numberOfLobsters)

  # Loops over each lobster and update their coordinate accordingly
  for( lobsterIndex in 1:numberOfLobsters ){
    
    xOld <- lobsterCoordinates[lobsterIndex,1]
    yOld <- lobsterCoordinates[lobsterIndex,2]
    trapped <- lobsterCoordinates[lobsterIndex,3]
    if(trapped==1){
      #repeats trap loc for a trapped lobster
      xNew[lobsterIndex] = xOld
      yNew[lobsterIndex] = yOld
      trappedLobster[lobsterIndex] = trapped
      next()
    }
    #AMC this needs to be made generic so that it is distance to closest trap
   # distanceToTrap <- distanceToTrapCalculator(Lobster = c(xOld, yOld),trap=c(xTrap, yTrap))
    distanceToTrap = distanceClosestTrap(lob_loc = c(xOld,yOld), trap_loc = trapCoordinates[,c(1,2)] )

    if( distanceToTrap[1] > radius_of_influence){
      temp <- randomMove(xLobster = xOld , yLobster = yOld , dStep = dstep)
      xNew[lobsterIndex] <- temp$EASTING
      yNew[lobsterIndex] <- temp$NORTHING
      trappedLobster[lobsterIndex] = 0
    }else{
      temp <- directionalMove(xLobster = xOld , yLobster = yOld , distanceToTrap = distanceToTrap[1], radius_of_influence = radius_of_influence, dStep = dstep, ZoI = currentZoI)
      xNew[lobsterIndex] <- temp$EASTING
      yNew[lobsterIndex] <- temp$NORTHING
    
      #AMC now need to check if trap is within path of lobster to be caught, i.e. does the lobster actually interact with the trap along its path from p1 to p2
      trappedQ = trapInPath(loc1 = c(xOld,yOld), loc2 = c(xNew,yNew), trap_loc = trapCoordinates[distanceToTrap[2],],how_close=how_close)
      if(trappedQ[3]==1) {
        
        #this means the lobster is close enough to be trapped and we need to apply the catchability criteria
        #But we need to know how many lobsters in the trap at this time using the trapCatch vector (which should be indexed the same as the trapCoordinates)
        #catchability returns a prob or being caught in the trap given the current catch  
          pC = catchability(q0= q0,qmin=qmin, saturationThreshold=saturationThreshold, Ct=trapCatch[distanceToTrap[2]] )
          caught = rbinom(n=1,size=1,prob=pC)
          if(caught==1){
            trapCatch[distanceToTrap[2]] =trapCatch[distanceToTrap[2]]+1
            xNew[lobsterIndex] <- trapCoordinates[distanceToTrap[2],1] 
            yNew[lobsterIndex] <- trapCoordinates[distanceToTrap[2],2] 
            trappedLobster[lobsterIndex] = 1
            }
          } 
        }
      }

  
  
  updatedGrid <- data.frame(EASTING = xNew, NORTHING = yNew, trapped = trappedLobster)
  return(list(updatedGrid, trapCatch))
  
}

#density per grid 100 grids
# each element in the vector represents the starting value in each grid

#One iteration 
initialGrid = rpoisD(n=100,lambda=.1, D = 3)
LobsterStart = data.frame(EASTING = rep(1:10,times=10), NORTHING = rep(1:10,each=10), Lobs = initialGrid)

LobsterStart <- subset(LobsterStart,Lobs>0)
replicateCoordinates <- function(d){ rep(d[1:2], d[3]) } #replicates coordinates for grids with more than 1 lobster
tt <- unlist( apply(X = LobsterStart, MARGIN = 1, FUN = replicateCoordinates) )
tt<- matrix(tt, ncol = 2, byrow = TRUE)
colnames(tt)<- c("EASTING","NORTHING")
initialxyCoordinate  = as.data.frame(tt)
initialxyCoordinate$trapped = 0 # this will update as a lobster gets caught and we don't need to update movements
coordinatesOverTime <- list()
coordinatesOverTime[[1]] <- initialxyCoordinate 
currentZoI<- 1
s =  0.993
ntraps = 1
trapCoordinates = data.frame(EASTING=5,NORTHING=5)
trapCatch = list()
trapCatch[[1]] = rep(0,length=ntraps)

for(t in 2:100){
  if(t>2) currentZoI<- currentZoI * s
  ko = updateGrid( lobsterCoordinates = coordinatesOverTime[[t-1]], trapCoordinates=trapCoordinates, trapCatch=trapCatch[[t-1]], currentZoI = currentZoI,saturationThreshold=5,how_close=1,dstep=5)
  coordinatesOverTime[[t]] <- ko[[1]]
  trapCatch[[t]] <- ko[[2]]
  par( mfrow=c(1,2) ) # create a plot with 1 row and 2 columns to show plots side by side
  plot( coordinatesOverTime[[t-1]], xlim = c(0,10), ylim = c(0,10), main = paste0('Time = ', t-1) )
  points(x = 5, y = 5, col = 'red')
  plot( coordinatesOverTime[[t]],   xlim = c(0,10), ylim = c(0,10), main = paste0('Time = ', t) )
  points(x = 5, y = 5, col = 'red')
  
}

#lobsters
outmove = do.call(rbind,coordinatesOverTime)
outmove$T = rep(0:99, each=nrow(tt))
outmove$I = rep(1:nrow(tt), times=100)

#traps
outtraps = as.data.frame(do.call(rbind, trapCatch))
outtraps$trapno = rep(ntraps,times=100) #if >1 trap this needs to be 1:ntraps

with(subset(outtraps, trapno==1),plot(1:100, V1, type='l',ylab='Catch', xlab='Time'))

#One lobster
  par( mfrow=c(2,2) ) # create a plot with 1 row and 2 columns to show plots side by side

with(subset(outmove, I==1),plot(EASTING, NORTHING, type='l'))
points(trapCoordinates, pch=16, col='red')
with(subset(outmove, I==1), text(x=EASTING[c(1,100)], y=NORTHING[c(1,100)], c('Start','End')))


with(subset(outmove, I==3),plot(EASTING, NORTHING, type='l'))
points(trapCoordinates, pch=16, col='red')
with(subset(outmove, I==3), text(x=EASTING[c(1,100)], y=NORTHING[c(1,100)], c('Start','End')))

with(subset(outmove, I==5),plot(EASTING, NORTHING, type='l'))
points(trapCoordinates, pch=16, col='red')
with(subset(outmove, I==5), text(x=EASTING[c(1,100)], y=NORTHING[c(1,100)], c('Start','End')))

with(subset(outmove, I==7),plot(EASTING, NORTHING, type='l'))
points(trapCoordinates, pch=16, col='red')
with(subset(outmove, I==7), text(x=EASTING[c(1,100)], y=NORTHING[c(1,100)], c('Start','End')))





SimulateLobsterMovement <- function(p=p,plot=F) {
with(p,{
  

initialGrid = rpoisD(n=ngrids,lambda=initlambda, D = initD)
LobsterStart = data.frame(EASTING = rep(1:ncolgrids,times=nrowgrids), NORTHING = rep(1:nrowgrids,each=ncolgrids), Lobs = initialGrid)

LobsterStart <- subset(LobsterStart,Lobs>0)
replicateCoordinates <- function(d){ rep(d[1:2], d[3]) } #replicates coordinates for grids with more than 1 lobster
tt <- unlist( apply(X = LobsterStart, MARGIN = 1, FUN = replicateCoordinates) )
tt<- matrix(tt, ncol = 2, byrow = TRUE)
colnames(tt)<- c("EASTING","NORTHING")
initialxyCoordinate  = as.data.frame(tt)
initialxyCoordinate$trapped = 0 # this will update as a lobster gets caught and we don't need to update movements
coordinatesOverTime <- list()
coordinatesOverTime[[1]] <- initialxyCoordinate 
currentZoI<- currentZoIInit
s =  smult
ntraps = ntrapsstart
trapCoordinates = data.frame(EASTING=trapEastStart,NORTHING=trapNorthStart)
trapCatch = list()
trapCatch[[1]] = rep(0,length=ntraps)

for(t in 2:niter){
  if(t>2) currentZoI<- currentZoI * s
  ko = updateGrid( lobsterCoordinates = coordinatesOverTime[[t-1]], trapCoordinates=trapCoordinates, trapCatch=trapCatch[[t-1]], currentZoI = currentZoI,saturationThreshold=saturationThresholdStart,how_close=how_closeStart,dstep=dstepstart)
  coordinatesOverTime[[t]] <- ko[[1]]
  trapCatch[[t]] <- ko[[2]]
  if(plot){
  par( mfrow=c(1,2) ) # create a plot with 1 row and 2 columns to show plots side by side
  plot( coordinatesOverTime[[t-1]], xlim = c(0,10), ylim = c(0,10), main = paste0('Time = ', t-1) )
  points(x = 5, y = 5, col = 'red')
  plot( coordinatesOverTime[[t]],   xlim = c(0,10), ylim = c(0,10), main = paste0('Time = ', t) )
  points(x = 5, y = 5, col = 'red')
  }
}

#lobsters
outmove = do.call(rbind,coordinatesOverTime)
outmove$T = rep(0:(niter-1), each=nrow(tt))
outmove$I = rep(1:nrow(tt), times=niter)

#traps
outtraps = as.data.frame(do.call(rbind, trapCatch))
if(ntraps==1) outtraps$trapno = rep(ntraps,times=niter) #if >1 trap this needs to be 1:ntraps
if(ntraps>1) outtraps$trapno = rep(1:ntraps,times=niter) #if >1 trap this needs to be 1:ntraps
outputs = list()
outputs$traps = outtraps
outputs$lobsters = outmove
return(outputs)
  })
}


#Lets run Multiple iterations

#initalize a parameter file to pass info into the code and then put all into a function

p = list()
p$nrowgrids = 10
p$ncolgrids = 10
p$ngrids=p$nrowgrids * p$ncolgrids
p$initlambda=.1
p$initD = 3
p$smult = 0.993
p$currentZoIInit = 1

p$trapEastStart = 5
p$trapNorthStart = 5
p$ntrapsstart = length(p$trapEastStart)

p$saturationThresholdStart = 5
p$how_closeStart = 1
p$dstepstart = 5 

p$niter =100




#run the model
    a = SimulateLobsterMovement(p=p)

      plot(1:p$niter,a$traps[,1],xlab='Time',ylab='N Caught',ylim=c(0,15))

#lets change a parameter
  p$saturationThresholdStart=10

# rerun
    b = SimulateLobsterMovement(p=p)

    lines(1:p$niter,b$traps[,1])

#or just run it a bunch of times since the model is stochastic
p$saturationThresholdStart = 5
time.to.max=c()
max.catch = c()
  realizations = 50
    plot(1:p$niter,xlab='Time',ylab='N Caught',ylim=c(0,15),type='n')

          for(i in 1:realizations){
                  a = SimulateLobsterMovement(p=p)
                   lines(1:p$niter,a$traps[,1])
                time.to.max = c(time.to.max , which.max(a$traps$V1))
                max.catch = c(max.catch , max(a$traps$V1))
                }

hist(time.to.max)
hist(max.catch)


p$saturationThresholdStart = 8
time.to.max8=c()
max.catch8 = c()
  realizations = 50
    
          for(i in 1:realizations){
                  a = SimulateLobsterMovement(p=p)
                   lines(1:p$niter,a$traps[,1],col='red')
                time.to.max8 = c(time.to.max8 , which.max(a$traps$V1))
                max.catch8 = c(max.catch8 , max(a$traps$V1))
                }

hist(time.to.max8)
hist(max.catch8)

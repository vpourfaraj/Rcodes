library(raster)
library(ggplot2)
library(reshape2)

# A function to get the distance of each lobster from the trap
# By default the trap is at (5,5)
# The function needs the x and y coordinate of the lobster
# adding trap location makes it generic
distanceToTrapCalculator<- function(xLobster, yLobster,xtrap=5, ytrap=5){
    distanceToTrap<- sqrt((xLobster - xtrap)^2 + (yLobster -ytrap)^2)
    return(distanceToTrap)
}

#when a lobster is outside the influence of bait, we need a A function
#to randomly move the lobster for a fixed length = dStep 
# The function needs the current x and y coordinate of lobster + dStep (here 5)

randomMove<- function(xLobster, yLobster, dStep){
  
    randomAngle<- runif(n=1, min = 0, max=360) #selects a random angle of direction 
    xNew<- dStep * sin(randomAngle * pi / 180) + xLobster #moves dStep based on the angle
    yNew<- dStep * cos(randomAngle * pi / 180) + yLobster
    
    # Check if the lobster is located outside of the grid
    # I think its ok for them to move outside the grid, so might not be necessary to have this step
    #if(xNew < 0 ){
    #  xNew <- 0
    #}
    #if(xNew > 10 ){
    #  xNew <- 10
    #}
    #if(yNew < 0 ){
    #  yNew <- 0
    #}
    #if(yNew > 10 ){
    #  yNew <- 10
    #}
    
    return( list(EASTING = xNew, NORTHING = yNew) )
}

# Anothe function is needed to  move the lobster for a fixed length = dStep toward the trap (+some randomness) 
#if it is within bait influence
# The function needs the current x and y coordinate of lobster + dStep + distance of the lobster from the bait 
#to determine if the condition for directionalMove is met
    #S = 0.993
    #atTime=0
    #ZoI = radius_of_influence =1
    #for each updated time step
    #ZoI = ZoI * S 

directionalMove<- function(xLobster, yLobster, dStep, distanceToTrap, xtrap=5, ytrap=5,radius_of_influence=15,ZoI=0){
  # Calculating tethaT, !!!!!!!! the concept of angle to trap (ask Adam!)
  #it looks like the quadrant of lobster matters here?
  #AMC - it shouldn't matter, I would just use the bearing function from geosphere package to get thetaT
  
  #if (xLobster < 5){  # for lobsters in 2nd and 3rd quarter
  #  degree<- (5 - xLobster)/ distanceToTrap
  #  tethaT<- acos(( (degree * pi) / 180 ) )
  #}
  #if (xLobster > 5){  # for lobsters in 1st and 4th quarter( not sure!!!)
  #  degree<- (xLobster - 5)/ distanceToTrap
  #  tethaT<- 180 - acos( (degree * pi) / 180 ) 
  #}
  thetaT = bearing(c(xLobster,yLobster),c(xtrap,ytrap))
  
  # Calculating tethaR
    #ZoI is adjusted each time step for decay (shrinkage) of radius of influence (equation 3)
    b = 1 + 0.9 * (distanceToTrap - ZoI) / radius_of_influence
    thetaR = -180:180
    P = 1/(180^b) * abs(thetaR) ^ b
    Prtheta_R = (1-P) / sum(1-P)
    theta_r = sample(thetaR,size = 1, prob=Prtheta_R)  
  
  
   tetha      <- thetaT + theta_r
  
  xNew   <- dStep * sin(tetha * pi / 180) + xLobster
  yNew   <- dStep * cos(tetha * pi / 180) + yLobster
  
  # Check if the lobster is going to be outside of the grid
  #AMC again I think its ok for lobsters to move outside grid
  #if(xNew < 0 ){
  #  xNew <- 0
  #}
  #if(xNew > 10 ){
  #  xNew <- 10
  #}
  #if(yNew < 0 ){
  #  yNew <- 0
  #}
  #if(yNew > 10 ){
  #  yNew <- 10
  #}
  
  
  return( list(EASTING = xNew, NORTHING = yNew) )
  
}
# A function to get the coordinate at time t, update the coordinate, and return the coordinate at time t + 1
updateGrid    = function(lobsterCoordinates,radius_of_influence=15, dstep = 5){
  
  # Takes the argument lobsterCoordinates as the x and y coordinates and updates it accordingly
  numberOfLobsters <- nrow(lobsterCoordinates)
  xNew <- vector(mode = 'numeric', length = numberOfLobsters)
  yNew <- vector(mode = 'numeric', length = numberOfLobsters)
  
  # Loops over each lobster and update their coordinate accordingly
  for( lobsterIndex in 1:numberOfLobsters ){
    
    xOld <- lobsterCoordinates[lobsterIndex,1]
    yOld <- lobsterCoordinates[lobsterIndex,2]
    distanceToTrap <- distanceToTrapCalculator(xLobster = xOld,yLobster = yOld)
    
    if( distanceToTrap > radius_of_influence){
      temp <- randomMove(xLobster = xOld , yLobster = yOld , dStep = dstep)
      xNew[lobsterIndex] <- temp$EASTING
      yNew[lobsterIndex] <- temp$NORTHING
    }else{
      temp <- directionalMove(xLobster = xOld , yLobster = yOld , distanceToTrap = distanceToTrap, dStep = dstep)
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

#density per grid 100 grids
# each element in the vector represents the starting value in each grid
initialGrid = rpoisD(n=100,lambda=.1, D = 3)
LobsterStart = data.frame(EASTING = rep(1:10,times=10), NORTHING = rep(1:10,each=10), Lobs = initialGrid)

# setting up Grid extents 
xmin=0 
ymin=0
xmax=10
ymax=10

grid_extent=extent(xmin, xmax, ymin, ymax)

# Create grid and get xy coordinates from which random points(lobsters) are extracted
simulated_raster <- raster(ext=grid_extent, resolution=1)
values(simulated_raster)= initialGrid # initial density of lobsters is given to raster  

# Extract cell centre coordinates
x_centres=xFromCol(simulated_raster)
y_centres=yFromRow(simulated_raster)

# Select some random points for traps
random_point_count=0
random_point_sample_number=7 # the number of points we want based on initial grid
random_points = NULL
while (random_point_count < random_point_sample_number){
  
  easting_random=sample(x_centres, 1)
  northing_random=sample(y_centres, 1)
  
  if (exists("random_points")==TRUE){
    # random_points dataframe already exists
    random_points=rbind(random_points, 
                        c("EASTING"=easting_random,
                          "NORTHING"=northing_random))
  } else {
    # random_points dataframe doesn't exist yet...")
    random_points=data.frame("EASTING"=easting_random,"NORTHING"=northing_random)	
  }
  
  random_point_count=random_point_count+1
}

#AMC need to include the catching process

#initialxyCoordinate <- subset(LobsterStart,Lobs>0)
initialxyCoordinate  = random_points
coordinatesOverTime <- list()
coordinatesOverTime[[1]] <- initialxyCoordinate #AMC not sure why this and not the output of simulated_raster
for(t in 2:100){
  coordinatesOverTime[[t]] <- updateGrid( lobsterCoordinates = coordinatesOverTime[[t-1]] )
  
  par( mfrow=c(1,2) ) # create a plot with 1 row and 2 columns to show plots side by side
  plot( coordinatesOverTime[[t-1]], pch = rownames(coordinatesOverTime[[t-1]]), xlim = c(0,10), ylim = c(0,10), main = paste0('Time = ', t-1) )
  points(x = 5, y = 5, col = 'red')
  plot( coordinatesOverTime[[t]],   pch = rownames(coordinatesOverTime[[t]]),   xlim = c(0,10), ylim = c(0,10), main = paste0('Time = ', t) )
  points(x = 5, y = 5, col = 'red')
  
}




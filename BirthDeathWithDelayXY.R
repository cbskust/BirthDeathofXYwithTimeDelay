## Stochastic Simulation
TimeDelayGillespie <- function(lambda1, mu1, lambda2, mu2, alpha1, beta1, alpha2, beta2, repnum = 1000){
  X <- 0
  XList <- rep(NA, repnum)
  Y <- 0
  YList <- rep(NA, repnum)
  currentTime <- 0
  TList <- rep(NA, repnum)
  stackTimeX <- c()
  stackTimeY <- c()
  
  for (i in 1:repnum){
    a1 <- lambda1
    a2 <- mu1 * X
    a3 <- lambda2 * X # a3 = lambda2 (X^n / (theta^n + X^n))
    a4 <- mu2 * Y
    a0 <- sum(a1,a2,a3,a4)
    # r2 <- runif(1)
    currentTime <- currentTime + rexp(1, rate = a0)
    
    stackTimeX <- sort(stackTimeX)  
    stackTimeY <- sort(stackTimeY)
    if(!(is.null(stackTimeX) & is.null(stackTimeY))){
       minStack <- min(stackTimeX, stackTimeY)
    } else {
       minStack <- Inf
    }
    if (currentTime < minStack){                                    
      r1 <- runif(1)
      if (r1 < a1/a0){
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- currentTime
        stackTimeX <- c(stackTimeX, currentTime + k*rgamma(n=1, shape = alpha1, scale = beta1))
      } else if (r1 < (a1+a2)/a0){
        X <- X-1;
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- currentTime
      } else if (r1 < (a1+a2+a3)/a0){
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- currentTime
        stackTimeY <- c(stackTimeY, currentTime + k*rgamma(n=1, shape = alpha2, scale = beta2))
      } else {
        Y <- Y-1;
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- currentTime
      }
    } else{
      if (min(stackTimeX) < min(stackTimeY)){
        X <- X+1;
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- minStack
        currentTime <- minStack
        stackTimeX <- stackTimeX[-1]
      } else {
        Y <- Y+1;
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- minStack
        currentTime <- minStack
        stackTimeY <- stackTimeY[-1]
      }
    }
  }
  my_list <- list("XList" = XList, "YList" = YList, "TList" = TList)
  return(my_list)
}

## Deterministic Simulation

TimeDelayDet <- function(t, y, parms){
  if (t<6){
    lagOne <- 0
    lagX <- 0
  }
  else {
    lagOne <- parms[1]
    lagX <- lagvalue(t-6, 1) #lagvalue(t-6, 1)
  }
  dy1 <- lagOne - parms[2] * y[1];
  dy2 <- parms[3] * lagX - parms[4] * y[2];
  list(c(dy1,dy2))
}


## Simulation Implementation

lambda1 <- 1 #Birth rate of X
mu1 <- 0.05 #Death rate of X
lambda2 <- 1.5#Birth rate of Y
mu2<- 0.05 #Birth rate of Y

alpha1 <- 18/5 #shape parameter of the gamma distribution for a time delay of X
beta1 <- 5/3 #scale parameter of the gamma distribution for a time delay of X
# mean 6, var 10
alpha2 <- 18/5 #shape parameter of the gamma distribution for a time delay of X
beta2 <- 5/3 #scale parameter of the gamma distribution for a time delay of X
# mean 6, var 10
k <- 1
theta <- 10
n <- 2

AA = TimeDelayGillespie(lambda1 = lambda1, mu1 = mu1, lambda2 = lambda2, mu2 = mu2, alpha1 = alpha1, beta1 = beta1, alpha2 = alpha2, beta2 = beta2, repnum = 100000)

parameters <- c(lambda1 = lambda1, mu1 = mu1, lambda2 = lambda2, mu2 = mu2)
state <- c(0, 0) 
times <- seq(0,1200, by = 0.1)
#install.packages("deSolve")
library(deSolve)
out <- dede(y = state,  times = times, func = TimeDelayDet, parms = parameters)

par(mfrow=c(2,1))
plotX <- plot(AA$TList, AA$XList, cex=0.1)
# lines(AA$TList, AA$XList)
lines(out[, 1], out[, 2], col = 'red')
plotY <- plot(AA$TList, AA$YList, cex=0.1)
# lines(AA$TList, AA$YList)
lines(out[, 1], out[, 3], col = 'red')

## Convert the time series X to time series which have equal time step size (default 0.1)
EqualTimeStep <- function(tspan, x, tstepsize = 0.1){
  # tspan: a vector of time points, x: a vector of data points corresponding 'tspan', tstepsize: time step size of output series.
  Tmax <- max(tspan)
  eqXList <- rep(0, floor(1/tstepsize*Tmax))
  currentT = 0;
  prevT = 0;
  for(tt in 1:floor(1/tstepsize*Tmax)){
    for (jj in 2:length(tspan)){
      currentT <- floor(tspan[jj]* 1/tstepsize)
      if (prevT != currentT){
        eqXList[(prevT+1):currentT] <- x[jj-1]
        prevT <- currentT
      }
    }
    if (eqXList[length(eqXList)] == 0){
      eqXList[length(eqXList)] <- x[length(x)]
    }
  }
  eqTList <- seq(0, Tmax, by = tstepsize)
  my_list <- list("eqXList" = eqXList, "eqTList" = eqTList)
  return(my_list)
}


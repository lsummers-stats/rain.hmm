exphmm <- function(numDays, numYears=1, numLoc=1, numStates, numMix, initDist, tmat, zeta,lambda){
  v.sim <- matrix(nrow = numDays,ncol = numYears)
  for(year in 1:numYears)
  {
    u <- runif(numDays)
    v.sim[,year] <- generate.mc(u,initDist,tmat)
  }
  v.sim.vector <- as.vector(v.sim)
  y.sim=array(dim = c(numDays,numYears,numLoc))
  for(year in 1:numYears)
    for(l in 1:numLoc){
      mix.unif <- matrix(runif(numDays),nrow = numDays)
      for(t in 1:numDays){
        which.mixture <- rej.state(u=mix.unif[t],alpha = zeta[v.sim[t,year],,l])
        y.sim[t,year,l] <- ifelse(which.mixture==1,0,rexp(1,lambda[v.sim[t,year],which.mixture-1,l]))
      }
    }
  y = list(y.sim,v.sim.vector)
  return(y)
}

rejstate <- function(u, pi){
  m <- length(pi)
  I <- 0
  for(i in 1:(m-1)){
    I <- I + pi[i]
    if(u < I){
      return(i)
    }
  }
  return(m)
}

gammahmm <- function(numDays, numYears=1, numLoc=1, numStates, numMix, initDist, tmat, zeta, lambda, omega){
  v.sim <- matrix(nrow = numDays,ncol = numYears)
  for(year in 1:numYears)
  {
    u <- runif(numDays)
    v.sim[,year] <- generate.mc(u,initDist,tmat)
  }
  v.sim.vector <- as.vector(v.sim)
  y.sim=array(dim = c(numDays,numYears,numLoc))
  for(year in 1:numYears)
    for(l in 1:numLoc){
      mix.unif <- matrix(runif(numDays),nrow = numDays)
      for(t in 1:numDays){
        which.mixture <- rej.state(u=mix.unif[t],alpha = zeta[v.sim[t,year],,l])
        y.sim[t,year,l] <- ifelse(which.mixture==1,0,rgamma(1,omega[v.sim[t,year],which.mixture-1,l],
                                                            lambda[v.sim[t,year],which.mixture-1,l]))
      }
    }
  y = list(y.sim,v.sim.vector)
  return(y)
}

#' HMM Simulator
#'
#' @description Creates synthetic data based on which distribution you choose.
#'
#' @export

generate.hmm <- function(D, Y = 1, L = 1, S, M = 2, initDist, tmat, zeta,dist = c("exp", "gamma"), lambda, omega){
  if(dist == "exp"){
    exphmm(numDays = D, numYears= Y, numLoc= L, numStates = S, numMix = M, initDist, tmat, zeta,lambda)
  }
  if(dist == "gamma"){
    gammahmm(numDays = D, numYears= Y, numLoc= L, numStates = S, numMix = M, initDist, tmat, zeta,lambda, omega)
  }
  print("please select a model type: exp for exponential or gamma for gamma")
}

#' Markov Chain Generator
#' @description
#' Used in generating a markov chain
#'
#' @export

generate.mc <- function(x, pi, transmat){
  Z <- c()
  n <- length(x)
  Z[1] <- rejstate(x[1], pi)
  for(i in 2:n){
    Z[i] <- rejstate(x[i], transmat[Z[i-1],])
  }
  return(Z)
}

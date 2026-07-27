#' Viterbi Decoding Function
#'
#' @description
#' The Viterbi Algorithm applied to Rain Models.
#'
#' @param b Matrix of Emission Probabilities
#'
#' @param a Matrix of Transition Probabilities
#'
#' @param initDist The Vector of Initial Probabilities for each state.
#'
#' @param numDays Number of Days in the Data
#'
#' @param numDays Number of Years recorded
#'
#' @param numDays Number of States decided
#'
#' @returns A vector that represents the best state sequence.
#'
#' @export


viterbi.decoding <- function(b, a, initDist, numDays, numYears, numStates){
  delta <- array(dim=c(numDays,numStates,numYears))
  psi <- array(dim=c(numDays,numStates,numYears))
  q_star <- matrix(nrow = numDays,ncol = numYears)
  delta[1,,] <- apply(as.matrix(b[1,,]),2,function(x)t(x)+log(initDist))
  psi[1,,] <- 0

  for(t in 2:numDays){
    for(j in 1:numStates){
      transition <- as.matrix(delta[t-1,,]) + log(a[,j])
      delta[t,j,] <- apply(transition, 2, max) + b[t,j,]
      for(n in 1:numYears)
        psi[t,j,n] <- min(which(transition[,n]==apply(transition, 2, max)[n]))
    }
  }
  P_star <- apply(as.matrix(delta[numDays,,]),2,max)
  tmp <- t(delta[numDays,,]) == P_star
  q_star[numDays,] <- apply(tmp, 1, function(x)min(which(x==1)))
  for(t in numDays:2)
    q_star[t-1,] <- psi[t,q_star[t],]
  return(as.vector(q_star))
}

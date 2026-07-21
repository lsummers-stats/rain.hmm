#' Forward Recursion Algorithm
#'
#'  @description
#'
#'  A function used to calculate forward variables.
#'
#'  @details
#'
#'  This function is left to show the programming done for the Forward step in the Forward-Backward Algorithm
#'  made for a Hidden Markov Model.
#'
#'  @param numObs Number of Observations
#'
#'  The total number of observations in your data set.
#'
#'  @param numStates Number of States
#'
#'  The total number of states planned for the Hidden Markov Model
#'
#'  @param initDist Initial Distribution
#'
#'  The initial probability kernel. If no value is given, it will automatically give
#'  the first state a probability of 1, and all states a probability of 0.
#'
#'  @param a State Transition Matrix
#'
#'  The current state transition matrix. Input can be either a matrix or an array.
#'
#'  @param b State Emission Matrix
#'
#'  The current state emission matrix. Input can be either a matrix or an array.
#'
#'  @returns A list holding two objects:
#'  * `fvar_tilde` : A matrix containing all your Forward Variables. Each row corresponds to each observation
#'  and each column corresponds to each state.
#'  * `ct` : A vector containing the constant used to normalize `F_tilde` at each point in time `t`.
#'
#'  @examples
#'  n = 100
#'  K = 3
#'  A = matrix(0,nrow = S, ncol = S)
#'  B = array(1, dim = c(n,K))
#'  fvars <- forward_recursion(numObs = n, numStates = K, a = A, b = B)
#'  print(fvars$fvar_tilde)
#'  print(fvars$ct)
#'
#'  @export


forward_recursion <- function(numObs, numStates, initDist = c(1,rep(0,K-1)), a = NULL, b = NULL){
  ct <- rep(1,numObs)
  fvar <- matrix(nrow = numObs,ncol = numStates)
  fvar_star <- matrix(nrow = numObs,ncol = numStates)
  fvar_tilde <- matrix(nrow = numObs,ncol = numStates)
  fvar[1,] <- initDist*b[1,]
  ct[1] <- 1/sum(fvar[1,])
  fvar_tilde[1,] <- fvar[1,]*ct[1]

  for(t in 2:numObs){
    for(j in 1:numStates){
      fvar[t,j] <- sum(fvar[t-1,]*a[,j])*b[t,j]
      fvar_star[t,j] <- sum(fvar_tilde[t-1,]*a[,j])*b[t,j]
      #print(t)
    }
    ct[t] <- 1/sum(fvar_star[t,])
    fvar_tilde[t,] <- fvar_star[t,]*ct[t]
  }
  Forward <- list(fvar_tilde, ct)
  return(Forward)
}

#' Backward Recursion Algorithm
#'
#'  @description
#' A function used to calculate backward variables.
#'
#'  @details
#' This function is left to show the programming done for the Backward step in the Forward-Backward Algorithm
#' made for a Hidden Markov Model.
#'
#'  @param numObs Number of Observations
#'
#'  The total number of observations in your data set.
#'
#'  @param numStates Number of States
#'
#'  The total number of states planned for the Hidden Markov Model
#'
#'  @param a State Transition Matrix
#'
#'  The current state transition matrix. Input can be either a matrix or an array.
#'
#'  @param b State Emission Matrix
#'
#'  The current state emission matrix. Input can be either a matrix or an array.
#'
#'  @param ct Normalizing Constant
#'
#'  A vector containing constants normalizing `F_tilde`, the Forward Variables. This vector is typically
#'  found from running the Forward part of the Forward Algorithm.
#'
#'  @returns A matrix of all your backward variables:
#'  * The rows correspond to the observations
#'  * The columns correspond to the states.
#'
#'  @examples
#'  n = 100
#'  K = 3
#'  A = matrix(0,nrow = S, ncol = S)
#'  B = array(1, dim = c(n,K))
#'  fvars <- forward_recursion(numObs = n, numStates = K, a = A, b = B, ct = fvars[[2]])
#'  backward_recursion(numObs = n, numStates = K, )
#'
#'  @export

backward_recursion <- function(numObs, numStates, a = NULL, b = NULL, ct = NULL){
  bvar_star <- matrix(nrow = numObs,ncol = numStates)
  bvar_tilde <- matrix(nrow = numObs,ncol = numStates)
  bvar_tilde[numObs,] <- ct[numObs]
  for(t in (numObs-1):1){
    for(j in 1:numStates){
      bvar_star[t,j] <- sum(a[j,]*bvar_tilde[t+1,]*b[t+1,])
    }
    bvar_tilde[t,] <- bvar_star[t,]*ct[t]
  }
  return(bvar_tilde)
}


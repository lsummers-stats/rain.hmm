#' ELBO Function
#'
#' @description
#' A function designed to calculate an approximation of the ELBO.
#'
#' @details
#' Since the ELBO is computationally difficult, or sometimes impossible, the function was designed to approximated it.
#' This function is used within the Variational Bayes EM loop to check for improvements.
#'
#' @param numStates Number of States laid out in the model.
#'
#' @param numMix Number of Mixtures laid out in the model.
#'
#' @param numLoc Number of locations collected in the data.
#'
#' @param stateProb Matrix of stationary probabilities.
#'
#' @param mixProb Matrix of values used to weight the mixtures.
#'
#' @param initProb Matrix of Initial Values for each state.
#'
#' @param jtTransMat Transition Probability Matrix.
#'
#' @param ct A matrix of values used to normalize the forward variables.
#'
#' @param xi The Current Matrix for the Initial States.
#'
#' @param alpha The Current Matrix for the Transitions.
#'
#' @param zeta The Current Matrix for the Mixtures.
#'
#' @param gamma_shape The Current matrix for the Shape Parameters
#'
#' @param gamma_rate The Current Matrix for the Rate Parameters
#'
#' @param obs The vector of observations.
#'
#' @param h A matrix of values used in calculating the ELBO.
#'
#' @return A single value representing the ELBO for the given model parameters.
#'
#' @export


ELBO = function(numStates, numMix, numLoc, stateProb, mixProb, initProb, jtTransMat, ct, xi, alpha, zeta, gamma_shape, gamma_rate, obs, h){
  kl_C <- 0
  kl_A <- 0
  kl_theta <- 0
  kl_pi <- 0
  for(j in 1:numStates){
    for(l in numLoc){
      kl_C <- kl_C + sum(stateProb[,j,]*mixProb[,j,1,,l])*(digamma(zeta[j,1,l]) - digamma(sum(zeta[j,,l]))) +
        lgamma(sum(zeta[j,,l])) - sum(lgamma(zeta[j,,l])) #-
      #lgamma(sum(zeta_0[j,])) + sum(lgamma(zeta_0[j,]))
      for(m in 2:numMix){
        kl_C <- kl_C + sum(stateProb[,j,]*mixProb[,j,m,,l])*(digamma(zeta[j,m,l]) - digamma(sum(zeta[j,,l])) )
        kl_theta <- kl_theta + sum(stateProb[,j,]*mixProb[,j,m,,l])*(digamma(gamma_shape[j,m-1,l]) - log(gamma_rate[j,m-1,l])) -
          sum(stateProb[,j,]*mixProb[,j,m,,l]*obs[,,l])*(gamma_shape[j,m-1,l]/gamma_rate[j,m-1,l] ) + h[j,m-1,l] #- h_j0[j,p-1]
      }
    }
    kl_A <- kl_A + sum(rowSums(jtTransMat,dims = 2)[j,]*( digamma(alpha[j,]) - digamma(sum(alpha[j,])) )) +
      lgamma(sum(alpha[j,])) - sum(lgamma(alpha[j,]))   #-
    #lgamma(sum(alpha_0[j,])) + sum(lgamma(alpha_0[j,]))
  }
  kl_pi <-  sum(apply(initProb, 2, sum)*(digamma(xi) - digamma(sum(xi)))) +
    lgamma(sum(xi)) - sum(lgamma(xi)) #- lgamma(sum(pi_0)) + sum(lgamma(pi_0))
  elboresult <- - sum(log(ct)) - kl_theta - kl_A  - kl_C  - kl_pi
}

#' DIC Function
#'
#' @description
#' A function designed to calculate an approximation of the DIC.
#'
#' @details
#' Since the DIC is computationally difficult, or sometimes impossible, the function was designed to approximated it.
#' This function is used within the Variational Bayes EM loop to check for improvements.
#'
#' @param numStates Number of States laid out in the model.
#'
#' @param numMix Number of Mixtures laid out in the model.
#'
#' @param numLoc Number of locations collected in the data.
#'
#' @param stateProb Matrix of stationary probabilities.
#'
#' @param mixProb Matrix of values used to weight the mixtures.
#'
#' @param initProb Matrix of Initial Values for each state.
#'
#' @param jtTransMat Transition Probability Matrix.
#'
#' @param ct A matrix of values used to normalize the forward variables.
#'
#' @param xi The Current Matrix for the Initial States.
#'
#' @param alpha The Current Matrix for the Transitions.
#'
#' @param zeta The Current Matrix for the Mixtures.
#'
#' @param gamma_shape The Current matrix for the Shape Parameters
#'
#' @param gamma_rate The Current Matrix for the Rate Parameters
#'
#' @return A list of two objects:
#' *`DIC`: A single value representing the DIC for the given model parameters.
#' *`pd`: A value used in calculating the DIC.
#'
#' @export

DIC <- function(numStates, numMix, numLoc, stateProb, mixProb, initProb, jtTransMat, ct, xi, alpha, zeta, gamma_shape, gamma_rate){
  pd0 <- 0
  pd1 <- 0
  pd2 <- 0
  for(j in 1:numStates){
    for(l in 1:numLoc){
      pd0 <- pd0 + sum(stateProb[,j,]*mixProb[,j,1,,l])*(log(zeta[j,1,l]) - log(sum(zeta[j,,l])) -
                                                           digamma(zeta[j,1,l]) + digamma(sum(zeta[j,,l])) )

      for(m in 2:numMix){
        pd2 <- pd2 + sum(stateProb[,j,]*mixProb[,j,m,,l])*( log(gamma_shape[j,m-1,l]) - digamma(gamma_shape[j,m-1,l]) +
                                                              log(zeta[j,m,l]) - log(sum(zeta[j,,l])) - digamma(zeta[j,m,l]) + digamma(sum(zeta[j,,l])) )
      }
    }
    for(k in 1:numStates)
      pd1 <- pd1 + rowSums(jtTransMat,dims = 2)[j,k]*( log(alpha[j,k]) - log(sum(alpha[j,])) - digamma(alpha[j,k]) + digamma(sum(alpha[j,])) )

  }
  pd3 <- sum(initProb*(log(xi) - log(sum(xi)) - digamma(xi) + digamma(sum(xi))))
  pd <- pd0+sum(pd1)+pd2+pd3
  dic <- 4*(pd)  + 2*sum(log(ct))
  result = list('pd' = pd, 'dic' = dic)
}

#' Stochastic ELBO Function
#'
#' @description
#' A function designed to calculate an approximation of the ELBO.
#'
#' @details
#' Since the ELBO is computationally difficult, or sometimes impossible, the function was designed to approximated it.
#' This function is used within the Variational Bayes EM loop to check for improvements, with SVI in mind.
#'
#' @param numStates Number of States laid out in the model.
#'
#' @param numMix Number of Mixtures laid out in the model.
#'
#' @param numLoc Number of locations collected in the data.
#'
#' @param stateProb Matrix of stationary probabilities.
#'
#' @param mixProb Matrix of values used to weight the mixtures.
#'
#' @param initProb Matrix of Initial Values for each state.
#'
#' @param jtTransMat Transition Probability Matrix.
#'
#' @param ct A matrix of values used to normalize the forward variables.
#'
#' @param xi The Current Matrix for the Initial States.
#'
#' @param alpha The Current Matrix for the Transitions.
#'
#' @param zeta The Current Matrix for the Mixtures.
#'
#' @param gamma_shape The Current matrix for the Shape Parameters
#'
#' @param gamma_rate The Current Matrix for the Rate Parameters
#'
#' @param obs The vector of observations.
#'
#' @param h A matrix of values used in calculating the ELBO.
#'
#' @return A single value representing the ELBO for the given model parameters.
#'
#' @export

StoELBO = function(numStates, numMix, numLoc, stateProb, mixProb, initProb, jtTransMat, ct, xi, alpha, zeta, gamma_shape, gamma_rate, obs, h){
  kl_C <- 0
  kl_A <- 0
  kl_theta <- 0
  kl_pi <- 0
  for(j in 1:numStates){
    for(l in numLoc){
      kl_C <- kl_C + sum(stateProb[,j]*mixProb[,j,1,l])*(digamma(zeta[j,1,l]) - digamma(sum(zeta[j,,l]))) +
        lgamma(sum(zeta[j,,l])) - sum(lgamma(zeta[j,,l]))
      for(m in 2:numMix){
        kl_C <- kl_C + sum(stateProb[,j]*mixProb[,j,m,l])*(digamma(zeta[j,m,l]) - digamma(sum(zeta[j,,l])) )
        kl_theta <- kl_theta + sum(stateProb[,j]*mixProb[,j,m,l])*(digamma(gamma_shape[j,m-1,l]) - log(gamma_rate[j,m-1,l])) -
          sum(stateProb[,j]*mixProb[,j,m,l]*obs[,l])*(gamma_shape[j,m-1,l]/gamma_rate[j,m-1,l] ) + h[j,m-1,l]
      }
    }
    kl_A <- kl_A + sum(rowSums(jtTransMat,dims = 2)[j,]*(digamma(alpha[j,]) - digamma(sum(alpha[j,])) )) +
      lgamma(sum(alpha[j,])) - sum(lgamma(alpha[j,]))
  }
  kl_pi <-  sum(apply(initProb, 2, sum)*(digamma(xi) - digamma(sum(xi)))) +
    lgamma(sum(xi)) - sum(lgamma(xi))
  elboresult <- - sum(log(ct)) - kl_theta - kl_A  - kl_C  - kl_pi
}


#' DIC Function
#'
#' @description
#' A function designed to calculate an approximation of the DIC.
#'
#' @details
#' Since the DIC is computationally difficult, or sometimes impossible, the function was designed to approximated it.
#' This function is used within the Variational Bayes EM loop to check for improvements, with SVI in mind.
#'
#' @param numStates Number of States laid out in the model.
#'
#' @param numMix Number of Mixtures laid out in the model.
#'
#' @param numLoc Number of locations collected in the data.
#'
#' @param stateProb Matrix of stationary probabilities.
#'
#' @param mixProb Matrix of values used to weight the mixtures.
#'
#' @param initProb Matrix of Initial Values for each state.
#'
#' @param jtTransMat Transition Probability Matrix.
#'
#' @param ct A matrix of values used to normalize the forward variables.
#'
#' @param xi The Current Matrix for the Initial States.
#'
#' @param alpha The Current Matrix for the Transitions.
#'
#' @param zeta The Current Matrix for the Mixtures.
#'
#' @param gamma_shape The Current matrix for the Shape Parameters
#'
#' @param gamma_rate The Current Matrix for the Rate Parameters
#'
#' @return A list of two objects:
#' *`DIC`: A single value representing the DIC for the given model parameters.
#' *`pd`: A value used in calculating the DIC.
#'
#' @export
StoDIC <- function(numStates, numMix, numLoc, stateProb, mixProb, initProb, jtTransMat, ct, xi, alpha, zeta, gamma_shape, gamma_rate){
  pd0 <- 0
  pd1 <- 0
  pd2 <- 0
  for(j in 1:numStates){
    for(l in 1:numLoc){
      pd0 <- pd0 + sum(stateProb[,j]*mixProb[,j,1,l])*(log(zeta[j,1,l]) - log(sum(zeta[j,,l])) -
                                                         digamma(zeta[j,1,l]) + digamma(sum(zeta[j,,l])) )

      for(m in 2:numMix){
        pd2 <- pd2 + sum(stateProb[,j]*mixProb[,j,m,l])*( log(gamma_shape[j,m-1,l]) - digamma(gamma_shape[j,m-1,l])+     log(zeta[j,m,l]) - log(sum(zeta[j,,l])) - digamma(zeta[j,m,l]) + digamma(sum(zeta[j,,l])) )
      }
    }
    for(k in 1:numStates)
      pd1 <- pd1 + rowSums(jtTransMat,dims = 2)[j,k]*( log(alpha[j,k]) - log(sum(alpha[j,])) - digamma(alpha[j,k]) + digamma(sum(alpha[j,])) )

  }
  pd3 <- sum(initProb*(log(xi) - log(sum(xi)) - digamma(xi) + digamma(sum(xi))))
  pd <- pd0+sum(pd1)+pd2+pd3
  dic <- 4*(pd)  + 2*sum(log(ct))
  result = list('pd' = pd, 'dic' = dic)
}

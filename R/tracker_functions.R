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

ELBO.gam = function(numStates, numMix, numLoc, stateProb, mixProb, initProb, jtTransMat, ct, xi, alpha, zeta, gamma_hyper, delta_hyper, theta_hyper, logbetaprior, logbetapost, exp_omega, exp_psi_omega, exp_lomega, obs, h){
  kl_C <- 0
  kl_A <- 0
  kl_theta <- 0
  kl_pi <- 0
  for(j in 1:numStates){
    for(l in numLoc){
      kl_C <- kl_C + sum(stateProb[,j,]*mixProb[,j,1,,l])*(digamma(zeta[j,1,l]) - digamma(sum(zeta[j,,l]))) +
        lgamma(sum(zeta[j,,l])) - sum(lgamma(zeta[j,,l]))
      for(m in 2:numMix){
        kl_C <- kl_C + sum(stateProb[,j,]*mixProb[,j,m,,l])*(digamma(zeta[j,m,l]) - digamma(sum(zeta[j,,l])) )
        kl_theta <- kl_theta +
          #gamma hyper diff
          sum(stateProb[,j,]*mixProb[,j,m,,l])*(exp_psi_omega[j,m-1,l] + exp_omega[j,m-1,l]*log(theta_hyper[j,m-1,l])) -
          #delta hyper diff
          sum(stateProb[,j,]*mixProb[,j,m,,l])*(exp_lomega[j,m-1,l]) -
          #theta hyper diff
          sum(stateProb[,j,]*mixProb[,j,m,,l]*obs[,,l])*gamma_hyper[j,m-1,l]/theta_hyper[j,m-1,l] +
          #beta hyper diff
          sum(logbetapost[j,m-1,l] - logbetaprior[j,m-1,l])*exp_omega[j,m-1,l]
      }
    }
    kl_A <- kl_A + sum(rowSums(jtTransMat,dims = 2)[j,]*( digamma(alpha[j,]) - digamma(sum(alpha[j,])) )) +
      lgamma(sum(alpha[j,])) - sum(lgamma(alpha[j,]))
  }
  kl_pi <-  sum(apply(initProb, 2, sum)*(digamma(xi) - digamma(sum(xi)))) +
    lgamma(sum(xi)) - sum(lgamma(xi))
  elboresult <- - sum(log(ct)) - kl_theta - kl_A  - kl_C  - kl_pi}

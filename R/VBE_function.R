VBE.exp <- function(numDays, numStates, numYears, numLoc, numMix, xi, alpha, zeta, gamma_shape, gamma_rate, var, obs){
  del_y0      <- var$del
  a_jk        <- var$a_jk # posterior state kernel
  b_tj        <- var$b_tj
  b_tjl       <- var$b_tjl
  a_1j        <- var$a_1j  # posterior initial probability kernel
  ct          <- var$ct
  fvar        <- var$fvar
  bvar        <- var$bvar
  b_star      <- var$b_star
  q_tj        <- var$q_tj # posterior probability of stationary distribution
  q_tjml      <- var$q_tjml
  q_jk        <- var$q_jk # posterior joint transition probability matrix
  q_1j        <- var$q_1j # posterior initial probability

  #### Update latent variables
  a_1j <- exp(digamma(xi) - digamma(sum(xi)))
  for(j in 1:numStates){
    a_jk[j,] <- exp(digamma(alpha[j,]) - digamma(sum(alpha[j,])))
    for(l in 1:numLoc)
      for(n in 1:numYears){
        tmp <- 0
        for(m in 2:numMix){
          tmp <- tmp + (1- del_y0[,n,l])*exp(digamma(zeta[j,m,l]) - digamma(sum(zeta[j,,l])) + digamma(gamma_shape[j,m-1,l]) - log(gamma_rate[j,m-1,l]) - obs[,n,l]*(gamma_shape[j,m-1,l]/gamma_rate[j,m-1,l] ))}
        b_tjl[,j,n,l] <- del_y0[,n,l]*exp(digamma(zeta[j,1,l]) - digamma(sum(zeta[j,,l])) ) + tmp
      }
  }
  for(l in 1:numLoc){
    b_tj[,,] <- b_tj[,,]*b_tjl[,,,l]
  }
  #### Forward and Backward recursions
  for(n in 1:numYears){
    fvar_out  <- forward_recursion(numObs = numDays, numStates = numStates, initDist = a_1j, a = a_jk, b = b_tj[,,n])
    fvar[,,n] <- fvar_out[[1]]
    ct[,n]    <- fvar_out[[2]]

    bvar[,,n] <- backward_recursion(numObs = numDays, numStates = numStates, a = a_jk, b = b_tj[,,n], ct=ct[,n])
  }

  #### Update the posterior state probabilities
  for(n in 1:numYears){
    q_tj[,,n] <- fvar[,,n]*bvar[,,n]/apply(bvar[,,n]*fvar[,,n],1,sum)
    for(j in 1:numStates)
      for(k in 1:numStates){
        q_jk[j,k,,n] <- fvar[-numDays,j,n]*a_jk[j,k]*b_tj[-1,k,n]*bvar[-1,k,n]
      }
    #q_jk[,,,m] <- q_jk[,,,m]/rep(colSums(q_jk[,,,m],dims = 2),each=K^2)
    for(t in 1:(numDays-1)){
      denom <- 0
      for(j in 1:numStates)
        for(k in 1:numStates){
          denom <- denom+q_jk[j,k,t,n]
        }
      q_jk[,,t,n] <- q_jk[,,t,n]/denom
    }
  }
  q_1j <- t(q_tj[1,,])

  #### Update the mixing probabilties
  for(l in 1:numLoc)
    for(n in 1:numYears)
      for(j in 1:numStates){
        b_star[,j,1] <- del_y0[,n,l]
        for(m in 2:numMix)
          b_star[,j,m] <- (1-del_y0[,n,l])*exp(digamma(zeta[j,m,l]) - digamma(sum(zeta[j,,l])) +
                                                 digamma(gamma_shape[j,m-1,l]) - log(gamma_rate[j,m-1,l]) - obs[,n,l]*(gamma_shape[j,m-1,l]/gamma_rate[j,m-1,l]) )
        sum_b <- rowSums(b_star[,j,],dims = 1)
        q_tjml[,j,-1,n,l] <- b_star[,j,-1]/sum_b
        q_tjml[,j,1,n,l] <- b_star[,j,1]
      }
  output = list('a_jk' = a_jk, 'b_tj' = b_tj, 'b_tjl' = b_tjl, 'ct' = ct,'q_1j' = q_1j, 'q_tj' = q_tj, 'q_tjml' = q_tjml, 'q_jk' = q_jk)
}

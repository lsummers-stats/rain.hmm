VBEmpfr.exp = function(numDays, numStates, numLoc, numMix, xi, alpha, zeta, gamma_shape, gamma_rate, var, obs){
  del_y0      <- var$del
  a_jk        <- var$a_jk # posterior state kernel
  b_tj        <- var$b_tj
  b_tjl       <- var$b_tjl
  a_1j        <- var$a_1j  # posterior initial probability kernel
  logct       <- var$ct
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
      for(m in 2:numMix){
        tmp <- 0
        tmp <- tmp + (1 - del_y0[,l])*exp(digamma(zeta[j,m,l]) - digamma(sum(zeta[j,,l])) + digamma(gamma_shape[j,m-1,l]) - log(gamma_rate[j,m-1,l]) - obs[,l]*(gamma_shape[j,m-1,l]/gamma_rate[j,m-1,l]))
        b_tjl[,j,l] <- del_y0[,l]*exp(digamma(zeta[j,1,l]) - digamma(sum(zeta[j,,l])) ) + tmp
      }
  }
  for(l in 1:numLoc){
    b_tj[,] <- b_tj[,]*b_tjl[,,l]
  }
  #### Forward and Backward recursions
  fvar_out  <- forward_recursion_mpfr(numObs = numDays, numStates = numStates, initDist = a_1j, a = a_jk, b = b_tj[,], precision = 120)
  fvar[,] <- fvar_out[[1]]
  logct <- fvar_out[[2]]
  bvar[,] <- backward_recursion_mpfr(numObs = numDays, numStates = numStates, a = a_jk, b = b_tj[,], ct=logct, precision = 120)
  q_tj[,] <- as.double(fvar[,]*bvar[,]/apply(bvar[,]*fvar[,],1,sum))
  for(j in 1:numStates)
    for(k in 1:numStates){
      q_jk[j,k,] <- fvar[-numDays,j]*a_jk[j,k]*exp(mpfr(b_tj[-1,k],120))*bvar[-1,k]
    }
  for(t in 1:(numDays-1)){
    denom <- mpfr(0,120)
    for(j in 1:numStates)
      for(k in 1:numStates){
        denom <- denom+q_jk[j,k,t]
      }
    q_jk[,,t] <- as.double(q_jk[,,t]/denom)
  }

  #### Update the posterior state probabilities

  q_1j <- t(q_tj[1,])

  #### Update the mixing probabilties
  for(l in 1:numLoc){
    for(j in 1:numStates){
      b_star[,j,1] <- del_y0[,l]
      for(m in 2:numMix){
        b_star[,j,m] <- (1-del_y0[,l])*exp(digamma(zeta[j,m,l]) - digamma(sum(zeta[j,,l])) +
                                             digamma(gamma_shape[j,m-1,l]) - log(gamma_rate[j,m-1,l]) - obs[,l]*(gamma_shape[j,m-1,l]/gamma_rate[j,m-1,l]) )}
      sum_b <- rowSums(b_star[,j,],dims = 1)
      fix = which(sum_b==0)
      sum_b[fix] = exp(-700)
      for(i in fix){
        b_star[fix,j,2] = exp(-700)
      }
      q_tjml[,j,-1,l] <- exp(log(b_star[,j,-1]) - log(sum_b))
      q_tjml[,j,1,l] <- b_star[,j,1]}
    fix = which(logct < 0)
    logct[fix] = exp(-700)
  }
  output = list('a_jk' = a_jk, 'b_tj' = b_tj, 'ct' = logct, 'q_1j' = q_1j, 'q_tj' = q_tj, 'q_tjml' = q_tjml, 'q_jk' = q_jk)
}

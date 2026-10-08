VBMS.exp = function(numStates, numLoc, numMix, numYears, numDays, xi, alpha, zeta, gamma_shape, gamma_rate, q_1j, q_tj, q_tjml, q_jk, obs, iter, priors){
  #Q's will be randomized already, so no need to edit the sums, just multiply by N
  N = numYears
  step = 1/iter
  gamma_jml    <- gamma_shape[[iter]] # posterior shape of exponential rate
  delta_jml    <- gamma_rate[[iter]] # posterior rate of exponential rate
  zeta_jl      <- zeta[[iter]] # posterior Dirichlet parameters for mixing probabilities
  alpha_j      <- alpha[[iter]] # posterior Dirichlet parameters for transition matrix rows
  xi_j         <- xi[[iter]] # posterior Dirichlet parameters for initial distribution
  gamma_prior      <- gamma_shape[[iter - 1]]
  delta_prior      <- gamma_rate[[iter - 1]]
  xi_prior         <- xi[[iter - 1]]
  alpha_prior      <- alpha[[iter - 1]]
  zeta_prior      <- xi[[iter - 1]]
  #### Update hyperparameters
  for(j in 1:numStates){
    xi_j[j] <- (1 - step)*xi_j[j] + (step)*(xi_prior[j] + sum(q_1j[j]))
    for(l in 1:numLoc){
      zeta_jl[j,1,l] <- (1 - step)*zeta_jl[j,1,l] + (step) * as.numeric(zeta_prior[j,1,l] + N*sum(q_tj[,j]*q_tjml[,j,1,l]))
      for(m in 2:numMix){
        zeta_jl[j,m,l] <- (1 - step)*zeta_jl[j,m,l] + (step) * as.numeric(zeta_prior[j,m,l] + N*sum(q_tj[,j]*q_tjml[,j,m,l]))
        gamma_jml[j,m-1,l] <- (1 - step)*gamma_jml[j,m-1,l]+(step)*as.numeric(gamma_prior[j,m-1,l] + N*sum(q_tj[,j]*q_tjml[,j,m,l]))
        delta_jml[j,m-1,l] <- (1 - step)*delta_jml[j,m-1,l] + (step)*as.numeric(delta_prior[j,m-1,l] + sum(q_tj[,j]*q_tjml[,j,m,l]*obs[,l]))
      }
    }
    for(k in 1:numStates)
      alpha_j[j,k] <- (1 - step)*(alpha_prior[j,k])+(step)*(alpha_prior[j,k] + sum(as.numeric(q_jk[j,k,])))
  }

  h_jml <- gamma_jml*log(delta_jml) - lgamma(gamma_jml)
  output = list('gamma_jml' = gamma_jml, 'delta_jml' = delta_jml, 'xi_j' = xi_j, 'alpha_j' = alpha_j, 'zeta_jl' = zeta_jl, 'h_jml' = h_jml)
}

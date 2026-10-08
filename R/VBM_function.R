VBM.exp = function(numStates, numLoc, numMix, xi, alpha, zeta, gamma_shape, gamma_rate, q_1j, q_tj, q_tjml, q_jk, obs, iter){
  gamma_jml    <- gamma_shape # posterior shape of exponential rate
  delta_jml    <- gamma_rate # posterior rate of exponential rate
  zeta_jl      <- zeta # posterior Dirichlet parameters for mixing probabilities
  alpha_j      <- alpha # posterior Dirichlet parameters for transition matrix rows
  xi_j         <- xi # posterior Dirichlet paramters for initial distribution
  #### Update hyperparameters
  for(j in 1:numStates){
    xi_j[j] <- xi[j] + sum(q_1j[,j])
    for(l in 1:numLoc){
      zeta_jl[j,1,l] <- zeta[j,1,l] + sum(q_tj[,j,]*q_tjml[,j,1,,l])
      for(m in 2:numMix){
        zeta_jl[j,m,l] <- zeta[j,m,l] + sum(q_tj[,j,]*q_tjml[,j,m,,l])
        gamma_jml[j,m-1,l] <- gamma_shape[j,m-1,l] + sum(q_tj[,j,]*q_tjml[,j,m,,l])
        delta_jml[j,m-1,l] <- gamma_rate[j,m-1,l] + sum(q_tj[,j,]*q_tjml[,j,m,,l]*obs[,,l])
      }
    }
    for(k in 1:numStates)
      alpha_j[j,k] <- alpha[j,k] + sum(q_jk[j,k,,])
  }

  h_jml <- gamma_jml*log(delta_jml) - lgamma(gamma_jml)
  output = list('gamma_jml' = gamma_jml, 'delta_jml' = delta_jml, 'xi_j' = xi_j, 'alpha_j' = alpha_j, 'zeta_jl' = zeta_jl, 'h_jml' = h_jml)
}

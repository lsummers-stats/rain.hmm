VBEM.gam = function(D, S, Y, L, M = 2, xi, alpha, zeta, gammah, deltah, thetah, logbetah, obs, maxiter) {
  elbo      <- rep(0,maxiter)
  elbo_old  <- -50000
  elbo[1]   <- -25000
  tol       <- 10^(-9)
  iter      <- 1
  improvement_elbo  <- (elbo_old-elbo[1])/elbo_old
  lambda.post = array(0,dim = c(S,M-1,L))
  omega.post = array(0,dim = c(S,M-1,L))
  zeta.post = array(0,dim = c(S,M,L))
  tmat.post = matrix(0,S,S)
  pi.post = rep(0,S)
  #Empty Variables
  del_y0      <- ifelse(obs==0,1,0)
  a_jk        <- matrix(0,nrow = S, ncol = S) # posterior state kernel
  b_tj        <- array(1, dim = c(D,S,Y)) # posterior emission kernel
  b_tjl       <- array(0,dim = c(D,S,Y,L))
  a_1j        <- matrix(0,nrow = Y,ncol = S) # posterior initial probability kernel

  ct          <- matrix(0,nrow = D,ncol = Y)
  fvar        <- array(0,dim = c(D,S,Y))
  bvar        <- array(0,dim = c(D,S,Y))
  b_star      <- array(0,dim = c(D,S,M))

  q_tj        <- array(0,dim = c(D,S,Y)) # posterior probability of stationary distribution
  q_tjml      <- array(0,dim = c(D,S,M,Y,L))
  q_jk        <- array(0,dim=c(S,S,D-1,Y)) # posterior joint transition probability matrix
  q_1j        <- matrix(0,nrow = Y, ncol = S) # posterior initial probability

  omega_norm <- array(1, c(S,M-1,L))
  omega_exp <- array(0, c(S,M-1,L))
  omega_lgamma <- array(0, c(S,M-1,L))
  omega_psi <- array(0, c(S,M-1,L))

  for(j in 1:S){
    for(l in 1:L){
      for(m in 2:M){
        omega_norm[j,m-1,l] <- omega_constant(gammah[j,m-1,l], deltah[j,m-1,l], thetah[j,m-1,l], logbetah[j,m-1,l])$value
        omega_exp[j,m-1,l] <- exp_omega(gammah[j,m-1,l], deltah[j,m-1,l], thetah[j,m-1,l], logbetah[j,m-1,l], con = omega_norm[j,m-1,l])
        omega_lgamma[j,m-1,l] <- exp_l_omega(gammah[j,m-1,l], deltah[j,m-1,l], thetah[j,m-1,l], logbetah[j,m-1,l], con = omega_norm[j,m-1,l])
        omega_psi[j,m-1,l] <- exp_psi_omega(gammah[j,m-1,l], deltah[j,m-1,l], thetah[j,m-1,l], logbetah[j,m-1,l], con = omega_norm[j,m-1,l])
      }
    }
  }
  xi.list    <- list('prior' = xi)
  alpha.list   <- list('prior' = alpha)
  gamma.list   <- list('prior' = gammah)
  delta.list   <- list('prior' = deltah)
  theta.list   <- list('prior' = thetah)
  logbeta.list   <- list('prior' = logbetah)
  mix.list   <- list('prior' = zeta)

  emptyfillers <- list('del' = del_y0, 'a_jk' = a_jk, 'b_tj' = b_tj, 'b_tjl' = b_tjl, 'a_1j' = a_1j, 'ct' = ct, 'fvar' = fvar, 'bvar' = bvar, 'b_star' = b_star, 'q_tj' = q_tj, 'q_tjml' = q_tjml, 'q_jk' = q_jk, 'q_1j' = q_1j)
  #Model Running
  while((abs(improvement_elbo) > tol | improvement_elbo <0) & iter<maxiter){
    VBEout <- VBE.gam(numDays = D, numStates = S, numYears = Y, numLoc = L, numMix = M, xi = xi, alpha = alpha, zeta = zeta, gamma = gammah, delta = deltah, logbeta = logbetah, theta = thetah, exp_omega = omega_exp, exp_psi_omega = omega_psi, exp_lomega = omega_lgamma, var = emptyfillers, obs)
    VBMout <- VBM.gam(numStates = S, numLoc = L, numMix = M, xi = xi, alpha = alpha, zeta = zeta, gamma_hyper = gammah, delta_hyper = deltah, theta_hyper = thetah, log_beta_hyper = logbetah, VBEout$q_1j, VBEout$q_tj, VBEout$q_tjml, VBEout$q_jk, obs)
    ## Update iter
    iter = iter+1
    ##Quick update
    xi <- VBMout$xi_j
    alpha <- VBMout$alpha
    zeta <- VBMout$zeta_jl
    gammah <- VBMout$gamma_jml
    deltah <- VBMout$delta_jml
    thetah <- VBMout$theta_jml
    logbetah <- VBMout$log_beta_jml
    for(j in 1:S){
      for(l in 1:L){
        for(m in 2:M){
          omega_norm[j,m-1,l] <- omega_constant(gammah[j,m-1,l], deltah[j,m-1,l], thetah[j,m-1,l], logbetah[j,m-1,l])$value
          omega_exp[j,m-1,l] <- exp_omega(gammah[j,m-1,l], deltah[j,m-1,l], thetah[j,m-1,l], logbetah[j,m-1,l], con = omega_norm[j,m-1,l])
          omega_lgamma[j,m-1,l] <- exp_l_omega(gammah[j,m-1,l], deltah[j,m-1,l], thetah[j,m-1,l], logbetah[j,m-1,l], con = omega_norm[j,m-1,l])
          omega_psi[j,m-1,l] <- exp_psi_omega(gammah[j,m-1,l], deltah[j,m-1,l], thetah[j,m-1,l], logbetah[j,m-1,l], con = omega_norm[j,m-1,l])
        }
      }
    }
    emptyfillers$a_jk <- VBEout$a_jk
    emptyfillers$b_tj <- VBEout$b_tj
    emptyfillers$b_tjl <- VBEout$b_tjl
    emptyfillers$q_tj <- VBEout$q_tj
    emptyfillers$q_tjml <- VBEout$q_tjml
    emptyfillers$q_jk <- VBEout$q_jk
    emptyfillers$q_1j <- VBEout$q_1j
    xi.list[[iter]]    <- xi
    alpha.list[[iter]]    <- alpha
    gamma.list[[iter]]    <- gammah
    delta.list[[iter]]    <- deltah
    theta.list[[iter]]    <- thetah
    logbeta.list[[iter]]    <- logbetah
    mix.list[[iter]]    <- zeta
    #ELBO
    elboresult <- ELBO.gam(numStates = S, numMix = M, numLoc = L, stateProb = VBEout$q_tj, mixProb = VBEout$q_tjml, initProb = VBEout$q_1j, jtTransMat = VBEout$q_jk, ct = VBEout$ct, xi = VBMout$xi_j, alpha = VBMout$alpha_j, zeta = VBMout$zeta_jl, gamma_hyper = VBMout$gamma_jml, delta_hyper = VBMout$delta_jml, theta_hyper = VBMout$theta_jml, logbetaprior = logbetah, logbetapost = VBMout$log_beta_jml, exp_omega = omega_exp, exp_psi_omega = omega_psi, exp_lomega = omega_lgamma, obs, h = VBMout$h_jml)
    logbetah <- VBMout$log_beta_jml
    elbo[iter] <- elboresult
    if(iter > 1){elbo_old <- elbo[iter-1]}
    improvement_elbo <- (elboresult - elbo_old)/elbo_old}
  for(l in 1:L){
    params <- post_param.gam(numStates = K,numMix = M, gamma.post = VBMout$gamma_jml[,,l], delta.post = VBMout$delta_jml[,,l], theta.post = VBMout$theta_jml[,,l], exp_alpha = omega_exp[,,l], zeta = VBMout$zeta_jl[,,l], alpha = VBMout$alpha_j, xi = VBMout$xi_j)
    zeta.post[,,l]   <- zeta.post[,,l] + params$MixProb
    lambda.post[,,l] <- lambda.post[,,l] + params$Rainfall_rate
    omega.post[,,l] <- omega.post[,,l] + params$Rainfall_shape}
  pi.post     <- params$InitDist
  tmat.post   <- params$TransMat
  xi.track <- array(0, dim = c(iter, S))
  alpha.track <- array(0, dim = c(S,S,iter))
  mix.track <- array(0, dim = c(S,M,L,iter))
  gamma.track <- array(0, dim = c(S,M - 1,L, iter))
  delta.track <- array(0, dim = c(S,M - 1,L, iter))
  theta.track <- array(0, dim = c(S,M - 1,L, iter))
  logbeta.track <- array(0, dim = c(S,M - 1,L, iter))
  for(i in 1:iter){
    xi.track[i,] <- xi.list[[i]]
    alpha.track[,,i] <- alpha.list[[i]]
    mix.track[,,,i] <- mix.list[[i]]
    gamma.track[,,,i] <- gamma.list[[i]]
    delta.track[,,,i] <- delta.list[[i]]
    theta.track[,,,i] <- theta.list[[i]]
    logbeta.track[,,,i] <- logbeta.list[[i]]
  }
  param.tracker = list('xi' = xi.track, 'alpha' = alpha.track, 'zeta' = mix.track, 'gamma' = gamma.track, 'delta' = delta.track, 'theta'= theta.track, 'logbeta'= logbeta.track)
  posteriors <- list('pi' = pi.post, 'transmat' = tmat.post, 'gamma_hyper' = gammah, 'delta_hyper' = deltah, 'theta_hyper' = thetah, 'beta_hyper' = logbetah, 'mix' = zeta.post, 'constants' = VBMout$h_jml, 'lambda' = lambda.post, 'omega' = omega.post)
  output = list('posteriors' = posteriors, 'ELBO' = elbo, 'iternum' = iter, 'param.tracker' = param.tracker)}

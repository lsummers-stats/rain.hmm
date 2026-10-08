StoVBEM.exp = function(D, S, Y, L, M, xi, alpha, zeta, gamma_shape, gamma_rate, obs, mix.samples = F, maxiter = 100) {
  mix       <- mix.samples
  data      <- array(obs,dim = c(D,Y,L))
  dic       <- rep(0,maxiter)
  dic_old   <- 50000
  dic[1]    <- 25000
  elbo      <- rep(0,maxiter)
  elbo_old  <- -50000
  elbo[1]   <- -25000
  tol       <- 10^(-9)
  iter      <- 1
  improvement_dic   <- (dic_old-dic[1])/dic_old
  improvement_elbo  <- (elbo_old-elbo[1])/elbo_old
  lambda.post = array(0,dim = c(S,M-1,L))
  zeta.post = array(0,dim = c(S,M,L))
  tmat.post = matrix(0,S,S)
  pi.post = rep(0,S)
  #Empty Variables
  a_jk        <- array(0, dim = c(S,S)) # posterior state kernel
  b_tj        <- array(1, dim = c(D,S)) # posterior emission kernel
  b_tjl       <- array(1, dim = c(D,S,L))
  a_1j        <- c(rep(0, times = D)) # posterior initial probability kernel

  logct       <- c(rep(0, times = D))
  fvar        <- array(0, dim = c(D,S))
  bvar        <- array(0, dim = c(D,S))
  b_star      <- array(0, dim = c(D,S,M))
  ct          <- c(rep(0, times = D))

  q_tj        <- array(0, dim = c(D,S)) # posterior probability of stationary distribution
  q_tjml      <- mpfrArray(0, precBits=120, dim = c(D,S,M,L))
  q_jk        <- mpfrArray(0, precBits=120, dim=c(S,S,D-1)) # posterior joint transition probability matrix
  q_1j        <- c(rep(0, times = D)) # posterior initial probability

  xi.list    <- list('prior' = xi)
  alpha.list   <- list('prior' = alpha)
  shape.list   <- list('prior' = gamma_shape)
  rate.list   <- list('prior' = gamma_rate)
  mix.list   <- list('prior' = zeta)

  emptyfillers <- list('a_jk' = a_jk, 'b_tj' = b_tj, 'b_tjl' = b_tjl, 'a_1j' = a_1j, 'ct' = ct, 'fvar' = fvar, 'bvar' = bvar, 'b_star' = b_star, 'q_tj' = q_tj, 'q_tjml' = q_tjml, 'q_jk' = q_jk, 'q_1j' = q_1j)
  #Model Running
  while((abs(improvement_elbo) > tol | improvement_elbo <0) & iter<maxiter){
    y_sample <- Data_Randomizer(numLoc = L, obs = data, mix.sample = mix)
    emptyfillers[["del"]] <- y_sample$del_y0
    ## Inputs need to be changed to accommodate for the Stochastic Change
    VBEout <- VBEmpfr.exp(numDays = D, numStates = S, numLoc = L, numMix = M, xi = xi, alpha = alpha, zeta = zeta, gamma_shape = gamma_shape, gamma_rate = gamma_rate, var = emptyfillers, obs = y_sample$data)
    VBMout <- VBMS.exp(numStates = S, numLoc = L, numYears = Y, numMix = M, xi = xi.list, alpha = alpha.list, zeta = mix.list, gamma_shape = shape.list, gamma_rate = rate.list, q_1j = VBEout$q_1j, q_tj = VBEout$q_tj, q_tjml = VBEout$q_tjml, q_jk=VBEout$q_jk, obs = y_sample$data, iter = iter)
    ## Update iter
    iter = iter+1
    ##Quick update
    xi <- VBMout$xi_j
    alpha <- VBMout$alpha
    zeta <- VBMout$zeta_jl
    gamma_shape <- VBMout$gamma_jml
    gamma_rate <- VBMout$delta_jml
    numericmix <- VBEout$q_tjml
    numerictrans <- VBEout$q_jk
    emptyfillers$a_jk <- VBEout$a_jk
    emptyfillers$b_tj <- VBEout$b_tj
    emptyfillers$b_tjl <- VBEout$b_tjl
    emptyfillers$q_tj <- VBEout$q_tj
    emptyfillers$q_tjml <- VBEout$q_tjml
    emptyfillers$q_jk <- VBEout$q_jk
    emptyfillers$q_1j <- VBEout$q_1j
    xi.list[[iter]]    <- xi
    alpha.list[[iter]]  <- alpha
    shape.list[[iter]]   <- gamma_shape
    rate.list[[iter]]   <- gamma_rate
    mix.list[[iter]]   <- zeta
    #ELBO
    elboresult <- StoELBO(numStates = S, numMix = M, numLoc = L, stateProb = VBEout$q_tj, mixProb = numericmix, initProb = VBEout$q_1j, jtTransMat = numerictrans, ct = VBEout$ct, xi = VBMout$xi_j, alpha = VBMout$alpha, zeta = VBMout$zeta_jl, gamma_shape = VBMout$gamma_jml, gamma_rate = VBMout$delta_jml, obs = y_sample$data, h = VBMout$h_jml)
    elbo[iter] <- elboresult
    elbo_old <- elbo[iter-1]
    #DIC
    DICout <- StoDIC(numStates = S, numMix = M, numLoc = L, stateProb = VBEout$q_tj, mixProb = numericmix, initProb = VBEout$q_1j, jtTransMat = numerictrans, ct = VBEout$ct, xi = VBMout$xi_j, alpha = VBMout$alpha, zeta = VBMout$zeta_jl, gamma_shape = VBMout$gamma_jml, gamma_rate = VBMout$delta_jml)
    dic[iter] <- DICout$dic
    dic_old <- dic[iter-1]
    improvement_elbo <- (elbo_old-elbo[iter])/elbo_old
    improvement_dic <- (dic_old-dic[iter])/dic_old
    }
  for(l in 1:L){
    params <- post_param(numStates = S,numMix = M,gamma.post = VBMout$gamma_jml[,,l],delta.post = VBMout$delta_jml[,,l], zeta = VBMout$zeta_jl[,,l], alpha = VBMout$alpha_j, xi = VBMout$xi_j)
    zeta.post[,,l]   <-zeta.post[,,l] + params$MixProb
    lambda.post[,,l] <-lambda.post[,,l] + params$RainRate}
  xi.track <- array(0, dim = c(iter, S))
  alpha.track <- array(0, dim = c(S,S,iter))
  mix.track <- array(0, dim = c(S,M,L,iter))
  shape.track <- array(0, dim = c(S,M - 1,L, iter))
  rate.track <- array(0, dim = c(S,M - 1,L, iter))
  for(i in 1:iter){
    xi.track[i,] <- xi.list[[i]]
    alpha.track[,,i] <- alpha.list[[i]]
    mix.track[,,,i] <- mix.list[[i]]
    shape.track[,,,i] <- shape.list[[i]]
    rate.track[,,,i] <- rate.list[[i]]
  }
  param.tracker = list('xi' = xi.track, 'alpha' = alpha.track, 'zeta' = mix.track, 'shape' = shape.track, 'rate' = rate.track)
  pi.post     <- params$InitDist
  tmat.post   <- params$TransMat
  posteriors <- list('pi' = pi.post, 'transmat' = tmat.post, 'gamma_shape' = gamma_shape, 'gamma_rate' = gamma_rate, 'mix' = zeta.post, 'constants' = VBMout$h_jml, 'lambda' = lambda.post)
  output = list('posteriors' = posteriors, 'ELBO' = elbo, 'DIC' = dic, 'iternum' = iter, 'param.tracker' = param.tracker)}



#' Stochastic Variational Bayes Expectation Formula
#'
#' @description
#' A function dedicated to the Variational Bayes EM process, simplified using Stochastic VI.
#'
#' @details
#' Used to run the Stochastic Variational Bayes EM Process for a Hidden Markov Model designed
#'  for rain data.
#'
#' @param D Number of Days per Year.
#'
#' @param S Number of States laid out in the model.
#'
#' @param Y Number of Years collected.
#'
#' @param L Number of locations data recorded in the data.
#'
#' @param M Number of mixtures predetermined by user.
#'
#' @param xi The Initial Matrix for the Initial States.
#'
#' The initial matrix of the hyperparameters for the Dirichlet distribution used to describe
#'  the initial probabilities each state.
#'
#' @param alpha The Initial Matrix for the Transitions.
#'
#' The initial matrix of the hyperparameters for the Dirichlet distribution used to describe
#'  the transition probabilities from state to state.
#'
#' @param zeta The Initial Matrix for the Mixtures.
#'
#' The initial matrix of the hyperparameters for the Dirichlet distribution used to describe
#'  the mixture probabilities at each location.
#'
#' @param gamma_shape The Initial matrix for the Shape Parameters.
#'
#' The initial matrix of the shape hyperparameters for each Gamma mixture component.
#'
#' @param gamma_rate The Initial Matrix for the Rate Parameters.
#'
#' The initial matrix of the rate hyperparameters for each Gamma mixture component.
#'
#' @param obs The dataset.
#'
#' @param mix.samples The indicator used to decide the data should come from different years.
#'
#' The initial value is set to FALSE.
#'
#' @return A list of objects used to describe the model:
#'
#'  * `priors`: A list containing all the starting matrices and probabilities.
#'  * `posteriors`: A list containing all the posterior matrices and probabilities.
#'  * `ELBO`: A vector tracking the ELBO as the model progresses (with spot `i` corresponding to iteration `i`)
#'  * `DIC`: A vector tracking the DIC as the model progresses (with spot `i` corresponding to iteration `i`)
#'
#' @export

StoVBEM = function(D, S, Y, L, M, xi, alpha, zeta, gamma_shape, gamma_rate, obs, mix.samples = F) {
  mix       <- mix.samples
  N         <- dim(obs)[1]
  y2        <- array(y,dim = c(D,Y,L))
  maxiter   <- 5  # number of iterations to run the code for
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
  priors    <- list('xi_prior' = xi, 'alpha_prior' = alpha, 'zeta_prior' = zeta, 'gamma_shape_prior' = gamma_shape, 'gamma_rate_prior' = gamma_rate)
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
  b_star      <- mpfrArray(0, precBits=120, dim = c(D,S,M))
  ct          <- c(rep(0, times = D))

  q_tj        <- array(0, dim = c(D,S)) # posterior probability of stationary distribution
  q_tjml      <- mpfrArray(0, precBits=120, dim = c(D,S,M,L))
  q_jk        <- mpfrArray(0, precBits=120, dim=c(S,S,D-1)) # posterior joint transition probability matrix
  q_1j        <- c(rep(0, times = D)) # posterior initial probability

  emptyfillers <- list('a_jk' = a_jk, 'b_tj' = b_tj, 'b_tjl' = b_tjl, 'a_1j' = a_1j, 'ct' = ct, 'fvar' = fvar, 'bvar' = bvar, 'b_star' = b_star, 'q_tj' = q_tj, 'q_tjml' = q_tjml, 'q_jk' = q_jk, 'q_1j' = q_1j)
  #Model Running
  while((abs(improvement_elbo) > tol | improvement_elbo <0) & iter<maxiter){
    y_sample <- Data_Randomizer(numLoc = L, obs = y2, mix.sample = mix)
    emptyfillers[["del"]] <- y_sample$del_y0
    ## Inputs need to be changed to accommodate for the Stochastic Change
    VBEout <- VBEmpfr.exp(numDays = D, numStates = S, numLoc = L, numMix = M, xi = xi, alpha = alpha, zeta = zeta, gamma_shape = gamma_shape, gamma_rate = gamma_rate, var = emptyfillers, obs = y_sample$data)
    VBMout <- VBMS.exp(numStates = S, numLoc = L, numYears = Y, numMix = M, xi = xi, alpha = alpha, zeta = zeta, gamma_shape = gamma_shape, gamma_rate = gamma_rate, q_1j = VBEout$q_1j, q_tj = VBEout$q_tj, q_tjml = VBEout$q_tjml, q_jk=VBEout$q_jk, obs = y_sample$data, iter = iter, priors = priors)
    ## Update iter
    iter = iter+1
    ##Quick update
    xi <- VBMout$xi_j
    alpha <- VBMout$alpha
    zeta <- VBMout$zeta_jl
    gamma_shape <- VBMout$gamma_jml
    gamma_rate <- VBMout$delta_jml
    #ELBO
    elboresult <- ELBO(numStates = S, numMix = M, numLoc = L, stateProb = VBEout$q_tj, mixProb = VBEout$q_tjml, initProb = VBEout$q_1j, jtTransMat = VBEout$q_jk, ct = VBEout$ct, xi = VBMout$xi_j, alpha = VBMout$alpha, zeta = VBMout$zeta_jl, gamma_shape = VBMout$gamma_jml, gamma_rate = VBMout$delta_jml, y_sample$data, h = VBMout$h_jml)
    elbo[iter] <- elboresult
    elbo_old <- elbo[iter-1]
    #DIC
    DICout <- DIC(numStates = S, numMix = M, numLoc = L, stateProb = VBEout$q_tj, mixProb = VBEout$q_tjml, initProb = VBEout$q_1j, jtTransMat = VBEout$q_jk, ct = VBEout$ct, xi = VBMout$xi_j, alpha = VBMout$alpha, zeta = VBMout$zeta_jl, gamma_shape = VBMout$gamma_jml, gamma_rate = VBMout$delta_jml)
    dic[iter] <- DICout$dic
    dic_old <- dic[iter-1]
    improvement_elbo <- (elbo_old-elbo[iter])/elbo_old
    improvement_dic <- (dic_old-dic[iter])/dic_old}
  for(l in 1:L){
    params <- post_param(numStates = K,numMix = M,gamma.post = VBMout$gamma_jml[,,l],delta.post = VBMout$delta_jml[,,l], zeta = VBMout$zeta_jl[,,l], alpha = VBMout$alpha_j, xi = VBMout$xi_j)
    zeta.post[,,l]   <-zeta.post[,,l] + params$MixProb
    lambda.post[,,l] <-lambda.post[,,l] + params$RainRate}
  pi.post     <- params$InitDist
  tmat.post   <- params$TransMat
  posteriors <- list('pi' = pi.post, 'transmat' = tmat.post, 'gamma_shape' = gamma_shape, 'gamma_rate' = gamma_rate, 'mix' = zeta.post, 'constants' = VBMout$h_jml, 'lambda' = lambda.post)
  output = list('priors' = priors, 'posteriors' = posteriors, 'ELBO' = elbo, 'DIC' = dic)
}



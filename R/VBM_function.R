#' Variational Bayes Maximization Formula
#'
#' @description
#' A function dedicated to the maximization step in a Variational Bayes EM algorithm
#'
#' @details
#' Used as part of the Variational Bayes EM Process for a Hidden Markov Model designed
#' for rain data.
#'
#' @param numStates Number of States laid out in the model
#'
#' @param numLoc Number of locations data recorded in the data
#'
#' @param numMix Number of mixtures predetermined by user
#'
#' @param xi The Current Matrix for the Initial States.
#'
#' The current matrix of the hyperparameters for the Dirichlet distribution used to describe
#' the initial probabilities each state.
#'
#' @param alpha The Current Matrix for the Transitions.
#'
#' The current matrix of the hyperparameters for the Dirichlet distribution used to describe
#' the transition probabilities from state to state.
#'
#' @param zeta The Current Matrix for the Mixtures.
#'
#' The current matrix of the hyperparameters for the Dirichlet distribution used to describe
#' the mixture probabilities at each location.
#'
#' @param gamma_shape The Current matrix for the Shape Parameters.
#'
#' The current matrix of the shape hyperparameters for each Gamma mixture component.
#'
#' @param gamma_rate The Current Matrix for the Rate Parameters.
#'
#' The current matrix of the rate hyperparameters for each Gamma mixture component.
#'
#' @param q_1j Posterior initial probability.
#'
#' This matrix is usually the result of the VBE.exp function.
#'
#' @param q_tj Posterior probability of stationary distribution.
#'
#' This matrix is usually the result of the VBE.exp function.
#'
#' @param q_tjml An augmented `q_tjml` array.
#'
#' This matrix is usually the result of the VBE.exp function.
#'
#' @param q_jk Posterior joint transition probability matrix.
#'
#' This matrix is usually the result of the VBE.exp function.
#'
#' @param obs The vector of observations.
#'
#' The vector containing the y-values (typically the precipitation amounts) from the data.
#'
#' @param iter The current iteration step.
#'
#' @return A list of objects built from empty objects in `var`:
#'
#'  * `gamma_jml`: Posterior shape of exponential rate.
#'  * `delta_jml`: Posterior rate of exponential rate.
#'  * `xi_j`: Posterior Dirichlet paramters for initial distribution.
#'  * `alpha_j`: Posterior Dirichlet parameters for transition matrix rows.
#'  * `zeta_jl`: Posterior Dirichlet parameters for mixing probabilities.
#'  * `h_jml`: A matrix of coefficients for normalizing the mixture distributions.
#'
#'  @export


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

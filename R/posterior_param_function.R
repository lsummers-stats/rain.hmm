#' Posterior Parameter Function
#'
#' @description
#' A function used to put the posterior distributions of important parameters into a single list.
#'
#' @details
#' A function used to convert different datasets from an EM algorithm step into a single
#' consolidated list, while normalizing any variables that need to be normalized. Used at the very end of
#' an EM Algorithm step for a Variational Bayes Hidden Markov Model. This function is particular to using
#' Exponential Priors for each mixture.
#'
#' @param numStates Number of States
#'
#' The total number of states planned for the Hidden Markov Model.
#'
#' @param numMix Number of Mixtures
#'
#' The total number of mixtures planned for the Hidden Markov Model.
#'
#' @param gamma.post The posterior matrix for the Shape Parameters
#'
#' The posterior matrix of the shape hyperparameters for each Gamma mixture component.
#'
#' @param delta.post The Posterior Matrix for the Rate Parameters
#'
#' The posterior matrix of the rate hyperparameters for each Gamma mixture component.
#'
#' @param zeta The Posterior Matrix for the Mixtures
#'
#' The posterior matrix of the hyperparameters for the Dirichlet distribution used to describe
#' the mixture probabilities at each location.
#'
#' @param alpha The Posterior Matrix for the Transitions.
#'
#' The posterior matrix of the hyperparameters for the Dirichlet distribution used to describe
#' the transition probabilities from state to state.
#'
#' @param xi The Posterior Matrix for the Initial States.
#'
#' The posterior matrix of the hyperparameters for the Dirichlet distribution used to describe
#' the initial probabilities each state.
#'
#' @returns A list holding four objects:
#'  * `InitDist` : A vector containing the posterior distribution of `pi_1`, the initial state probabilities.
#'  * `TransMat` : A matrix containing the posterior distribution of `A`, the transition state probabilities.
#'  * `MixProb` : A matrix containing the posterior distribution of `C`, the cluster mixture probabilities.
#'  * `RainRate` : A matrix containing the posterior distribution of `lambda`, the parameter of the exponential distribution
#'  for rainfall at any given location.
#'
#' @examples
#' K = 2
#' M = 2
#' gamma = c(1,2)
#' delta = c(1,4)
#' zeta = matrix(c(1,1,1,1), nrow = K, ncol = M)
#' alpha = matrix(c(1,0,0,1), nrow = K, ncol = M)
#' xi = c(1,1,1,1)
#' post_var <- post_param(K, M, gamma, delta, zeta, alpha, xi)
#' print(post_var$RainRate)
#'
#' @export


post_param <- function(numStates, numMix, gamma.post, delta.post, zeta, alpha, xi){
  tmat.post <- matrix(nrow = numStates,ncol = numStates)
  lambda.post <- gamma.post/delta.post
  zeta.post <- matrix(nrow = numStates,ncol = numMix)
  for(j in 1:numStates){
    zeta.post[j,] <- zeta[j,]/sum(zeta[j,])
    tmat.post[j,] <- alpha[j,]/sum(alpha[j,])
  }
  pi.post <- xi/sum(xi)

  output <- list('InitDist' = pi.post, 'TransMat' = tmat.post, 'MixProb' = zeta.post, 'RainRate' = lambda.post)
}

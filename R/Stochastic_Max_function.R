#' Stochastic Maximization Formula
#'
#' @description
#' A function dedicated to the maximization step in a Variational Bayes EM algorithm, simplified using Stochastic VI.
#'
#' @details
#' Used as part of the Stochastic Variational Bayes EM Process for a Hidden Markov Model designed
#'  for rain data.
#'
#' @param numStates Number of States laid out in the model.
#'
#' @param numLoc Number of locations data recorded in the data.
#'
#' @param numMix Number of mixtures predetermined by user.
#'
#' @param numYears Number of Years collected.
#'
#' @param numDays Number of Days per Year.
#'
#' @param numStates Number of States laid out in the model.
#'
#' @param xi The Current Matrix for the Initial States.
#'
#' The current matrix of the hyperparameters for the Dirichlet distribution used to describe
#'  the initial probabilities each state.
#'
#' @param alpha The Current Matrix for the Transitions.
#'
#' The current matrix of the hyperparameters for the Dirichlet distribution used to describe
#'  the transition probabilities from state to state.
#'
#' @param zeta The Current Matrix for the Mixtures.
#'
#' The current matrix of the hyperparameters for the Dirichlet distribution used to describe
#'  the mixture probabilities at each location.
#'
#' @param gamma_shape The Current matrix for the Shape Parameters.
#'
#' The current matrix of the shape hyperparameters for each Gamma mixture component.
#'
#' @param gamma_rate The Current Matrix for the Rate Parameters.
#'
#' The current matrix of the rate hyperparameters for each Gamma mixture component.
#'
#' @param q_1j Expectation of latent variable s_1j.
#'
#' The latent variable s_tj represents if, at time t, the state is j.
#'  s_1j is the initial states.
#'
#' @param q_tj Expectation of latent variable s_tj.
#'
#' The latent variable s_tj represents if, at time t, the state is j.
#'
#' @param q_tjml Expectation of latent variable r_tjml.
#'
#' The latent variable r_tjml represents an indicator variable with 1 if y_tl (the data point at time t and location l)
#'  comes from mixture m, and s_t = j.
#'
#' @param q_jk Expectation of state transistions
#'
#' The latent variables s_tj*s_t+1,k represent the movement from state j to k from time t to t+1.
#'
#' @param obs The dataset.
#'
#' @param iter The current step of the iteration process.
#'
#' @param priors A list containing the priors used at the start of the process.
#'
#' The priors are needed for xi, alpha, zeta, gamma_shape, and gamma_rate.
#'
#' @export
#'
#' @return A list of objects:
#'
#'  * `gamma_jml`: Posterior shape of exponential rate.
#'  * `delta_jml`: Posterior rate of exponential rate.
#'  * `xi_j`: Posterior Dirichlet paramters for initial distribution.
#'  * `alpha_j`: Posterior Dirichlet parameters for transition matrix rows.
#'  * `zeta_jl`: Posterior Dirichlet parameters for mixing probabilities.
#'  * `h_jml`: A matrix of coefficients for normalizing the mixture distributions.
#'
#' @details
#'
#' Stochastic VI assumes exchangeablility for the emission distributions used for each year. Thus
#'  rather than trying to compute N data points, we split the data with N = D*Y, days and years, then account for the number of years
#'  by multiplying parts of the update by Y. However, the user must decide some step size tau. The step size must
#'  fulfill the Robbins-Monro conditions.
#'
#'  In particular, this function uses stochastic gradient ascent with step size tau = 1/iter.
#'

VBMS.exp = function(numStates, numLoc, numMix, numYears, numDays, xi, alpha, zeta, gamma_shape, gamma_rate, q_1j, q_tj, q_tjml, q_jk, obs, iter, priors){
  #Q's will be randomized already, so no need to edit the sums, just multiply by N
  N = numYears
  step = 1/iter
  gamma_jml    <- gamma_shape # posterior shape of exponential rate
  delta_jml    <- gamma_rate # posterior rate of exponential rate
  zeta_jl      <- zeta # posterior Dirichlet parameters for mixing probabilities
  alpha_j      <- alpha # posterior Dirichlet parameters for transition matrix rows
  xi_j         <- xi # posterior Dirichlet parameters for initial distribution
  gamma_prior      <- priors$gamma_shape_prior
  delta_prior      <- priors$gamma_rate_prior
  xi_prior         <- priors$xi_prior
  alpha_prior      <- priors$alpha_prior
  zeta_prior      <- priors$zeta_prior
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

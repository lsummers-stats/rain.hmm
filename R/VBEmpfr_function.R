#' High Precision Variational Bayes Expectation Formula
#'
#' @description
#' A function dedicated to the expectation step in a Variational Bayes EM algorithm with higher precision.
#'
#' @details
#' Used as part of the Variational Bayes EM Process for a Hidden Markov Model designed
#'  for rain data.
#'
#' @param numDays Number of Days per Year.
#'
#' @param numStates Number of States laid out in the model.
#'
#' @param numYears Number of Years collected.
#'
#' @param numLoc Number of locations data recorded in the data.
#'
#' @param numMix Number of mixtures predetermined by user.
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
#' @param gamma.post The Current matrix for the Shape Parameters.
#'
#' The current matrix of the shape hyperparameters for each Gamma mixture component.
#'
#' @param delta.post The Current Matrix for the Rate Parameters.
#'
#' The current matrix of the rate hyperparameters for each Gamma mixture component.
#'
#' @param var A list of matrices and vectors objects:
#' * `del_y0`: Rain indicator function. Typically an `ifelse` function.
#' * `a_jk`: Posterior state kernel, a matrix filled with 0's with dimensions: `numStates` by `numStates`.
#' * `b_tj`: Posterior emission kernel, an array filled with 1's with dimensions: `numDays` by `numStates` by `numYears`.
#' * `b_tjl`: An augmented `b_tj` array, filled with 0's with dimensions: `numDays` by `numStates` by `numYears` by `numLoc`.
#' * `a_1j`: Posterior initial probability kernel, a matrix with dimensions: `numYears` by `numStates`.
#' * `ct`: A matrix filled with 0's with dimensions: `numDays` by `numYears`.
#' * `fvar`: An array filled with 0's with dimensions: `numDays` by `numStates` by `numYears`.
#' * `bvar`: An array filled with 0's with dimensions: `numDays` by `numStates` by `numYears`.
#' * `b_star`: An array filled with 0's with dimensions: `numDays` by `numStates` by `numMix`.
#' * `q_tj`: Posterior probability of stationary distribution, an array filled with 0's with dimensions: `numDays` by `numStates` by `numYears`.
#' * `q_tjml`: An augmented `q_tjml` array, filled with 0's with dimensions: `numDays` by `numStates` by `numMix` by `numStates` by `numLoc`.
#' * `q_jk`: Posterior joint transition probability matrix, an array filled with 0's with dimensions: `numStates` by `numStates` by `numDays` - 1 by `numYears`.
#' * `q_1j`: Posterior initial probability, a matrix filled with 0's with dimensions: `numYears` by `numStates`.
#'
#' @param obs The vector of observations.
#'
#' The vector containing the y-values (typically the precipitation amounts) from the data.
#'
#' @return A list of objects built from empty objects in `var`:
#'  * `a_jk`: The posterior state kernel.
#'  * `b_tj`: Posterior emission kernel.
#'  * `c_t` : A series of constants used to normalize the forward variables.
#'  * `q_1j`: Posterior initial probability.
#'  * `q_tj`: Posterior probability of stationary distribution.
#'  * `q_tjml`: An augmented `q_tjml` array.
#'  * `q_jk`: Posterior joint transition probability matrix.
#'
#' @export


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

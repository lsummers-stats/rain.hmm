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

post_param.gam <- function(numStates, numMix, gamma.post, delta.post, theta.post, exp_alpha, zeta, alpha, xi){
  tmat.post <- matrix(nrow = numStates,ncol = numStates)
  omega.post <- exp_alpha
  lambda.post <- omega.post*gamma.post/theta.post
  zeta.post <- matrix(nrow = numStates,ncol = numMix)
  for(j in 1:numStates){
    zeta.post[j,] <- zeta[j,]/sum(zeta[j,])
    tmat.post[j,] <- alpha[j,]/sum(alpha[j,])
  }
  pi.post <- xi/sum(xi)

  output <- list('InitDist' = pi.post, 'TransMat' = tmat.post, 'MixProb' = zeta.post, 'Rainfall_rate' = lambda.post, 'Rainfall_shape' = omega.post)}

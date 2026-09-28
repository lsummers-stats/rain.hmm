#' Omega Normalization Function
#'
#' @description
#' A function designed to find the normalizing constant of the marginal distribution of a shape parameter in
#' a GC2 distribution.
#'
#' @details
#' In Variational Bayes, when dealing with a Gamma distribution, having both the shape and the rate parameter unknown
#' forces a new conjugate prior known as a GC2 distribution. The marginal of the rate parameter is another gamma, but
#' the marginal of the omega is a distribution that requires numerical integration. This function is made to find the
#' constant that normalizes this distribution.
#'
#' @param gamh Value of the gamma hyperparameter of the GC2 distribution.
#'
#' @param delh Value of the delta hyperparameter of the GC2 distribution.
#'
#' @param thetah Value of the theta hyperparameter of the GC2 distribution.
#'
#' @param logbetah Value of the log of beta hyperparameter of the GC2 distribution.
#'
#' @return A single value representing the normalizing constant. Note that, to use in a function,
#' it would need to be 1/result.
#'
#'
#' @export
#'
#'

omega_constant <- function(gamh, delh, thetah, logbetah){
  fun <- function(x){
    out <- exp(lgamma(gamh*x)-delh*(lgamma(x))+(x)*logbetah - (gamh*x)*log(thetah))
  }
  con <- integrate(fun,0,Inf)
  return(con)
}

#' Omega Expectation Function
#'
#' @description
#' A function designed to find the expectation of the shape parameter in a GC2 distribution.
#'
#' @details
#' In Variational Bayes, when dealing with a Gamma distribution, having both the shape and the rate parameter unknown
#' forces a new conjugate prior known as a GC2 distribution. The marginal of the rate parameter is another gamma, but
#' the marginal of the omega is a distribution that requires numerical integration. This function is made to find the
#' expectation of the shape parameter in a GC2 distribution.
#'
#' @param gamh Value of the gamma hyperparameter of the GC2 distribution.
#'
#' @param delh Value of the delta hyperparameter of the GC2 distribution.
#'
#' @param thetah Value of the theta hyperparameter of the GC2 distribution.
#'
#' @param logbetah Value of the log of beta hyperparameter of the GC2 distribution.
#'
#' @param con The normalizing constant, typically found using the `omega_constant` function.
#'
#' @return A single value representing the expectation of the shape parameter.
#'
#'
#' @export

exp_omega <- function(gamh, delh, thetah, logbetah, con){
  fun <- function(x){
    out <- exp(log(x) + lgamma(gamh*x)-delh*(lgamma(x))+(x)*logbetah - (gamh*x)*log(thetah))
  }
  result <- integrate(fun,0,Inf, stop.on.error = FALSE)$value
  return(result/con)
}

#' Omega Log-Gamma Function
#'
#' @description
#' A function designed to find the expectation of the log-gamma of the shape parameter in a GC2 distribution.
#'
#' @details
#' In Variational Bayes, when dealing with a Gamma distribution, having both the shape and the rate parameter unknown
#' forces a new conjugate prior known as a GC2 distribution. The marginal of the rate parameter is another gamma, but
#' the marginal of the omega is a distribution that requires numerical integration. This function is made to find the
#' log-gamma of the shape parameter in a GC2 distribution.
#'
#' @param gamh Value of the gamma hyperparameter of the GC2 distribution.
#'
#' @param delh Value of the delta hyperparameter of the GC2 distribution.
#'
#' @param thetah Value of the theta hyperparameter of the GC2 distribution.
#'
#' @param logbetah Value of the log of beta hyperparameter of the GC2 distribution.
#'
#' @param con The normalizing constant, typically found using the `omega_constant` function.
#'
#' @return A single value representing the expectation of the log-gamma of the shape parameter.
#'
#'
#' @export

exp_l_omega <- function(gamh, delh, thetah, logbetah, con){
  fun <- function(x){
    out <- lgamma(x)*exp(lgamma(gamh*x)-delh*(lgamma(x))+(x)*logbetah - (gamh*x)*log(thetah))
  }
  result <- integrate(fun,0,Inf, stop.on.error = FALSE)$value
  return(result/con)
}

#' Omega Psi-Gamma Function
#'
#' @description
#' A function designed to find a special expression of shape parameter in a GC2 distribution.
#'
#' @details
#' In Variational Bayes, when dealing with a Gamma distribution, having both the shape and the rate parameter unknown
#' forces a new conjugate prior known as a GC2 distribution. The marginal of the rate parameter is another gamma, but
#' the marginal of the omega is a distribution that requires numerical integration. The expression E(omega times Psi(omega * gammah))
#' appears when calculating the ELBO (in particular the KL Divergence of the GC2 distribution) and the expectation step in
#' an EM algorithm. Omegga is the shape parameter, gamma is the gamma hyperparameter of the GC2 distribution and Psi refers to the
#' derivative of the log of the gamma distribution.
#'
#' @param gamh Value of the gamma hyperparameter of the GC2 distribution.
#'
#' @param delh Value of the delta hyperparameter of the GC2 distribution.
#'
#' @param thetah Value of the theta hyperparameter of the GC2 distribution.
#'
#' @param logbetah Value of the log of beta hyperparameter of the GC2 distribution.
#'
#' @param con The normalizing constant, typically found using the `omega_constant` function.
#'
#' @return A single value representing special expectation of the shape parameter.
#'
#'
#' @export

exp_psi_omega <- function(gamh, delh, thetah, logbetah, con){
  fun <- function(x){
    out <- psigamma(x * gamh)*exp(lgamma(x) -delh*(lgamma(x))+(x)*logbetah - (gamh*x)*log(thetah))
  }
  onetoinf <- integrate(fun,1,Inf, stop.on.error = FALSE)$value
  zerotoone <- integrate(fun,0,1, stop.on.error = FALSE)$value
  result <- zerotoone + onetoinf
  return(result/con)
}

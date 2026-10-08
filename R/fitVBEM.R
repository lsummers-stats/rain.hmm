#' Variational Bayes EM Algorithm for Rainfall
#'
#' @description
#'  A function dedicated to running the Variational Bayes EM algorithm with several different types of models.
#'
#' @param data The recorded rainfall amounts. The amount must be in a matrix/array with dimensions L by T.
#'
#' @param D Number of Days per Year.
#'
#' @param S Number of States to put in the model.
#'
#' @param Y Number of Years collected.
#'
#'  Default value is set to 1.
#'
#' @param L Number of locations data recorded in the data.
#'
#' @param M Number of mixtures predetermined by user.
#'
#'  Default value is set to 2, 1 mixture being the 0 rainfall mixture, the other being a single distribution of a
#'  specified type.
#'
#' @param dist The preferred distribution to describe rainfall.
#'
#'  As of now, there are only three options, "exp" for an exponential distribution, "stoch.exp" for a stochastic variational bayes of
#'  an exponential distribution and "gamma" for a gamma distribution.
#'
#' @param maxiter The maximum number of iterations the algorithm will do before stopping.
#'
#'  The default value is 1000.
#'
#' @param hypers A list containing all necessary hyperparameters in array form.
#'
#'  For any distribution choice, one must provide hyperparameters for xi, alpha and zeta. For exponential distributions,
#'  the gamma shape and rate must be provided. For gamma distributions, four hyperparameters must all be provided (gamma, delta, theta, logbeta).
#'  See "" for explanations. Please provide the hyperparameters in a list.
#'
#' @param mix.samples A boolean that is only used if "stoch.exp" is chosen.
#'
#'  In Stochastic Variational Bayes, rather than using all data points, the data is broken into years, each with an equal
#'  number of days. Using an exchangeability assumption, the model then only runs the model on a subset of the data. Keeping
#'  mix.samples as `FALSE` means the data used is all from one year, changing it to `TRUE` means that for each day, the model
#'  will randomly pick a year to sample from.
#'
#' @return A list of objects used to describe the model:
#'  * `posteriors`: A list containing all the posterior matrices and probabilities.
#'  * `ELBO`: A vector tracking the ELBO as the model progresses (with spot `i` corresponding to iteration `i`)
#'  * `DIC`: A vector tracking the DIC as the model progresses (with spot `i` corresponding to iteration `i`)
#'  * `iternum`: An integer denoting how many iterations took place.
#'
#' @export

fit.VBEM <- function(data, L, S, M = 2, D, Y = 1, dist = c("exp", "gamma", "stoch.exp"), hypers, mix.samples = F, maxiter = 1000){
  if(dist == "exp" & stochastic == F){
    VBEM.exp(D = D, S = S, Y = Y, L = L, M = M, xi = hypers$xi, alpha = hypers$alpha, zeta = hypers$zeta,
          gamma_shape = hypers$gamma_shape, gamma_rate = hypers$gamma_rate, obs = data, maxiter = maxiter)
  }
  if(dist == "stoch.exp"){
    StoVBEM.exp(D = D, S = S, Y = Y, L = L, M = M, xi = hypers$xi, alpha = hypers$alpha, zeta = hypers$zeta,
         gamma_shape = hypers$gamma_shape, gamma_rate = hypers$gamma_rate, obs = data, maxiter = maxiter, mix.samples = F)
  }
  if(dist == "gamma" & M > 2){
    print("The code for multiple gamma mixtures is unstable, and will likely result in errors. Please keep M = 2 for now.")
  }
  if(dist == "gamma"){
    VBEM.gam(D = D, S = S, Y = Y, L = L, M = M, xi = hypers$xi, alpha = hypers$alpha, zeta = hypers$zeta,
             gammah = hypers$gamma, deltah= hypers$delta, thetah = hypers$theta, logbetah= hypers$logbeta,
             obs = data, maxiter = maxiter)
  }
  print("please select a model type: exp for exponential, stoch.exp for stochastic exponential or gamma for gamma")
}

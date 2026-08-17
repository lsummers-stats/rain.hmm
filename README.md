HMMVI: Hidden Markov Models with Variational Bayes
================

This package is designed…

# Variational Bayes Basic Function

The basic function uses Variational Bayes assuming exponential mixture
components:

``` r
library(HMMVI)
#Pulling the Fake Data Set Faux 1
set.seed(1)
y <- Faux1[[1]]

#The VBEM function requires you to start with some potential priors, as well as pre-pick the number of days, years, mixtures, and number of locations you collected data from.

##Assuming 3 states, 3 mixtures and 3 locations, with 1800 days of data collected and not split into years

S = 3
M = 3
L = 3

##Initial estimates of the hyperparameters:
gamma_0     <- matrix(data = c(.5,2,1.5,9,2,16),nrow = S,ncol = M-1,byrow = T) # shape of exponential rate
gamma_0     <- array(gamma_0,dim = c(S,M-1,L))
delta_0     <- matrix(data = c(2,2,2,2,2,2),nrow = S,ncol = M-1,byrow = T) # rate of exponential rate
delta_0     <- array(delta_0,dim = c(S,M-1,L))
zeta_0      <- matrix(data = c(6,8,6,6,7,7,8,6,6),nrow = S,ncol = M,byrow = T)/(M-1) # Dirichlet prior parameters for mixing probabilities
zeta_0      <- array(zeta_0,dim = c(S,M,L))
alpha_0     <- matrix(rep(10,S^2),nrow = S,byrow = T)/S # Dirichlet prior for transition matrix rows
xi_0        <- rep(1,S)/S # Dirichlet prior for initial probabilities
h_j0        <- gamma_0*log(delta_0) - lgamma(gamma_0) # constant terms in prior (log)

Model <- VBEM(D = 1800, S = S, Y = 1, L = 3, M = 3, xi = xi_0, alpha = alpha_0, zeta = zeta_0, gamma_shape = gamma_0, gamma_rate = delta_0, y)
```

The output is a list, containing information such as the posterior
hyperparameters, or the ELBO/DIC.

``` r
print(Model$posteriors$pi)
```

    ## [1] 0.0006409779 0.0032459527 0.9961130694

``` r
Model$posteriors$transmat
```

    ##           [,1]      [,2]      [,3]
    ## [1,] 0.6383142 0.2254607 0.1362250
    ## [2,] 0.1998351 0.3932420 0.4069229
    ## [3,] 0.2188876 0.3725002 0.4086122

## Post-Formulation Functions

This package includes additional functions like the viterbi-encoding
algorithm, so it does not have to be programmed on the side:

``` r
#Viterbi Here
```

# Stochastic Version

As data gets more complex, with more locations and possibly adding in
more mixtures and states, it becomes necessary to create a less
intensive method. This is the purpose of the stochastic version of the
function.

``` r
library(HMMVI)
library(Rmpfr)
```

    ## Warning: package 'Rmpfr' was built under R version 4.5.3

    ## Loading required package: gmp

    ## Warning: package 'gmp' was built under R version 4.5.3

    ## 
    ## Attaching package: 'gmp'

    ## The following objects are masked from 'package:base':
    ## 
    ##     %*%, apply, crossprod, matrix, tcrossprod

    ## C code of R package 'Rmpfr': GMP using 64 bits per limb

    ## 
    ## Attaching package: 'Rmpfr'

    ## The following object is masked from 'package:gmp':
    ## 
    ##     outer

    ## The following objects are masked from 'package:stats':
    ## 
    ##     dbinom, dchisq, dgamma, dnbinom, dnorm, dpois, dt, pgamma, pnorm

    ## The following objects are masked from 'package:base':
    ## 
    ##     cbind, pmax, pmin, rbind

``` r
#Define some hyperparameter guesses
gamma_0     <- matrix(data = c(.5,2,1.5,5,2,10),nrow = 3,ncol = 2,byrow = T) # shape of exponential rate
gamma_0     <- array(gamma_0,dim = c(3,2,10))
delta_0     <- matrix(data = c(2,2,2,2,2,2),nrow = 3,ncol = 2,byrow = T) # rate of exponential rate
delta_0     <- array(delta_0,dim = c(3,2,10))
zeta_0      <- matrix(data = c(3,4,3,3,3.5,3.5,4,3,3),nrow = 3,ncol = 3,byrow = T) # Dirichlet prior parameters for mixing probabilities
zeta_0      <- array(zeta_0,dim = c(3,3,10))
alpha_0     <- matrix(rep(10,9),nrow = 3,byrow = T)/3 # Dirichlet prior for transition matrix rows
xi_0        <- c(.4,.3,.3) # Dirichlet prior for initial probabilities

StoModel <- StoVBEM(D = 92, S = 3, Y = 20, L = 10, M = 3, xi = xi_0, alpha = alpha_0, zeta = zeta_0, gamma_shape = gamma_0, gamma_rate = delta_0, obs = Faux2, mix.samples = F)
```

``` r
print(StoModel$posteriors$pi)
```

    ## [1] 0.6104147 0.1947920 0.1947933

``` r
print(StoModel$iternum)
```

    ## [1] 10

# Gamma Mixtures

#' A Fake Data set from simulating a HMM
#'
#' A dataset containing the rrainfall amounts and the state associated with said rainfall amount. The data simulates 36
#' days over 10 years, 3 mixtures, 4 hidden states, and 5 locations.
#'
#' @format A list containing:
#' \describe{
#'   \item{rainfall}{The fake rainfall amount recorded over 3600 days at five different locations. Contains several
#'   zero amounts.}
#'   \item{hidden}{The hidden state associated with each rainfall calculation.}
#'   \item{init_tmat}{The transition matrix used in generating the data.}
#'   \item{init_pi}{The starting state probability matrix used in generating the data.}
#'   \item{init_lambda}{The starting lambdas used in generating the rainfall amounts.}
#'   \item{init_zeta}{The starting probabilities for the mixtures used in generating the rainfall amounts.}
#' }
#'
#' @source Generated using a fake HMM mode.
"Faux3"

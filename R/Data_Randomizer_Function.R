#' Data Randomizer Function
#'
#' @description
#'
#' A function used to help sample during Stochastic VI.
#'
#' @details
#'
#' In Stochastic VI, we split the data length N into D*Y, days and years. There are then two options on how to pick
#'  data. Either, one can pick from the same year (all D data points come from the same year) or the data points
#'  can come from all different years (D data points pick from different years).
#'
#' @param obs The dataset.
#'
#' @param mix.samples The indicator used to decide the data should come from different years.
#'
#' Initial set to FALSE.
#'
#' @param numLoc The number of locations.
#'
#' @return A list of objects:
#'
#'  * `y_sample` : The resulting dataset with dimensions days by locations.
#'  * `del_y0_sample` : The indicator data set that uses 0 if no rain, 1 if there is any rain amount.
#'
#' @export
#'

Data_Randomizer <- function(numLoc, obs, mix.sample = F){
  minibatch.size = dim(obs)[1]
  total.batches = dim(obs)[2]
  y_sample = array(dim = c(minibatch.size, numLoc))
  if(mix.sample==T){
    years = sample(1:total.batches,minibatch.size,replace = T)
  } else {
    years = rep(sample(1:total.batches,1), minibatch.size)}
  for(day in 1:minibatch.size){
    y_sample[day,] = y2[day,years[day],]}
  del_y0_sample = ifelse(y_sample==0,1,0)
  output = list('data' = y_sample, 'del_y0' = del_y0_sample)
}

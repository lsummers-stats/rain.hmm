Data_Randomizer <- function(numLoc, obs, mix.sample = F){
  minibatch.size = dim(obs)[1]
  total.batches = dim(obs)[2]
  y_sample = array(dim = c(minibatch.size, numLoc))
  if(mix.sample==T){
    years = sample(1:total.batches,minibatch.size,replace = T)
  } else {
    years = rep(sample(1:total.batches,1), minibatch.size)}
  for(day in 1:minibatch.size){
    y_sample[day,] = obs[day,years[day],]}
  del_y0_sample = ifelse(y_sample==0,1,0)
  output = list('data' = y_sample, 'del_y0' = del_y0_sample)
}

#' @import Rmpfr

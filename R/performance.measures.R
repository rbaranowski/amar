#' Hit rate
#' @export 
#' 
hit.rate <- function(x.true, x.predicted, ...){
  
  #TODO verify parameters
  
  tmp <- x.true != 0 
  
  return(mean(sign(x.true[tmp]) == sign(x.predicted[tmp])))
  
  
}

#' MSE
#' @export 
 

mse <- function(x.true, x.predicted, ...){
  
  #TODO verify parameters
  
  return(mean((x.true-x.predicted)^2))
  
  
}


#' R squared
#' @export 
#' 


r.squared <- function(x.true, x.predicted,  ...){

  #TODO verify parameters
  return(1- mean(abs(x.true-x.predicted)^2) / mean(abs(x.true)^2))



}

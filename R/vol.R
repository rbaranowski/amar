#' Remove volatility
#' 
#' @export 

remove.volatility <- function(returns, 
                              train.id, 
                              criterion = c("jarque.bera"), ...){
  
  #TODO optimise parameters
  criterion <- match.arg(criterion, c("jarque.bera"))
  
  optim.fn <- function(x){
    
    residuals <- vol.exp.smoothing(returns, x)[train.id, "residuals"]
    jarque.bera.test(residuals)$statistic
    
  }
  
  tmp <- optimize(optim.fn, c(0,1))
  tmp2 <- vol.exp.smoothing(returns, tmp$minimum)
  
  return(list(residuals = tmp2[,2], volatility = tmp2[,1], min.lambda=tmp$minimum))
  
}


#' Volatility Exponential Smoothing
#' @export 
#' @useDynLib amar vol_exp_smoothing
#' @rdname amar

vol.exp.smoothing <- function(x, lambda){
  #TODO verify parameters
  res <- .Call("vol_exp_smoothing", x, lambda)
  colnames(res) <- c("volatility", "residuals")
  res[,1] <- sqrt(res[,1])
  
  return(res)
}
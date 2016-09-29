
#' Fit autoregressive with OLS
#' Fit an autoregressive time series model to the data using the Ordinary Least Squares method to estimate the autoregressive coefficients.
#' @param x A numeric vector with data points.
#' @param p Order of model to fit. 
#' @param use.RcppEigen A logical variable. If \code{use.RcppEigen = TRUE}, the RcppEigen package is used to perform OLS estimation (recommended).
#' @param RcppEigen.method An integer specifining method used to solve OLS, see the RcppEigen package.
#' @param ... Other parameters that can be passed to \code{lm.fit} (If \code{use.RcppEigen = FALSE}).
#' @return An object of class "ar.fit", which contains the following fields:
#' \item{response}{The response vector.}
#' \item{p}{Order of the fitted model.}
#' \item{coefficient}{OLS estimates of the autoregressive coefficients.}
#' \item{fitted}{Fitted values.}
#' \item{residuals}{The difference between the response and the fitted value}
#' @export 
#' @rdname ar

ar.fit <- function(x, p, use.RcppEigen = TRUE,  RcppEigen.method = 1, ...){
  
  if(use.RcppEigen == TRUE) if(require("RcppEigen") == FALSE){
    warning("Package RcppEigen not available, using standard lm.fit")
    use.RcppEigen <- FALSE
  }
  
  x <- as.numeric(x)
  p <- as.integer(p)[1]
  n <- length(x)
  
  use.RcppEigen <- as.logical(use.RcppEigen)

  if(n-p < 1) stop("p has to be lower than length(x)")

  design <- ar.design.matrix(x, p)
  response <- ar.response(x, p)

  
  if(use.RcppEigen == TRUE) tmp <- fastLmPure(design, response, method = RcppEigen.method)
  else tmp <- lm.fit(design, response, ...)
  
  res <- list()
  class(res) <- "ar.fit"
  
  res$response <- response
  res$p <- p
  res$coefficients <- tmp$coefficients
  
  if(p>1) res$fitted <- design %*% res$coefficients
  else res$fitted <- design * res$coefficients

  res$residuals <- response - res$fitted
  
  
  return(res)
  
}

#' Predict \code{ar.fit} object
#' \code{predict} is a function to obtain predictions from  an \link{\code{ar.fit}} object.
#' @param object An object returned by \link{\code{ar.fit}}.
#' @param newx A design matrix with p columns to be used for prediction.
#' @return Vector with the predicted values.
#' @export 
#' @rdname ar

predict.ar.fit <- function(object, newx){
  
  if(missing(newx)) return(object$fitted)
  else{
    newx <- as.matrix(newx)
    if(ncol(newx) != object$p) stop("The numbe of columns in newx should match the order of the fitted model.")
    
    if(length(object$coefficients) > 1) return(newx%*%object$coefficients)
    else return(newx * object$coefficients)
  }
  
}


#' Train AR model.
#' @export 
#' @rdname ar

ar.train <- function(x, max.order, 
                      split.pct = c(0.5,0.25, 0.25),
                       remove.volatility = TRUE, 
                       target = c("hit.rate", "r.squared"), 
                       verbose = TRUE, ...){
  
  #TODO verify parameters
  
  res <- list()
  res$target <- match.arg(target, c("hit.rate", "r.squared"))
  
  #find the indices of the burn in, train, validate and test part

  n <- length(x)
  n.train <- floor(split.pct[1] * (n-max.order))
  n.validate <- floor(split.pct[2] * (n-max.order))
  n.test <- n - n.train - n.validate - max.order;
  
  burn.in.id <- 1:(max.order)
  train.id <- max.order + 1:n.train
  validate.id <- train.id[length(train.id)] + 1:n.validate
  test.id <- (validate.id[length(validate.id)]+1):n
  
  ##first remove the volatility optimising lambda for volatility.criterion
  if(remove.volatility == TRUE){
    
    tmp <- remove.volatility(x, c(validate.id), criterion = volatility.criterion)
    
    res$x <- tmp$residuals
    res$min.lambda <- tmp$min.lambda
    res$volatility <- tmp$volatility
    
  }else res$x <- x
  
  

  
  ## iterate through all potential solutions
  ## TODO add more verbosity
  
  if (res$target == "hit.rate") criterion.fun <- hit.rate 
  else if (res$target == "r.squared") criterion.fun <- r.squared
  else stop("target unknown")
  
  
  res$criterion.train <- res$criterion.validate <- rep(0, length(max.order))
  
  for (j in 1:max.order){
    
    if(verbose) cat("Checking ",j,"th out of ", max.order ," orders...\n")
    
    ar.fit.train <- ar.fit(res$x[c(train.id[1] - (j:1), train.id)], j)
    response.train <- ar.fit.train$response
    predicted.train <- predict(ar.fit.train)

    res$criterion.train[j] <- criterion.fun(response.train, predicted.train, k=j)
    
    design.validate <- ar.design.matrix(res$x[c(validate.id[1] - (j:1), validate.id)], j)
    response.validate <- ar.response(res$x[c(validate.id[1] - (j:1), validate.id)], j)
    predicted.validate <-  predict(ar.fit.train, design.validate)
    
    res$criterion.validate[j] <- criterion.fun(response.validate, predicted.validate, k=j)
    
    if(verbose) cat("Criterion on the train set:", res$criterion.train[j], "\n")
    if(verbose) cat("Criterion on the validate set:", res$criterion.validate[j], "\n")
    
  }
  
  ## final model selection
  
  
  if (res$target == "hit.rate" || res$target == "r.squared")  res$opt.order <- which.max(res$criterion.validate)
  else if(res$target == "aic") res$opt.order <- which.min(res$criterion.validate)
  

  ar.fit.train <- ar.fit(res$x[c(train.id[1] - (res$opt.order:1), train.id)], res$opt.order)

  design.train <- ar.design.matrix(res$x[c(train.id[1] - (res$opt.order:1), train.id)], res$opt.order)
  res$response.train <- ar.response(res$x[c(train.id[1] - (res$opt.order:1), train.id)], res$opt.order)
  res$predicted.train <-  predict(ar.fit.train, design.train)
  res$criterion.train <- criterion.fun(res$response.train, res$predicted.train, k=res$opt.order)
  res$train.id <- train.id
  
  design.validate <- ar.design.matrix(res$x[c(validate.id[1] - (res$opt.order:1), validate.id)], res$opt.order)
  res$response.validate <- ar.response(res$x[c(validate.id[1] - (res$opt.order:1), validate.id)], res$opt.order)
  res$predicted.validate <-  predict(ar.fit.train, design.validate)
  res$criterion.validate <- criterion.fun(res$response.validate, res$predicted.validate, k=res$opt.order)
  res$validate.id <- validate.id
  
  # ar.fit.train.validate <- ar.fit(res$x[c(train.id[1] - (res$opt.order:1), train.id)], res$opt.order)
  ar.fit.train.validate <- ar.fit(res$x[c(validate.id[1] - (res$opt.order:1), validate.id)], res$opt.order)
  
  design.test <- ar.design.matrix(res$x[c(test.id[1] - (res$opt.order:1), test.id)], res$opt.order)
  res$response.test <- ar.response(res$x[c(test.id[1] - (res$opt.order:1), test.id)], res$opt.order)
  res$predicted.test <-  predict(ar.fit.train.validate, design.test)
  res$criterion.test <- criterion.fun(res$response.test, res$predicted.test, k=res$opt.order)
  res$test.id <- test.id
  
  return(res)
}



#' blablabl
#' @export 
#' @useDynLib amar ar_design_matrix
#' @rdname ar

ar.design.matrix <- function(returns, order){
  #TODO validate the input
  return(.Call("ar_design_matrix", as.numeric(returns), as.integer(order)))
  
}


#' blablabl
#' @export 
#' @useDynLib amar ar_response_vector
#' @rdname ar

ar.response <- function(returns, order){
  
  returns <- as.numeric(returns)
  order <- as.integer(order)[1]
  

  #TODO validate the input
  
  return(.Call("ar_response_vector", as.numeric(returns), as.integer(order)))
  
}


#' @title Fitting AMAR models
#' @description R package for fitting amar models
#' @docType package
#' @name amar-package

NULL


#' Fit AMAR.
#' Fit AMAR time series model.
#' @param x A numeric vector with data points.
#' @param ... Currently not in use.
#' @return An object of class "amar.fit", which contains the following fields:
#' \item{response}{The response vector.}
#' \item{p}{Order of the fitted model.}
#' \item{coefficient}{OLS estimates of the AMAR coefficients.}
#' \item{fitted}{Fitted values.}
#' \item{residuals}{The difference between the response and the fitted value}
#' @export 
#' @rdname ar

amar.fit <- function(x, scales,...){
  
  #TODO check if parameters are valid
  
  design <- amar.design.matrix(x, scales)
  response <- amar.response(x, scales)

  tmp <- fastLmPure(design, response, 0)
  
  res <- list()
  class(res) <- "amar.fit"
  
  res$response <- response
  res$scales <- scales
  
  res$residuals <- tmp$residuals
  res$coefficients <- tmp$coefficients
  
  if(length(scales)>1) res$fitted <- design %*% res$coefficients
  else res$fitted <- design * res$coefficients
  # 
  # res$fitted <- tmp$fitted
  
  return(res)
  
}

#' Predict AMAR object
#' 
#' @export 
#' @rdname amar

predict.amar.fit <- function(object, newx){
  
  #TODO check if parameters are valid
  if(missing(newx)) return(object$fitted)
  else{
    
    return(newx%*%object$coefficients)
    
  }
}

#' Amar fit
#' @export 
#' @rdname amar

amar <- function(x,
                 max.order,
                 method = c("bic", "threshold"), 
                 threshold,
                 threshold.const = 0.5,
                 max.scales = 10, verbose=TRUE, use.RcppEigen = TRUE, ...){
  
  ##TODO verify parameters
  
  res <- list()
  res$method <- match.arg(method, c("bic", "threshold"))
  
  #fit a large AR model
  #TODO replace with a quicker method
  if(verbose) cat("Computing LSE coefficients for a large AR model.\n")
  res$ar.coefficients <- ar.fit(x, max.order, ...)$coefficients
  
  #detect scales using Narrowest-Over-Threshold Wild Binary Segmentation, get the entire solution path
  if(verbose) cat("Finding scales candidates using NOT.\n")
  w <- not(res$ar.coefficients, contrast = "pcwsConstMean",  ...)
  
  if(res$method=="threshold"){
    
    n <- length(x)

    
    if(missing(threshold)) threshold <- threshold.const *((log(n))^(3/2) / n^(1/2))
    if(length(threshold) != 1) stop("threshold should be a scalar")
    
    res$scales <- features(w, method="threshold", th = threshold)$cpt
    
    res$scales <- sort(unique(c(res$scales[!is.na(res$scales)], 1)))
    
    
  }else if(res$method=="bic"){
    
    scales.candidates <- list()
    n.scales.candidates <- 0
    
    for (j in 1:length(w$solution.path$cpt)) if (w$solution.path$n.cpt[j] <= max.scales){
      
      n.scales.candidates <- n.scales.candidates+1
      scales.candidates[[n.scales.candidates]]   <- sort(unique(c(1,w$solution.path$cpt[[j]])))
      
    }
    
    ## for each scale candidate, find the likelihood
    
    ic <- sapply(scales.candidates, function(scales){
      
      tmp <- amar.fit(x, scales)
      n <- length(tmp$residuals)
      
      return(n*log(mean(tmp$residuals^2))+ length(scales) * log(n))
      
    })
    
    res$ic <- ic
    res$scales <- scales.candidates[[which.min(ic)]]
    res$scales.candidates <- scales.candidates
    
  }
  
  tmp <- amar.fit(x, res$scales)
  
  res$response <- tmp$response
  res$fitted <- tmp$fitted
  res$residuals <- tmp$residuals
  res$coefficients <- tmp$coefficients
  
  class(res) <- "amar"
  return(res)
  
}

#' Simulate AMAR model
#' @export 
#' @rdname amar

amar.sim <- function(n, coefficients, scales, sigma=1){
  
  #TODO check the parameters
  scales <- unique(as.integer(scales))
  max.scale <- max(scales)
  #TODO check parameters values
  x <- ts(rnorm(n+max.scale, 0, sigma), start = 1 - max.scale)
  ar.coefficients <- rep(0, max.scale)
  
  for(j in 1:length(scales)) ar.coefficients[1:scales[j]] <- ar.coefficients[1:scales[j]] + coefficients[j]/scales[j]
  x <- filter(x, ar.coefficients, method = "recursive")
  
  return(x)
  
}

#' Train AMAR model
#' @export 
#' @rdname amar

amar.train <- function(x, 
                       max.order, 
                       split.pct = c(0.5, 0.25, 0.25),
                       remove.volatility = TRUE, 
                       volatility.criterion = c("jarque.bera"),
                       target = c("hit.rate", "r.squared"),
                       max.scales = 10, 
                       verbose = TRUE, 
                       ...){
  
  #TODO verify parameters
  res <- list()
  res$target <- match.arg(target, c("hit.rate", "r.squared"))
  
  #find the indices of the burn in, train, validate and test partf
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
    
    tmp <- remove.volatility(x, c( validate.id), criterion = volatility.criterion)
  
    res$x <- tmp$residuals
    res$min.lambda <- tmp$min.lambda
    res$volatility <- tmp$volatility
    
  }else res$x <- x
  
  if(verbose) cat("Computing LSE coefficients for a large AR model.\n")
  res$ar.coefficients.train <- ar.fit(res$x[c(burn.in.id, train.id)], max.order, ...)$coefficients
  
  #detect scales using Narrowest-Over-Threshold Wild Binary Segmentation, get the entire solution path
  if(verbose) cat("Finding scales candidates using NOT.\n")
  w <- not(res$ar.coefficients.train, 
           contrast = "pcwsConstMean", 
           method="not", 
           ...)

  scales.candidates <- list()
  # scales.candidates[[1]] <- c(1)

  
  for (j in 1:length(w$solution.path$cpt)) 
    if (w$solution.path$n.cpt[j] <= max.scales & w$solution.path$n.cpt[j]>0){

    scales.candidates <- c(scales.candidates,  list(sort(c(w$solution.path$cpt[[j]]))))
    scales.candidates <- c(scales.candidates,  list(sort(c(1, w$solution.path$cpt[[j]]))))
    
    }



  res$scales.candidates <- unique(scales.candidates)

  ## iterate through all potential solutions
  ## TODO add more verbosity
  
  if (res$target == "hit.rate") criterion.fun <- hit.rate 
  else if (res$target == "r.squared") criterion.fun <- r.squared
  
  
  res$criterion.train <- res$criterion.validate <- rep(0, length(res$scales.candidates))
  
  for (j in 1:length(res$scales.candidates)){
    
    if(verbose) cat("Checking ",j,"th scales...\n")
    
    scale.max <- max(res$scales.candidates[[j]])
    res$scales.candidates[[j]] <- sort(res$scales.candidates[[j]])
    
    amar.fit.train <- amar.fit(res$x[c(train.id[1] - (scale.max:1), train.id)], res$scales.candidates[[j]])
    response.train <- amar.fit.train$response
    predicted.train <- predict(amar.fit.train)
    res$criterion.train[j] <- criterion.fun(response.train, predicted.train)
    
    design.validate <- amar.design.matrix(res$x[c(validate.id[1] - (scale.max:1), validate.id)], res$scales.candidates[[j]])
    response.validate <- amar.response(res$x[c(validate.id[1] - (scale.max:1), validate.id)], res$scales.candidates[[j]])
    predicted.validate <-  predict(amar.fit.train, design.validate)

    res$criterion.validate[j] <- criterion.fun(response.validate, predicted.validate)
    
    if(verbose) cat("Criterion on the train set:", res$criterion.train[j], "\n")
    if(verbose) cat("Criterion on the validate set:", res$criterion.validate[j], "\n")
    
  }
  
  ## final model selection
  
  res$scales <- res$scales.candidates[[which.max(res$criterion.validate)]]

  scale.max <- max(res$scales)

  # amar.fit.train <- amar.fit(res$x[c(train.id[1] - (scale.max:1), train.id, validate.id)], res$scales)
  amar.fit.train <- amar.fit(res$x[c(train.id[1] - (scale.max:1), train.id)], res$scales)
  
  design.train <- amar.design.matrix(res$x[c(train.id[1] - (scale.max:1), train.id)], res$scales)
  res$response.train <- amar.response(res$x[c(train.id[1] - (scale.max:1), train.id)], res$scales)
  res$predicted.train <-  predict(amar.fit.train, design.train)
  res$criterion.train <- criterion.fun(res$response.train, res$predicted.train)
  res$train.id <- train.id
  

  design.validate <- amar.design.matrix(res$x[c(validate.id[1] - (scale.max:1), validate.id)], res$scales)
  res$response.validate <- amar.response(res$x[c(validate.id[1] - (scale.max:1), validate.id)], res$scales)
  res$predicted.validate <-  predict(amar.fit.train, design.validate)
  res$criterion.validate <- criterion.fun(res$response.validate, res$predicted.validate)
  res$validate.id <- validate.id

  amar.fit.train.validate <- amar.fit(res$x[c(validate.id[1] - (scale.max:1), validate.id)], res$scales)
  # amar.fit.train.validate <- amar.fit(res$x[c(train.id[1] - (scale.max:1),  train.id)], res$scales)
  
  design.test <- amar.design.matrix(res$x[c(test.id[1] - (scale.max:1), test.id)], res$scales)
  res$response.test <- amar.response(res$x[c(test.id[1] - (scale.max:1), test.id)], res$scales)
  res$predicted.test <-  predict(amar.fit.train.validate, design.test)
  # res$predicted.test <-  predict(amar.fit.validate, design.test)
  res$criterion.test <- criterion.fun(res$response.test, res$predicted.test)
  res$test.id <- test.id
  
  return(res)
}

#' Create AMAR design matrix
#' @export 
#' @useDynLib amar amar_design_matrix
#' @rdname amar

amar.design.matrix <- function(x, scales){
  
  #TODO check parameters
  
  design <- .Call("amar_design_matrix", as.numeric(x), as.integer(sort(setdiff(scales,0))))
  
  if(0 %in% scales) design <- cbind(rep(1, nrow(design)), design)
  
  
  return(design)
  
}


#' Create AMAR response
#' @export 
#' @useDynLib amar amar_response_vector
#' @rdname amar

amar.response <- function(returns, scales){
  
  #TODO check parameters
  
  return(.Call("amar_response_vector", as.numeric(returns),as.integer(sort(scales))))
  
  
}


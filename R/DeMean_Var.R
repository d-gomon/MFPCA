#' @title De-Mean data according to specified variable
#' 
#' @description Given realizations from a square integrable process, make the mean
#' over the realizations zero adjusted for a certain covariate.
#' 
#' @details Currently only numeric variables are supported
#' 
#' @param mFData A multiFunData object 
#' 
#' @importFrom funData nObs
#' 
#' @export
#' @noRd
#' @keywords internal
#' 
#' 


DeMean_Var <- function(mFData, var = NULL, nbasis = NULL, extraInfo = TRUE){
  #Possibly usefull for future?
  if(is.null(var)){
    return(mFData)
  }
  
  if(!all(is.numeric(var), length(var) == nObs(mFData))){
    stop("Var must be passed as a numeric vector of length nObs(mFData).")
  }
  if(!all(inherits(mFData, "multiFunData"))){
    stop("mFData must be of class multiFunData")
  }
  
  
  
  p <- length(mFData)
  
  if(isTRUE(extraInfo)){
    meanFuns <- vector(mode = "list", length = p)  
  }
  
  X_final <- vector(mode = "list", length = p)
  for(i in 1:p){
    temp_deMeaned <- DeMean_Var_FunData(funData = mFData[[i]], var = var, nbasis = nbasis, extraInfo = extraInfo)
    X_final[[i]] <- temp_deMeaned$funData
    if(isTRUE(extraInfo)){
      meanFuns[[i]]$meanFun = temp_deMeaned$meanFun
      meanFuns[[i]]$gam_mod = temp_deMeaned$gam_mod
    }
  }
  
  out <- list(mFData = multiFunData(X_final))
  
  if(isTRUE(extraInfo)){
    out$meanFuns = meanFuns
  }
  
  return(out)
}


#' @export
#' @noRd
#' @keywords internal

DeMean_Var_FunData <- function(funData, var = NULL, nbasis = NULL, extraInfo = TRUE){
  if(!all(is.numeric(var), length(var) == nObs(funData))){
    stop("Var must be passed as a numeric vector of length nObs(mFData).")
  }
  
  if(!all(inherits(funData, "funData"))){
    stop("funData must be of class funData")
  }
  
  if(is.null(nbasis)){
    nbasis = 10
  }
  
  Y <- funData@X
  X <- funData@argvals[[1]]
  
  D = NCOL(Y)
  if(D != length(X)){ # check if number of observation points in X & Y are identical
    stop("different number of (potential) observation points differs in X and Y!")
  }
  I = NROW(Y)
  #CREATE NEW d.vec_adj, which represents var at baseline. Basically, substract var at baseline for each subject from X
  #var.bl must be vector containing var in same time unit as funData object!!!
  #We create an adjusted vector, where each patient has individual starting times.
  d.vec_adj <- rep(X, each = I) + var
  gam0 = mgcv::gam(as.vector(Y) ~ s(d.vec_adj, k = nbasis))
  mu = mgcv::predict.gam(gam0, newdata = data.frame(d.vec_adj = sort(unique(d.vec_adj))))
  mu_mat <- matrix(mu[match(d.vec_adj, sort(unique(d.vec_adj)))], I, D, byrow = FALSE)
  Y.tilde = Y - mu_mat
  
  out <- list(funData = funData(argvals = X, X = Y.tilde))
  
  if(isTRUE(extraInfo)){
    out$meanFun = funData(argvals = sort(unique(d.vec_adj)), X = matrix(mu, nrow = 1))
    out$gam_mod = gam0
  }
  return(out)
}


#' @export
#' @keywords internal
#' @noRd

DeMean_test <- function(mFData, var = NULL, meanFuns){
  if(is.null(var)){
    return(mFData)
  }
  
  if(!all(is.numeric(var), length(var) == nObs(mFData))){
    stop("Var must be passed as a numeric vector of length nObs(mFData).")
  }
  if(!all(inherits(mFData, "multiFunData"))){
    stop("mFData must be of class multiFunData")
  }
  
  p <- length(mFData)
  X_final <- vector(mode = "list", length = p)
  
  if(length(meanFuns) != p){
    stop("Age de-meaning model must be fitted on training data with same amount of variables as test data.")
  }
  for(i in 1:p){
    temp_deMeaned <- DeMean_test_funData(funData = mFData[[i]], var = var, meanFun = meanFuns[[i]])
    X_final[[i]] <- temp_deMeaned$funData
  }
  out <- list(mFData = multiFunData(X_final))
  return(out)
}


#' @export
#' @keywords internal
#' @noRd
DeMean_test_funData <- function(funData, var = NULL, meanFun){
  Y.pred <- funData@X
  X <- funData@argvals[[1]]
  
  
  D = NCOL(Y.pred)
  if(D != length(X)){ # check if number of observation points in X & Y are identical
    stop("different number of (potential) observation points differs in X and Y!")
  }
  I.pred = NROW(Y.pred)
  
  d.vec_adj_pred <- rep(X, each = I.pred) + var
  #predict values at all distinct ages
  mu.pred = mgcv::predict.gam(meanFun$gam_mod, newdata = data.frame(d.vec_adj = sort(unique(d.vec_adj_pred))))
  #create matrix with rows = subjects and columns = number of observations
  #match subject means with correct entry in prediction
  mu.pred_mat <- matrix(mu.pred[match(d.vec_adj_pred, sort(unique(d.vec_adj_pred)))], I.pred, D, byrow = FALSE)
  Y.tilde = Y.pred - mu.pred_mat
  #mu_mat <- mu.pred_mat
  return(list(funData = funData(argvals = X, X = Y.tilde)))
}

  
  
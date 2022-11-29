#' @title Landmark data (step 1)
#' 
#' @description Landmark data (remove information on patients with failures before landmark time).
#' 
#' @details mFData_train, X_baseline_train, Y_surv_train and age_train must have the same row order. Each row should represent the same patient.
#' The same is true for prediction data.
#' 
#' @param time Landmark time. Default = NULL (no landmarking).
#' @param mFData_train Training data in multiFunData format.
#' @param X_baseline_train Baseline covariates for training data
#' @param Y_surv_train Survival information (time, event) for training data
#' @param age_train (optional) Numeric vector of ages at baseline for training data
#' @param mFData_pred Data for which predictions are required in multiFunData format.
#' @param X_baseline_pred Baseline covariates for testing data
#' @param Y_surv_pred Survival information (time, event) for prediction data. Can be left empty,
#' in that case predictions will be made for the whole test data.
#' @param age_pred (optional) Numeric vector of ages at baseline for prediction data.
#' 
#' 
#' 
#' @export
#' 
#' @keywords internal
#' 

landmark_data <- function(time = NULL, mFData_train, X_baseline_train, Y_surv_train, age_train = NULL, mFData_pred = NULL, 
                          X_baseline_pred = NULL, Y_surv_pred = NULL, age_pred = NULL){
  
  
  
  n_train <- nrow(Y_surv_train)
  if(!is.null(Y_surv_pred)){
    n_pred <- nrow(Y_surv_pred)  
  }
  
  if(!is.null(time)){
    #
    #Landmark training data
    #
    #Determine which argvals to use for landmarking by cutting of at landmark_time
    max_id <- max(which(mFData_train[[1]]@argvals[[1]] <= time))
    #Determine which observations to use by eliminating people with failures before landmark_time (requires Y_train)
    which_alive_train <- which(Y_surv_train[, "time"] > time)
    X_baseline_train <- X_baseline_train[which_alive_train, , drop = FALSE]
    Y_surv_train <- Y_surv_train[which_alive_train,]
    #Subset the data according to above 2 conditions.
    mFData_train <- subset(mFData_train, obs = which_alive_train, argvals = list(mFData_train[[1]]@argvals[[1]][1:max_id]))
    message(paste("Retained ",length(which_alive_train)*100/n_train, "percent of original training observations."))
    
    #
    #Landmark prediction data
    #
    #Determine which argvals to use for landmarking by cutting of at landmark_time
    if(!is.null(mFData_pred)){
      max_id_pred <- max(which(mFData_pred[[1]]@argvals[[1]] <= time))
      #original num of observations
      nrow_orig_pred <- nrow(mFData_pred[[1]]@X)
      if(!is.null(Y_surv_pred)){
        #Determine which observations to use by eliminating people with failures before landmark_time (requires Y_train)
        which_alive_pred <- which(Y_surv_pred[, "time"] > time)
        X_baseline_pred <- X_baseline_pred[which_alive_pred, , drop = FALSE]
        Y_surv_pred <- Y_surv_pred[which_alive_pred,]
        mFData_pred <- subset(mFData_pred, obs = which_alive_pred, argvals = list(mFData_pred[[1]]@argvals[[1]][1:max_id_pred]))
      } else{
        mFData_pred <- subset(mFData_pred, argvals = list(mFData_pred[[1]]@argvals[[1]][1:max_id_pred]))
      }
      message(paste("Retained ",length(which_alive_pred)*100/n_pred, "percent of original prediction observations."))
    }
  } else{
    message("No landmarking performed. Retained 100 percent of original observations.")
  }
  
  #
  #When different de-meaning required.
  #
  #If we landmarked, and want to de-mean according to age, we should adjust the age.bl vars
  #to remove baseline covs from people removed in landmarking.
  if(!is.null(time)){
    if(!is.null(age_train)){
      age_train <- age_train[which_alive_train]
    } 
    if(!is.null(age_pred)){
      age_pred <- age_pred[which_alive_pred]
    }  
  }
  
  return(list(time = time, 
              mFData_train = mFData_train,
              X_baseline_train = X_baseline_train, 
              Y_surv_train = Y_surv_train, 
              age_train = age_train, 
              mFData_pred = mFData_pred, 
              X_baseline_pred = X_baseline_pred, 
              Y_surv_pred = Y_surv_pred, 
              age_pred = age_pred))
}



#' 
#' @description Landmarks training data by making all observations after landmark time equal to NA.
#' 
#' @details Only difference with landmark_data is that instead of removing observations after landmark time,
#' these observations are here set to NA so that when "fake" landmarking the scores can be determined for test data.
#' Only apply this to test data.
#' 
#' @export
#' 
#' @keywords internal
#' 

landmark_data_fake <- function(time = NULL, mFData, X_baseline, Y_surv, age = NULL){
  n <- nrow(Y_surv)
  n_col <- ncol(mFData[[1]]@X)
  
  if(!is.null(time)){
    #
    #Landmark data
    #
    #Determine which argvals to use for landmarking by cutting of at landmark_time
    max_id <- max(which(mFData[[1]]@argvals[[1]] <= time))
    #Determine which observations to use by eliminating people with failures before landmark_time (requires Y_train)
    which_alive <- which(Y_surv[, "time"] > time)
    X_baseline <- X_baseline[which_alive, , drop = FALSE]
    Y_surv <- Y_surv[which_alive,]
    #Subset the data according to above 2 conditions.
    mFData <- subset(mFData, obs = which_alive)
    message(paste("Retained ",length(which_alive)*100/n, "percent of original testing observations."))
    if(max_id < n_col){
      for(i in 1:length(mFData)){
        mFData[[i]]@X[, (max_id +1):n_col] <- NA
      }
      message("Set all observations after landmark time to NA for test data.")
    }
    
    
    
    if(!is.null(age)){
      age <- age[which_alive]
    } 
  }
  
  return(list(time = time, 
              mFData = mFData,
              X_baseline = X_baseline, 
              Y_surv = Y_surv, 
              age = age))
}
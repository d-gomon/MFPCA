#' @title Fit multivariate Cox model (outdated)
#' 
#' @param mscores output of predict_mscores function (contains scores for train and test data).
#' @param Y_train Survival information (time, event) for training data
#' @param X_baseline_train Baseline covariates data frame for training data
#' @param landmark_time Landmark time for predictions. See equation in Section 2.4 in Li/Luo
#' 
#' 
#' 
#' 
#' 
#' @import glmnet
#' 
#' 
#' 
#' 

fit_mcox <- function(mscores, Y_train, X_baseline_train = NULL, X_baseline_test = NULL,
                     Y_test = NULL, landmark_time = NULL){
  #We extract the scores from the input. Both for train and test data.
  if(!is.null(mscores)){
    scores_train = mscores$mscores_train
    scores_test = mscores$mscores_test
  }
  
  #Determine some variables which can be extracted from the input.
  p_baseline <- ncol(X_baseline_train)
  p_long <- ncol(scores_train)
  
  #Construct the data for train and test data (baseline & mscores)
  X_train <- cbind(X_baseline_train, scores_train)
  if(!is.null(X_baseline_test)){
    X_test <- cbind(X_baseline_test, scores_test)  
  }
  
  #Use CV in glmnet to determine optimal lambda parameter on train data
  #and fit the regularized cox model on training data
  cv_fit_train <- cv.glmnet(x = X_train, y = Y_train,
                            family = "cox", type.measure = "C",
                            penalty.factor = c(rep(0, p_baseline), rep(1, p_long))) 
  fit_surv_train <- glmnet(X_train, Y_train,
                           penalty.factor = c(rep(0, p_baseline), rep(1, p_long)),
                           family = "cox", lambda = cv_fit_train$lambda.min)
  
  
  
  
  return(list(fitted_cox = fit_surv_train,
              Y_train = Y_train,
              X_train = X_train))
  
  #Now use this model from training data to predict for test data.
  
  #final <- predict_surv(model = fit_surv_train, x = X_surv, y = Y_surv, 
  #                      newx = X_surv[1:100,], times.pred = c(3, 4, 5, 6, 7, 8), 
  #                      landmark_time = 1)
}






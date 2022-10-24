#' @title Predict survival probabilities for new patients (step 3)
#' 
#' @description Predict survival for prediction data
#' 
#' 
#' @param step2 output of get_mscores
#' @param times_pred Times at which prediction is required. Default = landmark_time (useless).
#' @param accuracy_train Should accuracy measures (AUC/Brier score) be calculated on the training set?
#' @param reg_baseline Should baseline covariates be regularized?
#' @param reg_long Should longitudinal covariates be regularized?
#' 
#' @importFrom pec cindex
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @import prodlim
#' @import riskRegression
#' 
#' 


predict_surv <- function(step2, times_pred = step2$landmark_time, 
                         accuracy_train = FALSE, reg_baseline = FALSE, reg_long = TRUE){
  #step2 contains the variables:
  #type = c("scores", "AUC", "pp", "uscores"), M = NULL, 
  #uniExpansions = NULL, verbose = FALSE
  
  #step2 also contains step1, which includes the variables:
  #time = NULL, mFData_train, X_baseline_train, Y_surv_train, age_train = NULL, 
  #mFData_pred = NULL, X_baseline_pred = NULL, Y_surv_pred = NULL, age_pred = NULL
  
  if(!missing(step2)){
    list2env(step2, envir = environment())
  }


  
  p_baseline <- ncol(step1$X_baseline_train)
  p_long <- ncol(X_train) - p_baseline
  if(isTRUE(reg_baseline)){
    reg_base <- rep(1, p_baseline)
  } else{
    reg_base <- rep(0, p_baseline)
  }
  if(isTRUE(reg_long)){
    reg_lon <- rep(1, p_long)
  } else{
    reg_lon <- rep(0, p_long)
  }
  


  #First we fit a cox model on the training data. 
  #Use CV in glmnet to determine optimal lambda parameter on train data
  #and fit the regularized cox model on training data
  if(isFALSE(reg_baseline) & isFALSE(reg_long)){
    #cv_fit_train <- list(lambda.min = 0)
    reg_lon <- rep(1, p_long)
  } 
  cv_fit_train <- cv.glmnet(x = X_train, y = Y_train,
                            family = "cox", type.measure = "C",
                            penalty.factor = c(reg_lon, reg_base))
  if(isFALSE(reg_baseline) & isFALSE(reg_long)){
    s_opt <- 0
  } else{
    s_opt <- cv_fit_train$lambda.min
  }
  
  
  lp_train <- predict(cv_fit_train, newx = X_train, type = "link", s = s_opt)
  lp_pred <- predict(cv_fit_train, newx = X_pred, type = "link", s = s_opt)
  
  if(is.null(landmark_time)){
    landmark_time <- argvals(step1$mFData_pred)[[1]][[1]][1]
  }

  #Use the hazard to output predicted survival probabilities. 
  #But use the linear predictor to evaluate AUC - less smoothing so more accurate.
  haz <- hdnom::glmnet_basesurv(time = Y_train[, "time"], event = Y_train[, "status"], 
                                lp = lp_train,
                                times.eval = c(landmark_time, times_pred), 
                                centered = FALSE)
  #Determine baseline survival probability at landmark time:
  base_land_surv <- exp(-haz$cumulative_base_hazard[which(haz$times == landmark_time)])
  #Determine baseline survival probability at each time
  surv_fit <- exp(-haz$cumulative_base_hazard)
  #Baseline survival ratio S_0(time)/S_0(landmark_time)
  surv_ratio <- surv_fit/base_land_surv
  #Now exponentiate each term with exp(lp) to obtain prediction
  prob_surv_pred <- matrix(surv_ratio^(rep(exp(lp_pred), each = length(surv_ratio))), 
                      ncol = length(surv_ratio), byrow = TRUE)
  prob_surv_train <- matrix(surv_ratio^(rep(exp(lp_train), each = length(surv_ratio))), 
                            ncol = length(surv_ratio), byrow = TRUE)
  colnames(prob_surv_pred) <- c(landmark_time, times_pred)
  rownames(prob_surv_pred) <- rownames(X_pred)
  colnames(prob_surv_train) <- c(landmark_time, times_pred)
  rownames(prob_surv_train) <- rownames(X_train)
  
  #Evaluate AUC for training data.
  if(isTRUE(accuracy_train)){
    AUC_train <- sapply(times_pred, function(t) survivalROC::survivalROC(
      Stime = Y_train[, "time"], # Event time or censoring time for subjects
      status = Y_train[, "status"], # Indicator of status, 1 if death or event, 0 otherwise
      marker = lp_train, # Predictor or marker value
      entry = NULL, # Entry time for the subjects, default is NULL
      predict.time = t, # Time point of the ROC curve
      cut.values = NULL, # marker values to use as a cut-off for calculation of sensitivity and specificity
      method = "NNE", 
      span = 0.25 * nrow(Y_train)^(-0.2) # small span yield moderate smoothing, how to select?
    )$AUC)
  } else{
    AUC_train = NULL
    Brier_train = NULL
  }
  
  
  #Evaluate AUC here for prediction data if Y_pred specified. Otherwise, don't. Just don't.
  if(!is.null(Y_pred)){
    
    traindat <- cbind(step2$Y_train, step2$traindat)
    preddat <- cbind(step2$Y_pred, step2$preddat)
    cox_train_mod <- coxph(Surv(time, status) ~ ., data = traindat, x = TRUE, y = TRUE)
    
    #https://www.jesseislam.com/post/brier-score/
    featureNames <- paste(colnames(step1$X_baseline_pred), collapse = " + ")
    Brier_temp <- suppressMessages(riskRegression::Score(object = list("Cox1" = cox_train_mod),
                                        #formula = cox_train_mod$formula,
                          formula = as.formula(paste0("Surv(time, status) ~", featureNames)),
                          data = preddat,
                          exact = FALSE, # Do not predict at event times
                          times = times_pred,
                          conf.int = FALSE,
                          cens.model = "cox", # Method for estimating inverse probability of censoring weights:
                          splitMethod = "none",
                          B = 0,
                          verbose = TRUE))
    Brier_pred <- subset(Brier_temp$Brier$score, model == "Cox1")$Brier
    names(Brier_pred) <- subset(Brier_temp$Brier$score, model == "Cox1")$times
    
    AUC_pred <- Brier_temp$AUC$score$AUC
    names(AUC_pred) <- Brier_temp$AUC$score$times
  }
  
  out <- list(prob_surv_pred = prob_surv_pred,
              prob_surv_train = prob_surv_train,
              cox_train = cv_fit_train,
              times_pred = times_pred,
              AUC_pred = AUC_pred,
              AUC_train = AUC_train,
              Brier_pred = Brier_pred,
              Brier_train = Brier_train,
              landmark_time = landmark_time,
              lp_pred = lp_pred,
              M = M,
              PVE = mPVE)
  return(out)
  
}



# Decided to no longer calculate C-index.
# if(isTRUE(C_pred)){
#   #C_pred_glmnet does not use IPCW, but you can make it use the IPCW by manually specifying weights in Cindex
#   C_pred_glmnet <- sapply(times_pred,
#                           function(x) glmnet::Cindex(lp_pred, y = Surv(Y_pred[, "time"], ifelse(Y_pred[, "time"] <= x & Y_pred[, "status"] == 1, 1, 0))))
#   names(C_pred_glmnet) <- times_pred
#   
#   #C_pred uses IPCW, but does not use shrunk coefficients. 
#   C_pred <- suppressWarnings(pec::cindex(list("Cox1" = cox_train_mod),
#                                          formula = Surv(time, status) ~ .,
#                                          data = preddat,
#                                          eval.times = times_pred)$AppCindex$Cox1)
#   names(C_pred) <- times_pred
# }

# Used to calculate AUC with survivalROC
# #Check if there are any survivals past the landmark time
# event.times <- Y_pred[Y_pred[, "status"] == 1, "time"]
# if(isTRUE(AUC_pred)){
#   AUC_pred <- sapply(times_pred, function(t){
#     if (all(!(event.times <= t))) {
#       mess <- paste(
#         "No event (surv.new$event == 1) is observed between landmark time", landmark_time, 
#         "and prediction time", t)
#       warning(mess)
#       # Store no result
#       return(NA)
#     } else{
#       survivalROC::survivalROC(
#         Stime = Y_pred[, "time"], # Event time or censoring time for subjects
#         status = Y_pred[, "status"], # Indicator of status, 1 if death or event, 0 otherwise
#         marker = lp_pred, # Predictor or marker value
#         entry = NULL, # Entry time for the subjects, default is NULL
#         predict.time = t, # Time point of the ROC curve
#         cut.values = NULL, # marker values to use as a cut-off for calculation of sensitivity and specificity
#         method = "NNE", 
#         span = 0.25 * nrow(Y_pred)^(-0.2) # small span yield moderate smoothing, how to select?
#       )$AUC
#     }
#   })
#   names(AUC_pred) <- times_pred
# } else{
#   AUC_pred <- NULL
# }






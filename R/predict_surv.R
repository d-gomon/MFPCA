#' @title Predict survival probabilities for new patients (step 3)
#' 
#' @description Predict survival for prediction data
#' 
#' 
#' @param step2 output of get_mscores
#' @param times_pred Times at which prediction is required. Default = landmark_time (useless).
#' @param AUC_train Should AUC for training data be determined?
#' 
#' 


predict_surv <- function(step2, times_pred = step2$landmark_time, 
                         AUC_train = FALSE, reg_baseline = FALSE, reg_long = TRUE){
  
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
  if(isTRUE(AUC_train)){
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
  }
  
  
  #Evaluate AUC here for prediction data if Y_pred specified. Otherwise, don't. Just don't.
  if(!is.null(Y_pred)){
    #Check if there are any survivals past the landmark time
    event.times <- Y_pred[Y_pred[, "status"] == 1, "time"]
    
    AUC_pred <- sapply(times_pred, function(t){
      if (all(!(event.times <= t))) {
        mess <- paste(
          "No event (surv.new$event == 1) is observed between landmark time", landmark_time, 
          "and prediction time", t)
        warning(mess)
        # Store no result
        return(NA)
      } else{
        survivalROC::survivalROC(
          Stime = Y_pred[, "time"], # Event time or censoring time for subjects
          status = Y_pred[, "status"], # Indicator of status, 1 if death or event, 0 otherwise
          marker = lp_pred, # Predictor or marker value
          entry = NULL, # Entry time for the subjects, default is NULL
          predict.time = t, # Time point of the ROC curve
          cut.values = NULL, # marker values to use as a cut-off for calculation of sensitivity and specificity
          method = "NNE", 
          span = 0.25 * nrow(Y_pred)^(-0.2) # small span yield moderate smoothing, how to select?
        )$AUC
      }
    })

    # brier_pred <- tryCatch({
    #   res.bs <- pec::pec(
    #     # A matrix with predicted probabilities, dimension of n subjects by m times 
    #     object = fit_surv_train,
    #     #          formula = Surv(time, event) ~ AGE,
    #     formula = Surv(time, status) ~ AGE + PTGENDER + PTEDUCAT + status.bl + APOE4,
    #     data = as.data.frame(cbind(Y_train,step1$X_baseline_train)),
    #     exact = FALSE, # Do not predict at event times
    #     times = event.times, 
    #     #times = 0:15, 
    #     cens.model = "cox", # Method for estimating inverse probability of censoring weights:
    #     splitMethod = "none",
    #     B = 0,
    #     verbose = TRUE
    #   )$AppErr$model[-1]
    # }, error = function(e) {
    #   message(e)
    #   return(NA)
    # }, finally = {
    # })
  }
  return(list(prob_surv_pred = prob_surv_pred,
              prob_surv_train = prob_surv_train,
              cox_train = cv_fit_train,
              times_pred = times_pred,
              AUC_pred = AUC_pred,
              AUC_train = AUC_train,
              #Brier_pred = brier_pred,
              landmark_time = landmark_time,
              lp_pred = lp_pred))
}






predict_surv_new <- function(fit_mcox = NULL, X_surv_test,
                             Y_surv_test = NULL,
                             times.pred, landmark_time){
  #Extract (cumulative) baseline hazard
  x <- fit_mcox$X_surv_train
  y <- fit_mcox$Y_surv_train
  newx <- X_surv_test
  model <- fit_mcox$fit_surv_train
  haz <- hdnom::glmnet_basesurv(time = y[,1], event = y[,2], 
                                lp = predict(object = model, newx = x),
                                times.eval = c(landmark_time, times.pred), 
                                centered = FALSE)
  #Determine baseline survival probability at landmark time:
  base_land_surv <- exp(-haz$cumulative_base_hazard[which(haz$times == landmark_time)])
  #Determine baseline survival probability at each time
  surv_fit <- exp(-haz$cumulative_base_hazard)
  #Baseline survival ratio S_0(time)/S_0(landmark_time)
  surv_ratio <- surv_fit/base_land_surv
  #Calculate exp(lp) for new observations
  lp <- predict(model, newx = newx, type = "link")
  #Now exponentiate each term with exp(lp) to obtain prediction
  prob_surv <- matrix(surv_ratio^(rep(exp(lp), each = length(surv_ratio))), 
                      ncol = length(surv_ratio), byrow = TRUE)
  colnames(prob_surv) <- c(landmark_time, times.pred)
  rownames(prob_surv) <- rownames(newx)
  
  if(!is.null(Y_surv_test)){
    temp_res <- sapply(times.pred, 
                       FUN = function(t) 
                         survivalROC::survivalROC(Stime = Y_surv_test$time,
                                                  status = ifelse(Y_surv_test$status == "censored", 0, 1), 
                                                  marker = lp,
                                                  entry = rep(0, nrow(Y_surv_test)), predict.time = t,
                                                  method = "NNE", span = 0.25*nrow(Y_surv_test)^(-0.2))$AUC)
  }
  
  return(list(prob_surv = prob_surv,
              AUC = temp_res))
}


# predict_surv <- function(model, x, y, newx, times.pred, landmark_time){
#   #Extract (cumulative) baseline hazard
#   haz <- hdnom::glmnet_basesurv(time = y[,1], event = y[,2], 
#                   lp = predict(object = model, newx = x),
#                   times.eval = c(landmark_time, times.pred), centered = FALSE)
#   #Determine baseline survival probability at landmark time:
#   base_land_surv <- exp(-haz$cumulative_base_hazard[which(haz$times == landmark_time)])
#   #Determine baseline survival probability at each time
#   surv_fit <- exp(-haz$cumulative_base_hazard)
#   #Baseline survival ratio S_0(time)/S_0(landmark_time)
#   surv_ratio <- surv_fit/base_land_surv
#   #Calculate exp(lp) for new observations
#   exp_lp <- predict(model, newx = newx, type = "response")
#   #Now exponentiate each term with exp(lp) to obtain prediction
#   prob_surv <- matrix(surv_ratio^(rep(exp_lp, each = length(surv_ratio))), 
#                       ncol = length(surv_ratio), byrow = TRUE)
#   colnames(prob_surv) <- c(landmark_time, times.pred)
#   rownames(prob_surv) <- rownames(newx)
#   return(prob_surv)
# }

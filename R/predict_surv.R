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
#' @param alpha Regularization scaling. Alpha = 1 (lasso), alpha = 0 (ridge), between 0 and 1 = elastic net. Default = 1 (lasso)
#' @param IPCW_vars Character vector indicating which variables to use for IPCW.
#' Extra options are "all" (all baseline variables) and "none" (don't perform IPCW). Default is no IPCW ("none").
#' 
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @importFrom survival coxph.control
#' @import prodlim
#' @import riskRegression
#' @import glmnet
#' @importFrom success extract_hazard
#' 
#' 
#' @export
#' @keywords internal
#' 
#' 


predict_surv <- function(step2, times_pred, 
                         reg_baseline = FALSE, 
                         reg_long = TRUE, alpha = 1, accuracy_train = NULL, IPCW_vars = c("none", "all")){
  
  if(identical(IPCW_vars, c("none", "all"))){
    IPCW_vars = "none"
  }
  
  if(min(times_pred < step2$landmark_time)){
    stop("Prediction times should be larger than landmark time!")
  }
  
  #Read step2 into active environment
  if(!missing(step2)){
    list2env(step2, envir = environment())
  }
  #step2 contains the variables:
  #type = c("scores", "AUC", "pp", "uscores"), M = NULL, 
  #uniExpansions = NULL, verbose = FALSE
  #step2 also contains step1, which includes the variables:
  #time, mFData_train, X_baseline_train, Y_surv_train, age_train, 
  #mFData_pred, X_baseline_pred, Y_surv_pred, age_pred

  #Create vector to feed to cv.glmnet determining which variables to regularize
  if(isTRUE(reg_baseline)){
    reg_base <- rep(1, p_base)
  } else{
    reg_base <- rep(0, p_base)
  }
  if(isTRUE(reg_long)){
    reg_lon <- rep(1, p_long)
  } else{
    reg_lon <- rep(0, p_long)
  }
  if(isFALSE(reg_baseline) & isFALSE(reg_long)){
    reg_base <- rep(0, p_base)
    reg_lon <- rep(0, p_long)
  } 
  
  if(isFALSE(reg_baseline) & isFALSE(reg_long)){
    #If we don't want to perform regularization we cannot use glmnet to fit the model
    #Setting lambda to 0 for some reason doesn't work in glmnet 
    featureNames_cox <- paste(colnames(step2$traindat), collapse = " + ")
    survival_from_glm <- coxph(as.formula(paste0("Surv(time, status) ~", featureNames_cox)), data = cbind(step2$Y_train, step2$traindat), x = TRUE, y = TRUE)
    #We predict the linear predictor using the fitted cox model for both training and test data.
    lp_train <- predict(survival_from_glm, type = "lp")
    lp_pred <- predict(survival_from_glm, newdata = step2$preddat, type = "lp")
    cv_fit_train <- survival_from_glm
    
  }else{
    #Then we fit a cox model on the training data. 
    #Use CV in glmnet to determine optimal lambda parameter on train data
    #and fit the regularized cox model on training data
    cv_fit_train <- cv.glmnet(x = X_train, y = Y_train,
                              family = "cox", type.measure = "deviance", alpha = alpha,
                              penalty.factor = c(reg_lon, reg_base))
    s_opt <- "lambda.min"
    
    #We predict the linear predictor using the fitted cox model for both training and test data.
    lp_train <- predict(cv_fit_train, newx = step2$X_train, type = "link", s = s_opt)
    lp_pred <- predict(cv_fit_train, newx = step2$X_pred, type = "link", s = s_opt)
    
    
    #Transfer glm model to survival model for use in riskRegression::Score()
    #We do this because coxph has a predictRisk method. ?predictRisk.coxph
    #Giving the (regularized) linear predictor to coxph() allows to use basehaz()
    #to extract the (cumulative) baseline hazard from fitted regularized model
    df.orig <- as.data.frame(cbind(Y_train, lp_train))
    colnames(df.orig)[3] <- "linpred"
    survival_from_glm <- coxph(Surv(time = time, event = status) ~ linpred,         
                               data = df.orig, init = 1, 
                               control = coxph.control(iter.max = 0), x = TRUE, 
                               y = TRUE)
  }

  
  
  
  
  
  #DOESNT REALLY DO ANYTHING YET. FUTURE?
  if(is.null(landmark_time)){
    landmark_time <- argvals(step1$mFData_pred)[[1]][[1]][1]
  }
  
  
  #We use the following formula to predict survival:
  #S(s'|s) = (S_0(s')/S_0(s))^(exp(linpred))
  #See Li/Luo (2019) Dynamic prediction of Alzheimer
  
  
  
  
  #Evaluate AUC for training data.
  if(isTRUE(accuracy_train)){
    #Not implemented yet. FUTURE? Rewrite riskRegression::Score but validate using training data.
    #Not very interesting probably.
  } else{
    AUC_train = NULL
    Brier_train = NULL
  }
  
  
  #Evaluate AUC here for prediction data if Y_pred specified. Otherwise, don't. Just don't.
  if(!is.null(Y_pred)){
    
    if(isFALSE(reg_baseline) & isFALSE(reg_long)){
      #if we didn't regularize
      preddat <- cbind(step2$Y_pred, step2$preddat)
    } else{
      #If we regularized
      #We need to have a data.frame with outcome, linpred (from glmnet) and baseline covariates (for IPCW)
      preddat <- as.data.frame(cbind(step2$Y_pred, lp_pred, step1$X_baseline_pred))
      colnames(preddat)[3] <- "linpred"
    }
    if(IPCW_vars == "all"){
      featureNames <- paste(colnames(step1$X_baseline_pred), collapse = " + ")  
    } else if( IPCW_vars == "none"){
      featureNames <- "1"  
    } else{
      if(!is.character(IPCW_vars)){
        stop("Please specify IPCW_vars as character vector of variable names.")
      }
      featureNames <- paste(IPCW_vars, collapse = " + ")
    }
    
    
    #See explanation about Brier score calculation:
    #https://www.jesseislam.com/post/brier-score/
    
    #First we predict the absolute risk for all relevant time points:
    predRisk <- predictRisk(survival_from_glm, preddat, c(landmark_time, times_pred))
    rownames(predRisk) <- rownames(X_pred)
    colnames(predRisk) <- c(landmark_time, times_pred)
    
    
    #Add lp_train to traindat to predictRisk in regularized case, do this above.
    # predRisk_train <- predictRisk(survival_from_glm, traindat, c(landmark_time, times_pred))
    # colnames(predRisk_train) <- c(landmark_time, times_pred)
    # rownames(predRisk_train) <- rownames(X_train)
    
    #First column is predicted Risk at landmark time, second at first landmark_time and so on
    #Sometimes, when the risk-adjustment factor becomes too large, predRisk[,1] might contain 1's (rounding problem in R)
    #This causes a division by 0 further on in the code. We negate this further on in current_risk
    
    
    #Add error checking to this later. When performing IPCW we might get errors again!
    Brier_pred <- rep(NA, length(times_pred))
    AUC_pred <- rep(NA, length(times_pred))
    for(i in 1:length(times_pred)){
      pred_time <- times_pred[i]
      #The current absolute risk with landmarking is given by: 
      #1 - S(t)/S(t') = 1 - (1-F(t))/(1-F(t')) with t' landmark time (always first column)
      current_risk <- 1 - (1-predRisk[,i+1])/(1-predRisk[,1])
      #Negate the problem when predRisk[,1] becomes equal to 1.
      current_risk[which(is.nan(current_risk))] <- 1
      #preddat$predrisk <- current_risk
      
      Brier_temp <- tryCatch(
        {riskRegression::Score(
        list("timepred" = current_risk),
        formula = as.formula(paste0("Surv(time, status) ~", featureNames)),
        data = preddat,
        exact = FALSE, # Do not predict at event times
        times = pred_time,
        conf.int = FALSE,
        cens.model = "km", # Method for estimating inverse probability of censoring weights:
        splitMethod = "none",
        B = 0,
        verbose = FALSE
      )},
      error = function(cond){
        return(NULL)
      }
      )
      if(!is.null(Brier_temp)){
        Brier_pred[i] <- subset(Brier_temp$Brier$score, model == "timepred")$Brier
        AUC_pred[i] <- Brier_temp$AUC$score$AUC
      } else{
        Brier_pred[i] <- NA
        AUC_pred[i] <- NA
      }
      
    }
    names(Brier_pred) <- times_pred
    
    
    
  }
  
  out <- list(prob_surv_pred = 1- predRisk,
              #prob_surv_train = 1- predRisk_train,
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


# Brier_temp <- tryCatch({suppressMessages(riskRegression::Score(object = list("Cox1" = survival_from_glm),
#                       formula = as.formula(paste0("Surv(time, status) ~", featureNames)),
#                       data = preddat,
#                       exact = FALSE, # Do not predict at event times
#                       times = times_pred,
#                       conf.int = FALSE,
#                       cens.model = "cox", # Method for estimating inverse probability of censoring weights:
#                       splitMethod = "none",
#                       B = 0,
#                       verbose = FALSE))},
#                       warning = function(cond){
#                         message("Warning from Score() function")
#                       },
#                       error = function(cond){
#                         message("Something went wrong in the riskRegression:Score() function.")
#                         message("Throwing away results from this iteration.")
#                         return(NULL)
#                       })
# Brieer_temp <<- Brier_temp
# if(!is.null(Brier_temp)){
#   #Extract Brier score from Brier_temp
#   Brier_pred <- subset(Brier_temp$Brier$score, model == "Cox1")$Brier
#   Brier_pred <- c(Brier_pred, rep(NA, length(times_pred) - length(Brier_pred)))
#   names(Brier_pred) <- times_pred
#   #Extract AUC from Brier_temp
#   AUC_pred <- Brier_temp$AUC$score$AUC
#   AUC_pred <- c(AUC_pred, rep(NA, length(times_pred) - length(AUC_pred)))
#   names(AUC_pred) <- times_pred
# } else{
#   Brier_pred <- rep(NA, length(times_pred))
#   names(Brier_pred) <- times_pred
#   AUC_pred <- rep(NA, length(times_pred))
#   names(AUC_pred) <- times_pred
# }

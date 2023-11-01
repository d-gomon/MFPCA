#' @title Validate prediction model (step 4)
#' 
#' @description Validate prediction model by calculating tdAUC, Brier Score & MSE (if true cdf specified)
#' 
#' 
#' @param predRisk matrix containing the predicted cdf of subjects (rows) over times_pred (columns)
#' @param preddat data frame containing survival time, status indicator en covariates
#' @param times_pred Prediction times at which the models should be validated
#' @param reg_baseline Should baseline covariates be regularized?
#' @param reg_long Should longitudinal covariates be regularized?
#' @param IPCW_vars Character vector indicating which variables to use for IPCW.
#' Extra options are "all" (all baseline variables) and "none" (don't perform IPCW). Default is no IPCW ("none").
#' @param truecdf matrix containing the true cdf of subjects (rows) over times_pred (columns)
#' @param landmark_time Landmark time
#' @param FakeLM Are we performing Relaxed landmarking? Default = FALSE
#' 
#' @import riskRegression
#' 
#' 
#' @export
#' @keywords internal
#' 
#' 



validate_model <- function(predRisk, preddat, times_pred, IPCW_vars = c("none", "all"), 
                           truecdf = NULL, landmark_time, FakeLM, de_bug = FALSE,
                           de_bug_score = FALSE){
  
  if(identical(IPCW_vars, c("none", "all"))){
    IPCW_vars = "none"
  }
  
  if(isTRUE(de_bug)){
    score_storage <- vector(mode = "list", length = length(times_pred))
  }
  
  #CANNOT USE IPCW_vars with validate_model at the moment!!
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
  
  
  #First column is predicted Risk at landmark time, second at first landmark_time and so on
  #Sometimes, when the risk-adjustment factor becomes too large, predRisk[,1] might contain 1's (rounding problem in R)
  #This causes a division by 0 further on in the code. We negate this further on in current_risk
  
  
  if(isTRUE(de_bug_score)){
    current_risk <- 1 - (1-predRisk)/(1-predRisk[,1])
    #Negate the problem when predRisk[,1] becomes equal to 1.
    current_risk[which(is.nan(current_risk))] <- 1
    
    Brier_temp <- riskRegression::Score(
      list("timepred" = current_risk[, -1]),
      formula = as.formula(paste0("Surv(time, status) ~ 1")),
      data = preddat,
      exact = FALSE, # Do not predict at event times
      times = times_pred,
      conf.int = TRUE,
      cens.model = "km", # Method for estimating inverse probability of censoring weights:
      splitMethod = "none",
      B = 0,
      verbose = FALSE,
      contrasts = FALSE,
      null.model = FALSE
    )
    Brier_pred <- subset(Brier_temp$Brier$score, model == "timepred")$Brier
    AUC_pred <- Brier_temp$AUC$score$AUC
    score_storage <- Brier_temp
  } else{
    #Code below is in a for loop because this speeds up computation!!
    #See code even more below, there we let riskRegression:Score calculate Brier & AUC
    #over all time points, but this is significantly slower.
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
          verbose = FALSE,
          null.model = FALSE
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
      
      if(isTRUE(de_bug)){
        score_storage[[i]] <- Brier_temp
      }
      
    }
  }
  
  names(Brier_pred) <- times_pred
  names(AUC_pred) <- times_pred
  out <- list(Brier_pred = Brier_pred,
              AUC_pred = AUC_pred)
  
  #If true absolute risk is specified, calculate MSE between Landmarked Survival probability and actual landmarked survival probability
  #We have: F(t) true and \widehat{S(t)}, given as truecdf and step3$prob_surv_pred respectively.
  #So we calculate S(t|t_LM) = 1 - (1-F(t))/(1-F(t_LM)) and \hat{S}(t|t_LM) = (1-\hat{S}(t))/(1 - \hat{S}(t_LM)) respectively
  if(!is.null(truecdf)){
    truecdf <- truecdf[, match(c(landmark_time, times_pred), colnames(truecdf))]

    #Extract correct TRUE survival probabilities.
    Ft_timespred <- truecdf[as.numeric(rownames(predRisk)), match(times_pred, colnames(truecdf))]
    Ft_lm <- truecdf[as.numeric(rownames(predRisk)), match(landmark_time, colnames(truecdf))]
    #Extract predicted survival probabilities
    Sthat_timespred <- 1-predRisk[, match(times_pred, colnames(predRisk))]
    if(isTRUE(FakeLM)){
      Sthat_lm <- 1-predRisk[, match(landmark_time, colnames(predRisk))]  
    } else{
      Sthat_lm <- rep(1, nrow(Sthat_timespred))
    }
    
    # for(l in seq_along(times_pred)){
    #   at_risk <- which(step1$Y_surv_pred[,1] > times_pred[l])
    #   MSE_surv[i,l] <- mean(((1-Ft_timespred[at_risk, l])/(1-Ft_lm[at_risk]) - (Sthat_timespred[at_risk, l])/(Sthat_lm[at_risk]))^2)
    # }
    
    St_true <- (1-Ft_timespred)/(1-Ft_lm)
    St_hat <- (Sthat_timespred)/(Sthat_lm)
    MSE_surv <- colMeans((St_true - St_hat)^2)
    names(MSE_surv) <- times_pred
    out$MSE <- MSE_surv
  }
  
  if(isTRUE(de_bug) | isTRUE(de_bug_score)){
    out$score_storage <- score_storage
  }
  
  return(out)
  
}
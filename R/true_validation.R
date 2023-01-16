#' @title Calculate true validation scores if survival probability of subjects is known
#' 
#' @description Calculate true tdAUC and Brier score given inputs.
#' 
#' 
#' @param trueCDF CDF of observing an event at pred_time (matrix), with rows indicating subjects and columns the prediction times.
#' Column names should be the prediction times of interest. First column must be landmark time.
#' @param Y_surv Matrix of survival and censoring times (time, status)
#' @param landmark_time Landmark time for validation. Default = NULL will use first available time.
#'
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @importFrom survival coxph.control
#' @import prodlim
#' @import riskRegression
#' @importFrom success extract_hazard
#' 
#' 
#' @export
#' @keywords internal
#' 
#' 




true_validation <- function(trueCDF, Y_surv, landmark_time = NULL){
  
  
  
  asdt <- data.frame(asd = rep(1, nrow(trueCDF)))
  preddat <- as.data.frame(cbind(asdt, Y_surv))
  times_pred <- as.numeric(colnames(trueCDF))
  
  if(is.null(landmark_time)){
    lmidx = 1
  } else{
    lmidx = which.max(times_pred >= landmark_time)
    times_pred <- times_pred[lmidx:length(times_pred)]
    trueCDF <- trueCDF[, lmidx:ncol(trueCDF)]
  }
  
  Brier_pred <- rep(NA, length(times_pred) - 1)
  AUC_pred <- rep(NA, length(times_pred) - 1)
  for(i in 2:length(times_pred)){
    pred_time <- times_pred[i]
    #The current absolute risk with landmarking is given by: 
    #1 - S(t)/S(t') = 1 - (1-F(t))/(1-F(t')) with t' landmark time (always first column)
    current_risk <- 1 - (1-trueCDF[,i])/(1-trueCDF[,1])
    #Negate the problem when predRisk[,1] becomes equal to 1.
    current_risk[which(is.nan(current_risk))] <- 1
    #preddat$predrisk <- current_risk
    
    Brier_temp <- tryCatch(
      {riskRegression::Score(
        list("timepred" = current_risk),
        formula = as.formula(paste0("Surv(time, status) ~", "1")),
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
      Brier_pred[i-1] <- subset(Brier_temp$Brier$score, model == "timepred")$Brier
      AUC_pred[i-1] <- Brier_temp$AUC$score$AUC
    } else{
      Brier_pred[i-1] <- NA
      AUC_pred[i-1] <- NA
    }
    
  }
  names(Brier_pred) <- times_pred[-1]
  names(AUC_pred) <- times_pred[-1]
  
  return(list(AUC_pred = AUC_pred,
              Brier_pred = Brier_pred))
}
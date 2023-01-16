


Evaluate_tdauc <- function(
    surv.new, 
    linpred, 
    T.start, 
    deltaT) {
  
  res.tp <- vector(mode = "list", length = length(deltaT))
  res.fp <- vector(mode = "list", length = length(deltaT))
  res.tdauc <- vector(mode = "numeric", length = length(deltaT))
  
  # This line is no more required after shifting the time scale
  # # Remove deltaT before landmark time
  # deltaT <- deltaT[deltaT > T.start]
  
  for (j in 1:length(deltaT)) {
    
    predict.time <- deltaT[j]
    
    #    print(predict.time)
    event.times <- surv.new$time[surv.new$event == 1]
    if (all(!(event.times <= predict.time))) {
      mess <- paste(
        "No event (surv.new$event == 1) is observed between landmark time", T.start, 
        "and prediction time", predict.time)
      warning(mess)
      # Store no result
      res.tdauc[j] <- NA
      res.tp[[j]] <- NA
      res.fp[[j]] <- NA
    } else {
      temp <- survivalROC::survivalROC(
        Stime = surv.new$time, # Event time or censoring time for subjects
        status = surv.new$event, # Indicator of status, 1 if death or event, 0 otherwise
        marker = linpred, # Predictor or marker value
        entry = NULL, # Entry time for the subjects, default is NULL
        predict.time = predict.time, # Time point of the ROC curve
        cut.values = NULL, # marker values to use as a cut-off for calculation of sensitivity and specificity
        method = "NNE", 
        span = 0.25 * nrow(surv.new)^(-0.2) # small span yield moderate smoothing, how to select?
      )
      # Store result
      res.tdauc[j] <- temp$AUC
      res.tp[[j]] <- temp$TP
      res.fp[[j]] <- temp$FP
    }
  }
  
  return(list(
    tp = res.tp,
    fp = res.fp,
    tdauc = res.tdauc
  ))
}
#' @title Cross validation procedure for MFPCCox
#' @description Cross validate MFPCCox procedure
#' 
#' 
#' @param mdat A multiFunData object containing longitudinal covariates
#' @param X_baseline A data frame containing baseline covariates
#' @param Y_surv A 2 column data frame containing (time, status) survival info
#' @param landmark_time Time at which to landmark the data
#' @param FakeLM Should Fake (Relaxed) landmarking be performed? Default = False.
#' @param times_pred Times at which to predict the survival prob/AUC
#' @param M Number of dimensions used to approximate the process in mFPCA.
#' @param uniExpansions Univariate expansions to use in MFPCA, see MFPCA().
#' @param age Variable to perform Age Adjustment with.
#' @param AgeDM Should an extra step of Age Adjustment be performed? (Outdated, keep at FALSE)
#' @param type What measure should be used to summarize the longitudinal processes? Choice between mFPCA scores (scores), Area under Curve (AUC), 
#' predicted process (pp) and univariate scores (uscores). Default = scores.
#' @param n_folds Number of folds for cross validation.
#' @param seed Seed for random number generation
#' @param verbose Should progress messages be displayed?
#' @param reg_baseline Should baseline variables be LASSO regularized? Default = FALSE.
#' @param reg_long Should longitudinal variables be LASSO regularized? Default = TRUE.
#' @param IPCW_vars Vector indicating which variables should be used for Inverse Probability of Censoring Weights (IPCW) in determining validation scores.
#' Default is "none", meaning no IPCW will be performed.
#' @param de_bug Development option. Returns intermediate steps in the output. Default = FALSE.
#' @param trueSurv A matrix indicating the true Survival probabilities. Allows for comparing predicted and true survival probability
#' 
#' 
#' @import pbapply
#' @import parallel
#' @import doSNOW
#' @import devtools
#' 
#' @export
#' @keywords internal
#' 

cv_mfpccox <- function(mFData, X_baseline, Y_surv, landmark_time = NULL, FakeLM = FALSE,
                       times_pred = NULL, M = 50, uniExpansions = NULL, 
                       age = NULL, AgeDM = FALSE, type = c("scores", "AUC", "pp", "uscores"), 
                       n_folds = 10, seed = 01041996, verbose = FALSE,
                       reg_baseline = FALSE, reg_long = TRUE, IPCW_vars = "none",
                       de_bug = FALSE, trueSurv){
  set.seed(seed)
  
  #When using "Fake" Landmarking, do not landmark training data, only landmark test data (see below)
  if(isTRUE(FakeLM)){
    landmark_time_temp <- NULL
  } else{
    landmark_time_temp <- landmark_time
  }
  
  #First we landmark the data
  step1 <- landmark_data(time = landmark_time_temp, mFData_train = mFData, 
                         X_baseline_train = X_baseline, Y_surv_train = Y_surv,
                         age_train = age)
  mFData <- step1$mFData_train
  X_baseline <- step1$X_baseline_train
  Y_surv <- step1$Y_surv_train
  age <- step1$age_train
  
  #We perform n_fold cross validation. 
  n_obs <- nObs(mFData)
  t_id <- 1:n_obs
  t_id <- sample(t_id, n_obs)
  folds <- cut(t_id,breaks = n_folds,labels=FALSE)
  
  AUC_temp <- matrix(NA, nrow = n_folds, ncol = length(times_pred))
  Brier_temp <- matrix(NA, nrow = n_folds, ncol = length(times_pred))
  
  if(!missing(trueSurv)){
    MSE_surv <- matrix(NA, nrow = n_folds, ncol = length(times_pred))
  }
  
  for(i in 1:n_folds){
    message(paste("Working on fold", i))
    testIndexes <- which(folds==i,arr.ind=TRUE)
    
    message("Step 1/3")
    #Combine Training and prediction data in a list to feed to get_mscores
    step1 <- list(time = landmark_time, 
                  mFData_train = mFData[-testIndexes,],
                  X_baseline_train = X_baseline[-testIndexes, , drop = FALSE], 
                  Y_surv_train = Y_surv[-testIndexes,], 
                  age_train = age[-testIndexes], 
                  mFData_pred = mFData[testIndexes,], 
                  X_baseline_pred = X_baseline[testIndexes, , drop = FALSE], 
                  Y_surv_pred = Y_surv[testIndexes,], 
                  age_pred = age[testIndexes])
    
    
    #If we want to "Fake" Landmark, we only landmark the test data, but do so by setting all observations 
    #that are landmarked to NA.
    if(isTRUE(FakeLM)){
      
      if(is.null(step1$X_baseline_pred)){
        temp_X_base_pred <- NULL
      } else{
        temp_X_base_pred <- step1$X_baseline_pred
      }
      
      lm_temp <- landmark_data_fake(time = landmark_time, mFData = step1$mFData_pred, 
                                    X_baseline = temp_X_base_pred, Y_surv = step1$Y_surv_pred,
                                    age = step1$age_pred)
      step1$mFData_pred <- lm_temp$mFData
      step1["X_baseline_pred"] <- list(lm_temp$X_baseline)
      step1$Y_surv_pred <- lm_temp$Y_surv
      step1["age_pred"] <- list(lm_temp$age)
    }
    #Here something like: if Correct_LM = FALSE
    #asd <- landmark_data(step1$mFData_pred, age_pred, Y_surv_pred)
    #step1$mFData_pred <- asd$..., step1$Y_surv_pred <- asd$.... etc.
    #To remove information for test patients after landmark time, and also of patients that had failure before landmark time.
    
    if(!is.null(age)){
      if(isFALSE(AgeDM)){
        for(j in 1:length(uniExpansions)){
          uniExpansions[[j]]$age.bl <- step1$age_train
        }
      }
    }
    
    if(isTRUE(AgeDM)){
      if(is.null(age)){
        stop("Please provide age for de-meaning step.")
      }
      message("Step 1.5/3")
      DeMean_temp <- DeMean_Var(mFData = step1$mFData_train, var = step1$age_train)
      step1$mFData_train <- DeMean_temp$mFData
      step1$mFData_pred <- DeMean_test(step1$mFData_pred, var = step1$age_pred, meanFuns = DeMean_temp$meanFuns)$mFData
      step1["age_pred"] <- list(NULL)
      step1["age_train"] <- list(NULL)
    }
    
    message("Step 2/3")
    #Perform Step 2 of procedure. 
    step2 <- get_mscores(step1 = step1, M = M, 
                         uniExpansions = uniExpansions,
                         type = type, verbose = verbose)
    
    message("Step 3/3")
    step3 <- predict_surv(step2 = step2, times_pred = times_pred, 
                          reg_baseline = reg_baseline, reg_long = reg_long, IPCW_vars = IPCW_vars)
    
    #Calculate AUC for current fold
    AUC_temp[i, ] <- step3$AUC_pred
    Brier_temp[i, ] <- step3$Brier_pred
    
    #If true survival matrix is specified, calculate MSE between predicted and actual
    if(!missing(trueSurv)){
      MSE_surv[i,] <- sapply(1:ncol(step3$prob_surv_pred), function(x) mean((step3$prob_surv_pred[,x] - trueSurv[as.numeric(rownames(step3$prob_surv_pred)),x])^2))[-1]
    }
  }
  #Average AUC over folds:
  AUC <- colMeans(AUC_temp, na.rm = TRUE)
  names(AUC) <- times_pred
  Brier <- colMeans(Brier_temp, na.rm = TRUE)
  names(Brier) <- times_pred
  out <- list(AUC_pred = AUC,
              Brier_pred = Brier,
              landmark_time = step3$landmark_time)
  
  
  if(isTRUE(de_bug)){
    out$step1 <- step1
    out$step2 <- step2
    out$step3 <- step3
  }
  if(!missing(trueSurv)){
    MSE <- colMeans(MSE_surv, na.rm = TRUE)
    names(MSE) <- times_pred
    out$MSE <- MSE
  }
  return(out)
}


#' @title Repeated Cross Validated Regularized AA MFPCACox procedure.
#' 
#' @description Same as cv_mfpccox(), but repeated instead of single cross validation.
#' 
#' @inheritParams cv_mfpccox
#' @param n_reps Number of repeats of the cross validation procedure
#' @param displaypb Should a progress bar be displayed?
#' @param n_cores How many cores should be used?
#' 
#' @export
#' 
#' 
#' @keywords internal

rcv_mfpccox <- function(mFData, X_baseline, Y_surv, landmark_time = NULL, FakeLM = FALSE,
                       times_pred = NULL, M = 50, uniExpansions = NULL, AgeDM = FALSE,
                       age = NULL, type = c("scores", "AUC", "pp", "uscores"), n_reps = 10, 
                       n_folds = 10, seed = 01041996, displaypb = FALSE, 
                       n_cores = 1, verbose = FALSE,
                       reg_baseline = FALSE, reg_long = TRUE, IPCW_vars = "none", trueSurv){
  
  if(n_cores > 1){
    #Create a cluster for parallel computing and error check
    real_cores <- detectCores()
    if(n_cores > real_cores){
      warning(paste0("More cores requested (", n_cores ,") than detected (", real_cores,") by R. \n Proceed at own risk."))
    }
    cl <- makeCluster(n_cores)
    parallel::clusterEvalQ(cl, library(MFPCA))
    # clusterExport(cl, c("mFData", "X_baseline", "Y_surv", "landmark_time", "times_pred",
    #                     "M", "uniExpansions", "age", "AgeDM", "type", "n_folds", "seed",
    #                     "verbose", "reg_baseline", "reg_long", "predict_surv", "IPCW_vars", "trueSurv"),
    #               envir=environment())
    # clusterExport(cl, c("MFPCA", "univDecomp", "fpcaBasis", "PACE",
    #                     "cv_mfpccox", "landmark_data", "landmark_data_fake",
    #                     "DeMean_Var", "DeMean_Var_FunData", "DeMean_test", "DeMean_test_funData",
    #                     "get_mscores", "predict_uscores"))
  } else{
    cl <- NULL
  }
  if(displaypb){
    pboptions(type = "timer")
  } else{
    pboptions(type = "none")
  }
  set.seed(seed)
  out <- pblapply(1:n_reps, function(x){
    devtools::load_all()
    #suppressMessages(
    cv_mfpccox(mFData, X_baseline, Y_surv, landmark_time = landmark_time, FakeLM = FakeLM,
    times_pred = times_pred, M = M, uniExpansions = uniExpansions, 
    age = age, AgeDM = AgeDM, type = type, n_folds = n_folds,
    seed = seed + x, verbose = verbose, reg_baseline = reg_baseline,
    reg_long = reg_long, IPCW_vars = IPCW_vars, trueSurv = trueSurv)
    #)
    }, cl = cl)
  if(!is.null(cl)){
    stopCluster(cl)  
  }
  AUC_meas <- matrix(NA, nrow = length(times_pred), ncol = n_reps)
  Brier_meas <- matrix(NA, nrow = length(times_pred), ncol = n_reps)
  for(i in 1:length(out)){
    AUC_meas[,i] <- out[[i]]$AUC_pred
    Brier_meas[,i] <- out[[i]]$Brier_pred
  }
  AUC <- rowMeans(AUC_meas)
  names(AUC) <- times_pred
  Brier <- rowMeans(Brier_meas)
  names(Brier) <- times_pred
  if(!missing(trueSurv)){
    MSE_meas <- matrix(NA, nrow = length(times_pred), ncol = n_reps)
    for(i in 1:length(out)){
      MSE_meas[,i] <- out[[i]]$MSE
    }
    MSE <- rowMeans(MSE_meas)
    names(MSE) <- times_pred
  }
  final <- list(AUC_pred = AUC,
                Brier_pred = Brier)
  if(!missing(trueSurv)){
    final$MSE <- MSE
  }
  class(final) = "rcv_mfpccox"
  final
}










#' @title Cross validation procedure for MFPCCox
#' 
#' 
#' @param mdat A multiFunData object containing longitudinal covariates
#' @param X_baseline A data frame containing baseline covariates
#' @param Y_surv A 2 column data frame containing (time, status) survival info
#' @param landmark_time Time at which to landmark the data
#' @param times_pred Times at which to predict the survival prob/AUC
#' @param M Number of dimensions to use in mFPCA
#' @param uniExpansions Univariate expansions to use in MFPCA
#' @param age Variable to perform alternative de-meaning with.
#' @param n_folds Number of folds for cross validation.
#' @param seed Seed for random number generation

cv_mfpccox <- function(mFData, X_baseline, Y_surv, landmark_time = 0, 
                       times_pred = NULL, M = 50, uniExpansions = NULL, 
                       age = NULL, n_folds = 10, seed = 01041996){
  set.seed(seed)
  
  #First we landmark the data
  step1 <- landmark_data(time = landmark_time, mFData_train = mFData, 
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
  
  for(i in 1:n_folds){
    message(paste("Working on fold", i))
    testIndexes <- which(folds==i,arr.ind=TRUE)
    
    message("Step 1/3")
    #Combine Training and prediction data in a list to feed to get_mscores
    step1 <- list(time = landmark_time, 
                  mFData_train = mFData[-testIndexes,],
                  X_baseline_train = X_baseline[-testIndexes,], 
                  Y_surv_train = Y_surv[-testIndexes,], 
                  age_train = age[-testIndexes], 
                  mFData_pred = mFData[testIndexes,], 
                  X_baseline_pred = X_baseline[testIndexes,], 
                  Y_surv_pred = Y_surv[testIndexes,], 
                  age_pred = age[testIndexes])
    
    if(!is.null(age)){
      for(j in 1:length(uniExpansions)){
        uniExpansions[[j]]$age.bl <- step1$age_train
      }
    }
    message("Step 2/3")
    #Perform Step 2 of procedure. 
    step2 <- get_mscores(step1 = step1, M = M, 
                         uniExpansions = uniExpansions, verbose = TRUE)

    
    message("Step 3/3")
    step3 <- predict_surv(step2 = step2, times_pred = times_pred)
    
    #Calculate AUC for current fold
    AUC_temp[i, ] <- step3$AUC_pred
    #Brier_temp[i, ] <- step3$Brier_pred
  }
  #Average AUC over folds:
  AUC <- colMeans(AUC_temp, na.rm = TRUE)
  names(AUC) <- times_pred
  Brier <- colMeans(Brier_temp, na.rm = TRUE)
  names(Brier) <- times_pred
  return(list(AUC_pred = AUC,
              Brier_pred = Brier))
}




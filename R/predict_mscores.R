#' @title Predict multivariate scores for multiFunData (step 2)
#' 
#' @param mFData_train Training data in multiFunData format. Can also supply an 
#' MFPCAfit already constructed on training data. In that case, uniExpansions 
#' does not have to be specified, but will instead be extracted from mFData_train.
#' @param mFData_pred Testing data in multiFunData format
#' @param age_train Age of patients at baseline for training data.
#' @param age_pred Age of patients at baseline for testing data
#' @param M Number of dimensions to use in mFPCA
#' @param uniExpansions Univariate expansions to use in mFPCA
#' @param verbose Should intermediary results be printed?
#' @param landmark_time Time at which to landmark
#' @param Y_train Survival info (time, event) for train data. Only required in case of landmarking!
#' @param type What summary should be returned? Default is the multivariate scores ("scores"),
#' other choices are Area under curve of predicted process ("AUC") or predicted process at
#' landmark time ("pp"). Default is "scores".
#' 
#' 
#' @param step1 Output of landmark_data()
#' @param M Number of dimensions to use in mFPCA
#' @param uniExpansions Univariate expansions to use in mFPCA
#' @param verbose Should intermediary results be printed?
#' 
#' @description Time units in Y_train and mFData_train$argvals have to be the same, otherwise cannot do landmarking!!!
#'


get_mscores <- function(step1, type = c("scores", "AUC", "pp", "uscores"), M = NULL, 
                        uniExpansions = NULL, verbose = FALSE){
  if(!missing(step1)){
    list2env(step1, envir = environment())
    landmark_time <- time
    Y_train <- Y_surv_train
  }
  
  if(length(type) > 1){
    type <- "scores"
  }
  
  #Fit MFPCA on Training data
  if(inherits(mFData_train, "multiFunData")){
    if(length(mFData_train) != length(mFData_pred)){
      stop("Training and prediction data should be of same dimensions.")
    }
    #First we fit the mFPCA model to calculate scores and vectors for training_data
    if(type != "scores"){
      MFPCAfit <- MFPCA(mFData_train, M = M, uniExpansions = uniExpansions, verbose = verbose,
                        age.bl = age_train, fit = TRUE)
    } else{
      MFPCAfit <- MFPCA(mFData_train, M = M, uniExpansions = uniExpansions, verbose = verbose,
                        age.bl = age_train)  
    }
    uniExpansions <- MFPCAfit$uniExpansions
  } 
  
  #Use fitted MFPCA on training data to predict scores for prediction data.
  if(!is.null(mFData_pred)){
    #Data preparation and checks
    c <- MFPCAfit$vectors
    p <- length(mFData_pred)
    #Old number of functional principal components
    M <- ncol(c)
    #Then we need to predict the uFPCA scores for mFData_pred from mFData_train fit.
    #For this we need to re-fit uFPCA for mFData_train
    UFPCA_pred <- lapply(1:p, FUN = function(x){predict_uscores(UFPCAfit = MFPCAfit$uniBasis[[x]],
                                                                uData_pred = mFData_pred[[x]],
                                                                age.bl.pred = age_pred,
                                                                age.bl = age_train)})
    
    
    allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{UFPCA_pred[[j]]$scores}
    #Columns of c contain eigenvectors. See ?eigen for confirmation
    mscores <- allScores %*% c
    #normalization factor
    mscores <- as.matrix(mscores %*% diag(sqrt(MFPCAfit$values) * MFPCAfit$normFactors, nrow = M, ncol = M)) # normalization
    namesList <- lapply(mFData_pred, names)
    rownames(mscores) <- namesList[[1]]
  } else{
    mscores = NULL
  }
  
  
  #Finish later
  if(type == "scores"){
    #Construct the predictors for Step 3 of our method.
    long_sum_train <- MFPCAfit$scores
    long_sum_pred <- mscores
  } else if(type == "AUC"){
    #We need to calculate area under the curve.
    #MFPCAfit will have the truncated Karhuhen Loeve in MFPCAfit$fit for training data
    #Calculate it for test data as well, using predict.MFPCAfit
    trunc_KL_pred <- predict(MFPCAfit, scores = mscores)
    trunc_KL_train <- MFPCAfit$fit
    
    #Integrate to get AUC
    long_sum_train <- sapply(trunc_KL_train, integrate)
    long_sum_pred <- sapply(trunc_KL_pred, integrate)
    
  } else if(type == "pp"){
    trunc_KL_pred <- predict(MFPCAfit, scores = mscores)
    trunc_KL_train <- MFPCAfit$fit
    
    fObs_pred <- subset(trunc_KL_pred, 
                        argvals = list(tail(trunc_KL_pred[[1]]@argvals[[1]], 1)))
    fObs_train <- subset(trunc_KL_train, 
                         argvals = list(tail(trunc_KL_train[[1]]@argvals[[1]], 1)))
    
    long_sum_train <- matrix(NA, nrow = nObs(fObs_train), ncol = length(fObs_train))
    long_sum_pred <- matrix(NA, nrow = nObs(fObs_pred), ncol = length(fObs_pred))
    for(i in 1:length(fObs_pred)){
      long_sum_train[,i] <- fObs_train[[i]]@X
      long_sum_pred[,i] <- fObs_pred[[i]]@X
    }
  } else if(type == "uscores"){
    allScores_train <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{MFPCAfit$uniBasis[[j]]$scores}
    long_sum_pred <- allScores
    long_sum_train <- allScores_train
  }
  
  #Add the summaries to the baseline covariates
  X_train <- cbind(long_sum_train, X_baseline_train)
  X_train <- model.matrix(~.-1, X_train)
  
  X_pred <- cbind(long_sum_pred, X_baseline_pred)
  X_pred <- model.matrix(~.-1, X_pred)

  
  return(list(X_train = X_train,
              Y_train = Y_train,
              X_pred = X_pred,
              Y_pred = Y_surv_pred,
              landmark_time = landmark_time,
              scores_train = MFPCAfit$scores,
              scores_pred = mscores,
              MFPCAfit = MFPCAfit,
              step1 = step1))
}

predict_mscores <- function(mFData_train, mFData_pred = NULL, 
                            age.bl = NULL, age.bl.pred = NULL, 
                            M = NULL, uniExpansions = NULL, verbose = FALSE,
                            landmark_time = NULL, Y_train = NULL){
  
  #Some input checks
  #If mFData_train is specified as DATA, instead of mFPCA output: run MFPCA
  if(inherits(mFData_train, "multiFunData")){
    if(length(mFData_train) != length(mFData_pred)){
      stop("Training and prediction data should be of same dimensions.")
    }
    if(!is.null(landmark_time)){
      if(is.null(Y_train)){
        stop("Cannot landmark without survival information. Please specify 'Y_train'.")
      }
      #Determine which argvals to use for landmarking by cutting of at landmark_time
      max_id <- max(which(mFData_train[[1]]@argvals[[1]] <= landmark_time))
      #Determine which observations to use by eliminating people with failures before landmark_time (requires Y_train)
      which_alive <- which(Y_train[, "time"] > landmark_time)
      #original num of observations
      nrow_orig <- nrow(mFData_train[[1]]@X)
      #Subset the data according to above 2 conditions.
      mFData_train <- subset(mFData_train, obs = which_alive, argvals = list(mFData_train[[1]]@argvals[[1]][1:max_id]))
      warning(paste("Retained ",nrow(mFData_train[[1]]@X)*100/nrow_orig, "percent of original observations due to landmarking."))
      #If we landmarked, and want to de-mean according to age, we should adjust the age.bl vars
      #to remove baseline covs from people removed in landmarking.
      if(!is.null(age.bl)){
        age.bl <- age.bl[which_alive]
      }
    }
    #First we fit the mFPCA model to calculate scores and vectors for training_data
    MFPCAfit <- MFPCA(mFData_train, M = M, uniExpansions = uniExpansions, verbose = verbose,
                      age.bl = age.bl)
    uniExpansions <- MFPCAfit$uniExpansions
  } else if(inherits(mFData_train, "MFPCAfit")){ #If already MFPCAfit, extract info
    if(!is.null(landmark_time)){
      warning("Cannot landmark a MFPCAfit, make sure the fit is landmarked.")
      #Use subset(mFData_train, failtimes < landmark, argvals < landmark)
    }
    MFPCAfit <- mFData_train
    uniExpansions <- MFPCAfit$uniExpansions
  }
  if(!is.null(mFData_pred)){
    if(!is.null(landmark_time)){
      mFData_pred <- subset(mFData_pred, argvals = list(mFData_train[[1]]@argvals[[1]][1:max_id]))
    }
    #Data preparation and checks
    c <- MFPCAfit$vectors
    p <- length(mFData_pred)
    #Old number of functional principal components
    M <- ncol(c)
    
    #Then we need to predict the uFPCA scores for mFData_pred from mFData_train fit.
    #For this we need to re-fit uFPCA for mFData_train
    UFPCA_pred <- lapply(1:p, FUN = function(x){predict_uscores(UFPCAfit = MFPCAfit$uniBasis[[x]],
                                                                           uData_pred = mFData_pred[[x]],
                                                                           age.bl.pred = age.bl.pred,
                                                                age.bl = age.bl)})
    
    
    allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{UFPCA_pred[[j]]$scores}
    #Columns of c contain eigenvectors. See ?eigen for confirmation
    mscores <- allScores %*% c
    #normalization factor
    mscores <- as.matrix(mscores %*% diag(sqrt(MFPCAfit$values) * MFPCAfit$normFactors, nrow = M, ncol = M)) # normalization
    namesList <- lapply(mFData_pred, names)
    rownames(mscores) <- namesList[[1]]
  } else{
    mscores = NULL
  }
  return(list(mscores_test = mscores,
              mscores_train = MFPCAfit$scores,
              MFPCAfit = MFPCAfit))
}



predict_uscores <- function(UFPCAfit, uData_pred, age.bl.pred = NULL, age.bl = NULL, nbasis_mu = 10){
  #Preliminaries
  Y.pred <- uData_pred@X
  X <- uData_pred@argvals[[1]]
  I.pred = NROW(Y.pred)
  D = NCOL(Y.pred)
  
  #Extract data from UFPCAfit
  npc <- dim(UFPCAfit$scores)[2]
  Z <- efunctions <- t(UFPCAfit$functions@X)
  sigma2 <- UFPCAfit$sigma2
  evalues <- UFPCAfit$evalues
  gam0 <- UFPCAfit$gam0
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  
  if(!is.null(age.bl.pred) ){
    age.bl_order_pred <- age.bl.pred - min(age.bl)
    d.vec_adj_pred <- rep(X, each = I.pred) + age.bl_order_pred
    #predict values at all distinct ages
    mu.pred = mgcv::predict.gam(gam0, newdata = data.frame(d.vec_adj = sort(unique(d.vec_adj_pred))))
    #create matrix with rows = subjects and columns = number of observations
    #match subject means with correct entry in prediction
    mu.pred_mat <- matrix(mu.pred[match(d.vec_adj_pred, sort(unique(d.vec_adj_pred)))], I.pred, D, byrow = TRUE)
    Y.tilde = Y.pred - mu.pred_mat
    mu_mat <- mu.pred_mat
  } else{
    #length(mu) = length(X), because the range of values is equal everywhere.
    mu = mgcv::predict.gam(gam0, newdata = data.frame(d.vec = X))
    mu_mat <- t(as.matrix(mu))
    Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
  }
  fit = matrix(0, nrow = I.pred, ncol = D)
  scores = matrix(NA, nrow = I.pred, ncol = npc)
  #Create variables to keep track of scores which could not be estimated
  scores_fail <- c()
  scores_fail_partial <- c()
  # no calculation of confidence bands, no variance matrix
  for (i.subj in seq_len(I.pred)) {
    obs.points = which(!is.na(Y.pred[i.subj, ]))
    if (sigma2 == 0 & length(obs.points) < npc) {
      if(length(obs.points) == 0){
        scores_fail <- c(scores_fail, i.subj)
        scores[i.subj, ] <- rep(0, npc)
        if(!is.null(age.bl.pred)){
          fit[i.subj, ] <- mu_mat[i.subj,]
        } else{
          fit[i.subj, ] <- t(as.matrix(mu))
        }
        next  
      } else{
        scores_fail_partial <- c(scores_fail_partial, i.subj)
        npc_temp <- length(obs.points)
        Zcur = matrix(Z[obs.points, 1:npc_temp], nrow = length(obs.points), 
                      ncol = npc_temp)
        ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv[1:npc_temp, 1:npc_temp])
        scores_temp = ZtZ_sD.inv %*% crossprod(Zcur, Y.tilde[i.subj, obs.points])
        scores[i.subj, ] = c(scores_temp, rep(0, npc - npc_temp))
        if(!is.null(age.bl.pred)){
          fit[i.subj, ] <- mu_mat[i.subj,] + tcrossprod(scores[i.subj, ], efunctions)
        } else{
          fit[i.subj, ] = t(as.matrix(mu)) + tcrossprod(scores[i.subj, ], efunctions)
        }
        next
      }
    }
    Zcur = matrix(Z[obs.points, ], nrow = length(obs.points),
                  ncol = dim(Z)[2])
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
    scores[i.subj, ] = ZtZ_sD.inv %*% crossprod(Zcur, Y.tilde[i.subj, obs.points])
    if(!is.null(age.bl.pred)){
      fit[i.subj, ] <- mu_mat[i.subj,] + tcrossprod(scores[i.subj, ], efunctions)
    } else{
      fit[i.subj, ] = t(as.matrix(mu)) + tcrossprod(scores[i.subj, ], efunctions)
    }
  }
  if(!is.null(age.bl.pred)){
    mu <- mu_mat
  }
  ret.objects = c("fit", "scores", "mu", "efunctions", "evalues",
                  "npc", "sigma2", "scores_fail", "gam0") # add sigma2 to output
  if(length(scores_fail) != 0){
    warning(paste("Scores could not be estimated for", length(scores_fail), "subjects and were set to zero instead.
                  Indices can be found in scores_fail. Lowering pve or npc can alleviate this problem."))
  }
  if(length(scores_fail_partial) != 0){
    warning(paste("Scores could only be partially estimated for ", length(scores_fail_partial), "subjects. Scores which could not be estimated were set to zero.
                  Indices can be found in scores_fail_partial. Lowering pve or npc can alleviate this problem."))
  }
  ret = lapply(seq_len(length(ret.objects)), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  return(ret)
}

predict_mscores_old <- function(MFPCAfit, mFData_pred, uniExpansions = NULL, M = NULL){
  
  #mdat, M = 50, uniExpansions = uniexp,
  #fit = TRUE, verbose = TRUE
  
  
  #Data preparation and checks
  c <- MFPCAfit$vectors
  p <- length(mFData_pred)
  #Old number of functional principal components
  M_old <- ncol(c)
  
  #Use the same uniExpansions for new data as for old data
  #If this is not done, you cannot consistently generate predictions for new data
  #as the dimension of the principal components can be different
  if(is.null(uniExpansions)){
    if(!is.null(MFPCAfit$uniExpansions)){
      uniExpansions <- MFPCAfit$uniExpansions
    } else{
      uniExpansions <- rep(list(list(type = "uFPCA", pve = 0.99)), p) 
    }
  }
  #calculate univariate FPCA for the new data.
  #Returns a list of length(mFData_pred) with each component a uFPCA output (.PACE)
  uniBasis <- mapply(function(expansion, data){do.call(univDecomp, c(list(funDataObject = data), expansion))},
                     expansion = uniExpansions, data = mFData_pred, SIMPLIFY = FALSE)
  #determine maximal number of principal components
  npc <- vapply(uniBasis, function(x){dim(x$scores)[2]}, FUN.VALUE = 0) # get number of univariate basis functions
  if(sum(npc) < M_old){
    stop("Not enough principal components calculated during uFPCA, please increase pve or specify 'uniExpansions = NULL'.")
  }
  if(is.null(M)){
    M <- M_old
  } else if(M > sum(npc)){
    M <- sum(npc)
    warning("Function MFPCA: total number of univariate basis functions is smaller than given M. M was set to ", sum(npc), ".")
  } 
  
  allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{uniBasis[[j]]$scores}
  #Columns of c contain eigenvectors. See ?eigen for confirmation
  mscores <- allScores %*% c
  #normalization factor
  mscores <- as.matrix(mscores %*% diag(sqrt(MFPCAfit$values) * MFPCAfit$normFactors, nrow = M, ncol = M)) # normalization
  
  return(mscores)
}
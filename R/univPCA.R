## functional PCA: consider X-values as additional input
# slightly adapted version of fpca.sc-function in refund-package
# Internal function, called by PACE in context of funData objects

#' Calculate univariate functional PCA
#' 
#' This function is a slightly adapted version of the 
#' \code{fpca.sc} function in the \strong{refund} package for 
#' calculating univariate functional principal components based on a smoothed 
#' covariance function. The smoothing basis functions are penalized splines.
#' 
#' @param X A vector of xValues.
#' @param Y A matrix of observed functions (by row).
#' @param Y.pred A matrix of functions (by row) to be approximated using the 
#'   functional principal components. Defaults to \code{NULL}, i.e. the 
#'   prediction is made for the functions in \code{Y}.
#' @param nbasis An integer, giving the number of B-spline basis to use. 
#'   Defaults to \code{10}.
#' @param pve A value between 0 and 1, giving the percentage of variance 
#'   explained in the data by the functional principal components. This value is
#'   used to choose the number of principal components. Defaults to \code{0.99}
#' @param npc The number of principal components to be estimated. Defaults to 
#'   \code{NULL}. If given, this overrides \code{pve}.
#' @param makePD Logical, should positive definiteness be enforced for the 
#'   covariance estimate? Defaults to \code{FALSE}.
#' @param cov.weight.type The type of weighting used for the smooth covariance
#'   estimate. Defaults to \code{"none"}, i.e. no weighting. Alternatively, 
#'   \code{"counts"} (corresponds to \code{fpca.sc} in \strong{refund}) weights the pointwise estimates of the covariance function
#'   by the number of observation points.
#'   
#' @return \item{fit}{The approximation of \code{Y.pred} (if \code{NULL}, the 
#'   approximation of \code{Y}) based on the functional principal components.} 
#'   \item{scores}{A matrix containing the estimated scores (observations by 
#'   row).} \item{mu}{The estimated mean function.} \item{efunctions}{A matrix 
#'   containing the estimated eigenfunctions (by row).} \item{evalues}{The 
#'   estimated eigenvalues.} \item{npc}{The number of principal comopnents that 
#'   were calculated.} \item{sigma2}{The estimated variance of the measurement 
#'   error.}  \item{estVar}{The estimated smooth variance function of the data.}
#'   
#' @seealso \code{\link{PACE}}
#'   
#' @references Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009). 
#'   Multilevel functional principal component analysis. Annals of Applied 
#'   Statistics, 3, 458--488. Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). 
#'   Functional data analysis for sparse longitudinal data. Journal of the 
#'   American Statistical Association, 100, 577--590.
#'   
#' @importFrom mgcv gam predict.gam s te
#'   
#' @keywords internal
.PACE <- function(X, Y, Y.pred = NULL, age.bl = NULL, age.bl.pred = NULL, nbasis = 10, nbasis_mu = NULL, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
{
  if(is.null(nbasis_mu)){
    nbasis_mu <- nbasis
  }
  if (is.null(Y.pred)){
    Y.pred = Y
    Ypred_mis <- TRUE
  } else{
    Ypred_mis <- FALSE
  }
    
  D = NCOL(Y)
  if(D != length(X)) # check if number of observation points in X & Y are identical
    stop("different number of (potential) observation points differs in X and Y!")
  I = NROW(Y)
  if(is.null(age.bl)){
    I.pred = NROW(Y.pred)
    d.vec = rep(X, each = I) # use given X-values for estimation of mu
    gam0 = mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis_mu))
    #length(mu) = length(X), because the range of values is equal everywhere.
    mu = mgcv::predict.gam(gam0, newdata = data.frame(d.vec = X))
    mu_mat <- t(as.matrix(mu))
    Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)
  } else{
    I.pred = NROW(Y.pred)
    #CREATE NEW d.vec_adj, which represents age at baseline. Basically, substract age at baseline for each subject from X
    #age.bl must be vector containing age in same time unit as funData object!!!
    age.bl_order <- age.bl - min(age.bl)
    #We create an adjusted vector, where each patient has individual starting times.
    d.vec_adj <- rep(X, each = I) + age.bl_order
    d.vec = rep(X, each = I) # use given X-values for estimation of mu
    gam0 = mgcv::gam(as.vector(Y) ~ s(d.vec_adj, k = nbasis_mu))
    #mu = mgcv::predict.gam(gam0, newdata = data.frame(d.vec_adj = d.vec_adj))
    mu = mgcv::predict.gam(gam0, newdata = data.frame(d.vec_adj = sort(unique(d.vec_adj))))
    mu_mat <- matrix(mu[match(d.vec_adj, sort(unique(d.vec_adj)))], I, D, byrow = FALSE)
    Y.tilde = Y - mu_mat
  }
  cov.sum = cov.count = cov.mean = matrix(0, D, D)
  for (i in seq_len(I)) {
    obs.points = which(!is.na(Y[i, ]))
    cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
    cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
  }
  G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
  diag.G0 = diag(G.0)
  diag(G.0) = NA
  row.vec = rep(X, each = D) # use given X-values
  col.vec = rep(X, D) # use given X-values
  cov.weights <- switch(cov.weight.type,
                        none = rep(1, D^2),
                        counts = as.vector(cov.count),
                        stop("cov.weight.type ", cov.weight.type, " unknown in smooth covariance estimation"))
  
  npc.0 = matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0)~te(row.vec, col.vec, k = nbasis), weights = cov.weights),
                                   newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
  npc.0 = (npc.0 + t(npc.0))/2
  # no extra-option (useSymm) as in fpca.sc-method
  if (makePD) { # see fpca.sc
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE,
                            do2eigen = TRUE, trace = options()$verbose)
      as.matrix(tmp$mat)
    }
  }
  
  ### numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Chapter 8)
  w <- funData::.intWeights(X, method = "trapezoidal")
  Wsqrt <- diag(sqrt(w))
  Winvsqrt <- diag(1/(sqrt(w)))
  V <- Wsqrt %*% npc.0 %*% Wsqrt
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values
  ###
  evalues = replace(evalues, which(evalues <= 0), 0)
  npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
  efunctions = matrix(Winvsqrt%*%eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = D, ncol = npc)
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values[seq_len(npc)]  # use correct matrix for eigenvalue problem
  cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
  ### numerical integration for estimation of sigma2
  T.len <- X[D] - X[1] # total interval length
  T1.min <- min(which(X >= X[1] + 0.25*T.len)) # left bound of narrower interval T1
  T1.max <- max(which(X <= X[D] - 0.25*T.len)) # right bound of narrower interval T1
  DIAG = (diag.G0 - diag(cov.hat))[T1.min :T1.max] # function values
  # weights
  w <- funData::.intWeights(X[T1.min:T1.max], method = "trapezoidal")
  sigma2 <- max(1/(X[T1.max]-X[T1.min]) * sum(DIAG*w, na.rm = TRUE), 0) #max(1/T.len * sum(DIAG*w), 0)
  ####
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  Z = efunctions

  if(!is.null(age.bl.pred) ){
    age.bl_order_pred <- age.bl.pred - min(age.bl.pred)
    d.vec_adj_pred <- rep(X, each = I.pred) + age.bl_order_pred
    mu.pred = mgcv::predict.gam(gam0, newdata = data.frame(d.vec_adj = sort(unique(d.vec_adj_pred))))
    mu.pred_mat <- matrix(mu.pred[match(d.vec_adj_pred, sort(unique(d.vec_adj_pred)))], I.pred, D, byrow = FALSE)
    Y.tilde = Y.pred - mu.pred_mat
    mu_mat <- mu.pred_mat
  } else if(!is.null(age.bl) & isTRUE(Ypred_mis)){
    Y.tilde = Y.pred - mu_mat
  } else{
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
        if(!is.null(age.bl)){
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
        if(!is.null(age.bl)){
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
    if(!is.null(age.bl)){
      fit[i.subj, ] <- mu_mat[i.subj,] + tcrossprod(scores[i.subj, ], efunctions)
    } else{
      fit[i.subj, ] = t(as.matrix(mu)) + tcrossprod(scores[i.subj, ], efunctions)
    }
  }
  if(!is.null(age.bl)){
    mu <- funData(argvals = sort(unique(d.vec_adj)) + min(age.bl), X = matrix(mu, nrow = 1))
    #mu <- mu_mat
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
  ret$estVar <- diag(cov.hat)
  return(ret)
}







#' Univariate functional principal component analysis by smoothed covariance
#' 
#' This function calculates a univariate functional principal components 
#' analysis by smoothed covariance based on code from 
#' \code{fpca.sc} in package \strong{refund}.
#' 
#' @section Warning: This function works only for univariate functional data 
#'   observed on one-dimensional domains.
#'   
#' @param funDataObject An object of class \code{\link[funData]{funData}} or 
#'   \code{\link[funData]{irregFunData}} containing the functional data 
#'   observed, for which the functional principal component analysis is 
#'   calculated. If the data is sampled irregularly (i.e. of class 
#'   \code{\link[funData]{irregFunData}}), \code{funDataObject} is transformed 
#'   to a \code{\link[funData]{funData}} object first.
#' @param predData  An object of class \code{\link[funData]{funData}}, for which
#'   estimated trajectories based on a truncated Karhunen-Loeve representation 
#'   should be estimated. Defaults to \code{NULL}, which implies prediction for 
#'   the given data.
#' @param nbasis An integer, representing the number of  B-spline basis 
#'   functions used for estimation of the mean function and bivariate smoothing 
#'   of the covariance surface. Defaults to \code{10} (cf. 
#'   \code{fpca.sc} in \strong{refund}).
#' @param age.bl Age at baseline for age adjustment. Default = NULL (no age adjustment)
#' @param age.bl.pred Age at baseline for age adjustment for prediction set. Default = NULL (no age adjustment)
#' @param nbasis_mu Number of basis functions to use for mean smoothing. Default = 10.
#' @param pve A numeric value between 0 and 1, the proportion of variance 
#'   explained: used to choose the number of principal components. Defaults to 
#'   \code{0.99} (cf. \code{fpca.sc} in \strong{refund}).
#' @param npc An integer, giving a prespecified value for the number of 
#'   principal components. Defaults to \code{NULL}. If given, this overrides 
#'   \code{pve} (cf. \code{fpca.sc} in \strong{refund}).
#' @param makePD Logical: should positive definiteness be enforced for the 
#'   covariance surface estimate? Defaults to \code{FALSE} (cf. 
#'   \code{fpca.sc} in \strong{refund}).
#' @param cov.weight.type The type of weighting used for the smooth covariance 
#'   estimate. Defaults to \code{"none"}, i.e. no weighting. Alternatively, 
#'   \code{"counts"} (corresponds to \code{fpca.sc} in \strong{refund}) weights the
#'   pointwise estimates of the covariance function by the number of observation
#'   points.
#'   
#' @return \item{mu}{A \code{\link[funData]{funData}} object with one 
#'   observation, corresponding to the mean function.} \item{values}{A vector 
#'   containing the estimated eigenvalues.} \item{functions}{A 
#'   \code{\link[funData]{funData}} object containing the estimated functional 
#'   principal components.} \item{scores}{An matrix of estimated scores for the 
#'   observations in \code{funDataObject}. Each row corresponds to the scores of
#'   one observation.} \item{fit}{A \code{\link[funData]{funData}} object 
#'   containing the estimated trajectories based on the truncated Karhunen-Loeve
#'   representation and the estimated scores and functional principal components
#'   for \code{predData} (if this is not \code{NULL}) or \code{funDataObject} 
#'   (if \code{predData} is \code{NULL}).} \item{npc}{The number of functional 
#'   principal components: either the supplied \code{npc}, or the minimum number
#'   of basis functions needed to explain proportion \code{pve} of the variance 
#'   in the observed curves (cf. \code{fpca.sc} in \strong{refund}).} 
#'   \item{sigma2}{The estimated measurement error variance (cf. 
#'   \code{fpca.sc} in \strong{refund}).} \item{estVar}{The estimated smooth
#'   variance function of the data.}
#'   
#' @seealso \code{\link[funData]{funData}}, 
#'   \code{\link{fpcaBasis}}, \code{\link{univDecomp}}
#'   
#' @export PACE
#'   
#' @examples
#' \donttest{
#'   oldPar <- par(no.readonly = TRUE)
#' 
#'   # simulate data
#'   sim <- simFunData(argvals = seq(-1,1,0.01), M = 5, eFunType = "Poly",
#'                     eValType = "exponential", N = 100)
#' 
#'   # calculate univariate FPCA
#'   pca <- PACE(sim$simData, npc = 5)
#' 
#'   # Plot the results
#'   par(mfrow = c(1,2))
#'   plot(sim$trueFuns, lwd = 2, main = "Eigenfunctions")
#'   # flip estimated functions for correct signs
#'   plot(flipFuns(sim$trueFuns,pca$functions), lty = 2, add = TRUE)
#'   legend("bottomright", c("True", "Estimate"), lwd = c(2,1), lty = c(1,2))
#' 
#'   plot(sim$simData, lwd = 2, main = "Some Observations", obs = 1:7)
#'   plot(pca$fit, lty = 2, obs = 1:7, add = TRUE) # estimates are almost equal to true values
#'   legend("bottomright", c("True", "Estimate"), lwd = c(2,1), lty = c(1,2))
#' 
#'   par(oldPar)
#' }
PACE <- function(funDataObject, predData = NULL, age.bl = NULL, age.bl.pred = NULL, nbasis = 10, nbasis_mu = NULL, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
{
  if(is.null(nbasis_mu)){
    nbasis_mu <- nbasis_mu
  }
  # check inputs
  if(! class(funDataObject) %in% c("funData", "irregFunData"))
    stop("Parameter 'funDataObject' must be a funData or irregFunData object.")
  if(dimSupp(funDataObject) != 1)
    stop("PACE: Implemented only for funData objects with one-dimensional support.")
  if(methods::is(funDataObject, "irregFunData")) # for irregular functional data, use funData representation
    funDataObject <- as.funData(funDataObject)
  
  if(is.null(predData))
    Y.pred = NULL # use only funDataObject
  else
  {
    if(!isTRUE(all.equal(funDataObject@argvals, predData@argvals)))
      stop("PACE: funDataObject and predData must be defined on the same domains!")
    
    Y.pred = predData@X
  }
  
  
  if(!all(is.numeric(nbasis), length(nbasis) == 1, nbasis > 0))
    stop("Parameter 'nbasis' must be passed as a number > 0.")
  if(!all(is.numeric(pve), length(pve) == 1, 0 <= pve, pve <= 1))
    stop("Parameter 'pve' must be passed as a number between 0 and 1.")
  
  if(!is.null(npc) & !all(is.numeric(npc), length(npc) == 1, npc > 0))
    stop("Parameter 'npc' must be either NULL or passed as a number > 0.")
    
  if(!is.logical(makePD))
    stop("Parameter 'makePD' must be passed as a logical.")
  
  if(!is.character(cov.weight.type))
    stop("Parameter 'cov.weight.type' must be passed as a character.")
  
  
  res <- .PACE(X = funDataObject@argvals[[1]], funDataObject@X, age.bl = age.bl, age.bl.pred = age.bl.pred, Y.pred = Y.pred,
               nbasis = nbasis, nbasis_mu = nbasis_mu, pve = pve, npc = npc, makePD = makePD,
               cov.weight.type = cov.weight.type)
  if(is.null(age.bl)){
    mu_ret <- funData(funDataObject@argvals, matrix(res$mu, nrow = 1))
  } else{
    #mu_ret <- funData(, mu_ret)
    mu_ret <- res$mu
  }
  return(list(mu = mu_ret,
              values = res$evalues,
              functions = funData(funDataObject@argvals, t(res$efunctions)),
              scores = res$scores,
              fit = funData(funDataObject@argvals, res$fit),
              npc = res$npc,
              sigma2 = res$sigma2,
              estVar = funData(funDataObject@argvals, matrix(res$estVar, nrow = 1)),
              scores_fail = res$scores_fail,
              gam0 = res$gam0
              
  ))
}

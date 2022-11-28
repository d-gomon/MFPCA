#' @title Predict absolute risk
#' 
#' @description Predict absolute risk for new data
#' 
#' 
#' @param object A 'cv.glmnet' object which was used to evaluate a 'coxnet' object.
#' @param x Input matrix, see glmnet
#' @param y Response variable, see glmnet
#' @param newdata A data frame containing predictor variable combinations for which to compute predicted event probabilities
#' @param times Times at which to predict the absolute risk
#' 
#'
#' @import prodlim
#' @import glmnet 
#'
#' @export





predict.GLMnet <- function(object, newdata, type = type, s, ...){
  rest <- list(...)
  lambda=cv=NULL
  # library(glmnet)
  # requireNamespace(c("prodlim","glmnet"))
  # predict.cv.glmnet <- utils::getFromNamespace("predict.cv.glmnet","glmnet")
  # predict.glmnet <- utils::getFromNamespace("predict.glmnet","glmnet")
  rhs <- as.formula(delete.response(object$terms))
  newdata$dummy.time=rep(1,NROW(newdata))
  newdata$dummy.event=rep(1,NROW(newdata))
  dummy.formula=stats::update.formula(rhs,"prodlim::Hist(dummy.time,dummy.event)~.")
  EHF <- prodlim::EventHistory.frame(formula=dummy.formula,
                                     data=newdata,
                                     specials = NULL,
                                     unspecialsDesign=TRUE)
  asd <<- EHF
  newdata$dummy.time = NULL
  newdata$dummy.event = NULL
  if(missing(s)){
    s <- "lambda.min"
  }
  # blank Cox object obtained with riskRegression:::coxModelFrame
  info <- object$surv_info
  if (is.null(rest$lambda) && object$cv){
    p <- predict(object$fit,newx=EHF$design,type = type, s=s)
  }
  else if (is.null(rest$lambda) && !object$cv){
    if (length(object$lambda) == 1){
      p <- predict(object$fit,newx=EHF$design,type = type, s=s)
    }
    else {
      stop("Object fitted with multiple lambdas. You must pick one lambda for predict!")
    }
  }
  else {
    if (all(rest$lambda %in% object$lambda)){
      p <- predict(object$fit,newx=EHF$design,type = type, s=s)
    }
    else {
      stop("The fitted model was not fitted with one of the lambdas that was specified in predict. ")
    }
  }
  p
}



# predict.GLMnet <- function(object, newdata, type = type, s = NULL, ...){
#   if(!inherits(newdata, "data.frame")){
#     stop("Please provide newdata as data.frame.")
#   }
#   if(is.null(s)){
#     s <- object$fit$lambda.min
#   }
#   rest <- list(...)
#   rhs <- as.formula(delete.response(object$terms))
#   xnew <- model.matrix(rhs,data=newdata)
#   print(head(xnew))
#   if (is.null(rest$lambda) && object$cv){
#     p <- predict(object$fit,newx=xnew,type = type, s=s)
#   }
#   else if (is.null(rest$lambda) && !object$cv){
#     if (length(object$lambda) == 1){
#       p <- predict(object$fit,newx=xnew,type = type, s=object$lambda)
#     }
#     else {
#       stop("Object fitted with multiple lambdas. You must pick one lambda for predictRisk!")
#     }
#   }
#   else {
#     if (all(rest$lambda %in% object$lambda)){
#       p <- predict(object$fit,newx=xnew,type = type, s=rest$lambda)
#     }
#     else {
#       stop("The fitted model was not fitted with one of the lambdas that was specified in predict.")
#     }
#   }
#   
#   p
# }








# predictRisk.cv.glmnet <- function(object,
#                                   x,
#                                   y,
#                                   newdata,
#                                   times,
#                                   s,
#                                   landmark,
#                                   ...){
#   
#   #object is a cv.glmnet object
#   #Check if Cox model was built, glmnet can also be used for other models, which are not implemented yet.
#   if(inherits(object$glmnet.fit, "coxnet")){
#     fit <- object$glmnet.fit
#   } else{
#     stop(paste0("Fitted ", class(object$glmnet.fit)[1], " object. Please fit a coxnet object using cv.glmnet. (family = 'cox')"))
#   }
#   #Want as output out, the predicted failure probability at times 1-S(times)
#   #S(t) = 1 - exp(H(t)exp(\beta X))
#   
#   
#   #For this we first need to calculate the cumulative baseline hazard.
#   #Use baseHaz_cpp from riskRegression for this.
#   
#   #Is no optimal parameter specified, use 
#   if(missing(s)){
#     s <- object$lambda.min
#   }
#   
#   #If stratified model was built, read the strata from input.
#   if(inherits(y, "stratifySurv")){
#     stop("stratifySurv models not yet supported by predictRisk.cv.glmnet.")
#     strata <- attributes(y)$strata
#   } else{
#     strata <- NULL
#   }
#   
#   #stoptimes (y), status (y), eXb (derived from x) and strata (not implemented) must be sorted by strata, stoptimes and status
#   #See https://rdrr.io/github/tagteam/riskRegression/man/baseHaz_cpp.html
#   
#   #Check how Surv() was specified for response
#   #check for 2 or 3 inputs
#   #Ordering depends on strata
#   if(all(c("start","stop") %in% colnames(y))){
#     start_stop <- TRUE
#     if(!is.null(strata)){
#       ordering <- order(strata, y[, "stop"], y[, "status"])
#     }
#     else{
#       ordering <- order(y[, "stop"], y[, "status"])
#     }
#     y <- y[ordering,]
#     x <- x[ordering,]
#     start_times <- y[, "start"]
#     stop_times <- y[, "stop"]
#   } else if("time" %in% colnames(y)){
#     start_stop <- FALSE
#     if(!is.null(strata)){
#       ordering <- order(strata, y[, "time"], y[, "status"])
#     } else{
#       ordering <- order(y[, "time"], y[, "status"])
#     }
#     y <- y[ordering,]
#     x <- x[ordering,] 
#     start_times <- rep(0, nrow(y))
#     stop_times <- y[, "time"]
#   }
#   
#   #Determine number of strata, and maximum times in each of these strata
#   if(is.null(strata)){
#     strata <- rep(1, fit$nobs)
#     n_Strata <- 1
#     if(isFALSE(start_stop)){
#       emax_times <- max(y[, "time"])  
#     } else{
#       emax_times <- max(y[, "stop"])
#     }
#   } else{
#     unique_strata <- unique(strata)
#     n_Strata <- length(unique_strata)
#     emax_times <- vector(mode = "numeric", length = n_Strata)
#     for(i in seq_along(unique_strata)){
#       if(isFALSE(start_stop)){
#         emax_times[i] <- max(y[which(strata == unique_strata[i]), "time"])  
#       } else{
#         emax_times[i] <- max(y[which(strata == unique_strata[i]), "stop"])
#       }
#     }
#   }
#   
#   
#   #Use baseHaz_cpp to quickly estimate (cum)hazard
#   Lambda0 <- baseHaz_cpp(
#     starttimes = start_times,
#     stoptimes = stop_times,
#     status = y[, "status"],
#     eXb = predict(object, newx = x, type = "response", s = s),
#     strata = strata, #strata not implemented yet
#     predtimes = times,
#     emaxtimes = emax_times,
#     nPatients = fit$nobs,
#     nStrata = n_Strata, #not implemented yet
#     cause = 1,
#     Efron = FALSE
#   )
#   #Linear predictor for new data
#   lp <- predict(object, newx = newdata, type = "response", s = s)
#   #Estimate absolute risk at specified times
#   #F(t) = 1- exp(-H0(t) * exp(beta * x))
#   out <- t(sapply(lp, function(x) 1- exp(-x * Lambda0$cumhazard)))
# }
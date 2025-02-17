L0 <- riskRegression::baseHaz_cpp(starttimes = info$start,
                                  stoptimes = info$stop,
                                  status = info$status,
                                  eXb = coxnet_pred,
                                  strata = 1,
                                  nPatients = NROW(info$stop),
                                  nStrata = 1,
                                  emaxtimes = max(info$stop),
                                  predtimes = sort(unique(info$stop)),
                                  cause = 1,
                                  Efron = TRUE)
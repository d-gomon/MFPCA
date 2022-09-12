create_x <- function(var_name, scale = TRUE){
  dattemp <- dat
  if(isTRUE(scale)){
    dattemp[, var_name] <- scale(dattemp[, var_name])
  }
  X.temp <- as.matrix(reshape(dattemp[, c("RID", "times", var_name)],
                              idvar = "RID", timevar = "times",
                              direction = "wide"))[, -1]
  X.temp <- X.temp[,order(unique(dattemp$times))]
  rownames(X.temp) <- unique(dattemp$RID)
  return(funData(argvals = sort(unique(dattemp$times)), X = X.temp))
}


long_to_mfdat <- function(dat, id_name, var_name, time_var, scale = TRUE){
  df <- dat
  if(isTRUE(scale)){
    df[, var_name] <- scale(df[, var_name])
  }
  X.temp <- as.matrix(reshape(df[, c(id_name, time_var, var_name)],
                              idvar = id_name, timevar = time_var,
                              direction = "wide"))[, -1]
  X.temp <- X.temp[,order(unique(df[, time_var]))]
  rownames(X.temp) <- unique(df[, id_name])
  return(funData(argvals = sort(unique(df[, time_var])), X = X.temp))
}


spaghetti <- function(var_name){
  var <- enquo(var_name)
  p <- ggplot(data = dat, aes(x = times, y = !!var, group = RID))
  p <- p + geom_line()
  p
}
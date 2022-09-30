#' @title Plot performance measures of models
#' @description Can plot output of cv and rcv
#' 
#' @describeIn plot Plot a CGR-CUSUM
#' @import RColorBrewer
#' @export

plot.rcv_mfpccox <- function(x, ...){
  plot(names(x$AUC_pred), x$AUC_pred)
}


plot_rcvlist <- function(x, ...){
  cols <- brewer.pal(length(x), "Dark2")
  ymax <- max(sapply(x, function(t) max(t$AUC_pred) ))
  ymin <- min(sapply(x, function(t) min(t$AUC_pred) ))
  xlims <- as.numeric(names(x[[1]]$AUC_pred))
  plot(NULL, ylim=c(ymin, ymax), ylab="tdAUC", xlab="Time since baseline (Years)",
       xlim = c(min(xlims), max(xlims)))
  for(i in 1:length(x)){
    lines(names(x[[i]]$AUC_pred), x[[i]]$AUC_pred, col = cols[i])
  }
}

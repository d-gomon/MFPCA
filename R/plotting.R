#' @title Plot performance measures of models
#' @description Can plot output of cv and rcv
#' 
#' @describeIn plot Plot a CGR-CUSUM
#' @import RColorBrewer
#' @export

plot.rcv_mfpccox <- function(x, ...){
  plot(names(x$AUC_pred), x$AUC_pred)
}


plot_rcvlist <- function(x, legend, ...){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  ymax <- max(sapply(x, function(t) max(t$AUC_pred) ))
  ymin <- min(sapply(x, function(t) min(t$AUC_pred) ))
  xlims <- as.numeric(names(x[[1]]$AUC_pred))
  plot(NULL, ylim=c(ymin, ymax), ylab="tdAUC", xlab="Time since baseline (Years)",
       xlim = c(min(xlims), max(xlims)), ...)
  lines(names(x[[1]]$AUC_pred), x[[1]]$AUC_pred, col = cols[1], lwd = 3, type = "b")
  for(i in 2:length(x)){
    lines(names(x[[i]]$AUC_pred), x[[i]]$AUC_pred, col = cols[i], lwd = 2)
  }
  legend("bottomleft", legend = legend, lty = 1, lwd = 2, col = cols[1:length(x)])
}

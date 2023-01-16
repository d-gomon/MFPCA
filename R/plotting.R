#' @title Plot performance measures of models
#' @description Can plot output of rcv_mfpccox() functions.
#' 
#' @param x A list containing "rcv_mfpccox" objects
#' @param legend A list of length(x) containing a description of each object in x.
#' @param legend_pos Legend position. See plot()
#' @param title Title for legend.
#' @param main_title Main title for plot.
#' @param col A vector of length(x) specifying which color to use for each object in x. Default uses qual pallette from brewer.pal.
#' @param ... Further parameters to plot()
#' 
#' @describeIn plot Plot prediction performance of rcv_mfpccox() output
#' @import RColorBrewer
#' @export


plot_rcvlist <- function(x, legend = "", legend_pos = NULL, title = NULL, main_title = "", col = NULL, ...){
  if(is.null(col)){
    #Automatically choose colors if none specified
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))  
  } else{
    col = cols
  }
  if("MSE" %in% names(x[[1]])){
    par(mfrow = c(1,3))
  } else{
    par(mfrow = c(1,2))
  }
  
  
  #First we plot AUC
  ymax <- max(sapply(x, function(t) max(t$AUC_pred, na.rm = TRUE) ))
  ymin <- min(sapply(x, function(t) min(t$AUC_pred, na.rm = TRUE) ))
  xlims <- as.numeric(names(x[[1]]$AUC_pred))
  plot(NULL, ylim=c(ymin, ymax), ylab="tdAUC (higher is better)", xlab="Time since baseline (Years)",
       xlim = c(min(xlims), max(xlims)), ...)
  lines(names(x[[1]]$AUC_pred), x[[1]]$AUC_pred, col = cols[1], lwd = 2, type = "l")
  if(length(x) >= 2){
    for(i in 2:length(x)){
      lines(names(x[[i]]$AUC_pred), x[[i]]$AUC_pred, col = cols[i], lwd = 2)
    }
  }
  
  
  #Then we plot Brier
  ymax <- max(sapply(x, function(t) max(t$Brier_pred, na.rm = TRUE) ))
  ymin <- min(sapply(x, function(t) min(t$Brier_pred, na.rm = TRUE) ))
  xlims <- as.numeric(names(x[[1]]$Brier_pred))
  plot(NULL, ylim=c(ymin, ymax), ylab="Brier Score (lower is better)", xlab="Time since baseline (Years)",
       xlim = c(min(xlims), max(xlims)), ...)
  lines(names(x[[1]]$Brier_pred), x[[1]]$Brier_pred, col = cols[1], lwd = 2, type = "l")
  if(length(x) >= 2){
    for(i in 2:length(x)){
      lines(names(x[[i]]$Brier_pred), x[[i]]$Brier_pred, col = cols[i], lwd = 2)
    }  
  }
  
  
  #Optionally we plot MSE
  if("MSE" %in% names(x[[1]])){
    ymax <- max(sapply(x, function(t) max(t$MSE, na.rm = TRUE) ))
    ymin <- min(sapply(x, function(t) min(t$MSE, na.rm = TRUE) ))
    xlims <- as.numeric(names(x[[1]]$MSE))
    plot(NULL, ylim=c(ymin, ymax), ylab="MSE (lower is better)", xlab="Time since baseline (Years)",
         xlim = c(min(xlims), max(xlims)), ...)
    lines(names(x[[1]]$MSE), x[[1]]$MSE, col = cols[1], lwd = 2, type = "l")
    if(length(x) >= 2){
      for(i in 2:length(x)){
        lines(names(x[[i]]$MSE), x[[i]]$MSE, col = cols[i], lwd = 2)
      }  
    }
  }
  
  
  
  if(is.null(legend_pos)){
    legend_pos = "bottomright"
  }
  mtext(main_title,                   # Add main title
        side = 3,
        line = - 2,
        outer = TRUE)
  legend(x = legend_pos, xpd = TRUE, inset = c(0, -0.1), legend = legend, lty = 1, lwd = 2, col = cols[1:length(x)], title = title)
  title(...)
}

#' @title Plot performance measures of models
#' @description Can plot output of rcv_mfpccox() functions.
#' 
#' @param x A list containing "rcv_mfpccox" objects
#' @param legend A list of length(x) containing a description of each object in x.
#' @param legend_pos Legend position. See plot()
#' @param title Title for legend.
#' @param main_title Main title for plot.
#' @param sub_title Subtitle for plot.
#' @param col A vector of length(x) specifying which color to use for each object in x. Default uses qual pallette from brewer.pal.
#' @param lty A vector of length(x)
#' @param ... Further parameters to plot()
#' 
#' @describeIn plot Plot prediction performance of rcv_mfpccox() output
#' @import RColorBrewer
#' @export


plot_rcvlist <- function(x, legend = "", legend_pos = NULL, title = NULL, main_title = "", sub_title = "", col = NULL, lty = NULL, ...){
  
  N <- length(x)
  print(N)
  
  #----------------Graphical Parameters-----------------
  if(is.null(col)){
    #Automatically choose colors if none specified
    #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    #cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))  
    cols <- brewer.pal(8, "Dark2")
  } else{
    cols = col
  }
  
  if(is.null(lty)){
    ltys <- rep(1, N)
  } else{
    ltys <- lty
  }
  
  
  
  if("MSE" %in% names(x[[1]])){
    par(mfrow = c(1,3))
    par(mar = c(5.1, 4.1, 7.1, 2.1))
  } else{
    par(mfrow = c(1,2))
    par(mar = c(5.1, 4.1, 7.1, 2.1))
  }
  
  
  #First we plot AUC
  ymax <- max(sapply(x, function(t) max(t$AUC_pred, na.rm = TRUE) ))
  ymin <- min(sapply(x, function(t) min(t$AUC_pred, na.rm = TRUE) ))
  xlims <- as.numeric(names(x[[1]]$AUC_pred))
  plot(NULL, ylim=c(ymin, ymax), ylab="tdAUC (higher is better)", xlab="Years since baseline",
       xlim = c(min(xlims), max(xlims)), cex.lab = 1.3, ...)
  lines(names(x[[1]]$AUC_pred), x[[1]]$AUC_pred, col = cols[1], lwd = 2, type = "l", lty = ltys[1])
  if(length(x) >= 2){
    for(i in 2:length(x)){
      lines(names(x[[i]]$AUC_pred), x[[i]]$AUC_pred, col = cols[i], lwd = 2, lty = ltys[i])
    }
  }
  
  
  #Then we plot Brier
  ymax <- max(sapply(x, function(t) max(t$Brier_pred, na.rm = TRUE) ))
  ymin <- min(sapply(x, function(t) min(t$Brier_pred, na.rm = TRUE) ))
  xlims <- as.numeric(names(x[[1]]$Brier_pred))
  plot(NULL, ylim=c(ymin, ymax), ylab="Brier Score (lower is better)", xlab="Years since baseline",
       xlim = c(min(xlims), max(xlims)), cex.lab = 1.3, ...)
  lines(names(x[[1]]$Brier_pred), x[[1]]$Brier_pred, col = cols[1], lwd = 2, type = "l", lty = ltys[1])
  if(length(x) >= 2){
    for(i in 2:length(x)){
      lines(names(x[[i]]$Brier_pred), x[[i]]$Brier_pred, col = cols[i], lwd = 2, lty = ltys[i])
    }  
  }
  
  
  #Optionally we plot MSE
  if("MSE" %in% names(x[[1]])){
    ymax <- max(sapply(x, function(t) max(t$MSE, na.rm = TRUE) ))
    ymin <- min(sapply(x, function(t) min(t$MSE, na.rm = TRUE) ))
    xlims <- as.numeric(names(x[[1]]$MSE))
    plot(NULL, ylim=c(ymin, ymax), ylab="MSE (lower is better)", xlab="Years since baseline",
         xlim = c(min(xlims), max(xlims)), cex.lab = 1.3, ...)
    lines(names(x[[1]]$MSE), x[[1]]$MSE, col = cols[1], lwd = 2, type = "l", lty = ltys[1])
    if(length(x) >= 2){
      for(i in 2:length(x)){
        lines(names(x[[i]]$MSE), x[[i]]$MSE, col = cols[i], lwd = 2, lty = ltys[i])
      }  
    }
  }
  
  
  
  if(is.null(legend_pos)){
    legend_pos = "bottomright"
  }
  mtext(main_title,                   # Add main title
        side = 3,
        line = - 3,
        outer = TRUE,
        cex = 1.3)
  mtext(sub_title,                   # Add sub title
        side = 3,
        line = - 5,
        outer = TRUE,
        cex = 1.2,
        col = "grey")
  legend(x = legend_pos, xpd = TRUE, legend = legend, lty = ltys, lwd = 2, col = cols[1:length(x)], title = title, inset = c(0, 1))
  title(...)
}

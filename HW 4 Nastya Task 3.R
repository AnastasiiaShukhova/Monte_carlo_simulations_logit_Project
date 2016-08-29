###  3 ROC STUFF
setwd("/Users/nastyashukhova/Dropbox/UniMannheim/PoliSci/SS_16/AQM/HW/HW4")
library(MASS)
library(ggplot2)
library(reshape)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
load("HW4Example.RData")
ROC  <- function(data, y, X){
  model  <- glm(y ~ X,data = data, family=binomial())
  covar <- cbind(1,X)
  # Predicted Probability
  mu <- covar %*% model$coef
  probability <- 1/(1 + exp(-mu))
  predictedValue <- NULL
  i = 1
  truePositiveRate  <- c()
  trueNegativeRate  <- c()
  for (threshold in seq(0, 1, by = 0.001)){
    predictedValue[probability >= threshold] <- 1
    predictedValue[probability < threshold] <- 0
    truePositiveRate[i] <- sum(y==1 & predictedValue==1, na.rm =T)/sum(y==1, na.rm =T) # True Negative Rate
    trueNegativeRate[i] <- sum(y==0 & predictedValue==0, na.rm =T)/sum(y==0, na.rm =T)
    i = i + 1
  }
  return(data.frame("threshold" = seq(0, 1, by = 0.001), truePositiveRate, trueNegativeRate))
}

rocResults  <-ROC(data = HW4, y =HW4$y, X = cbind(HW4$X1,HW4$X2,HW4$X3))

# plot
roc_Curve  <- function(rocResults){
  plot  <- ggplot(data = rocResults, aes(x = 1- trueNegativeRate, y=truePositiveRate)) + 
    geom_line() + geom_abline(intercept = 0, color = "blue")  + scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1))
  return(plot)
}
roc_Curve.task1  <- roc_Curve(rocResults)
ggsave("rocExample.pdf",roc_Curve.task1, width = 6, height = 6 )
# 2 Laitin Fearon vs Collier and Hoeffler

fearonRoc  <- ROC(data = fl.three, y = fl.three$onset, X = X.fl[,-1])
fearonPlot  <- roc_Curve(fearonRoc)
fearonPlot  <- fearonPlot +  labs(title = "ROC Curve for the Fearon and Laitin Model")

load("ch.RData")
ch  <- na.omit(ch)
X.Col  <- cbind(ch$sxp, ch$sxp2, ch$secm, ch$gy1, ch$peace, ch$geogia, ch$lnpop, ch$frac, ch$etdo4590)
ncol(X.Col)

CollierRoc  <- ROC(data = ch, y = ch$warsa, X = X.Col)
CollierPlot  <- roc_Curve(CollierRoc)
CollierPlot  <- CollierPlot + labs(title = "ROC Curve for the Collier Hoeffler Model")

pdf("rocComparison.pdf", width = 15)
multiplot(fearonPlot, CollierPlot, cols = 2)
dev.off()

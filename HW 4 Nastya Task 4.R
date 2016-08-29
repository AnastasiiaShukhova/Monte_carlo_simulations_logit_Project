#### HW 4
setwd("/Users/nastyashukhova/Dropbox/UniMannheim/PoliSci/SS_16/AQM/HW/HW4")
library(MASS)
library(ggplot2)
library(reshape)


ll.logit <- function(theta, y, X) {
  # theta consists merely of beta (dim is ncol(X))
  beta <- theta[1:ncol(X)]
  # linear predictor; make sure that X is stored as.matrix
  mu <- X %*% beta
  # link function
  p <- 1/(1 + exp(-mu))
  # log-likelihood
  ll <- y * log(p) + (1 - y) * log(1 - p)
  # sum
  ll <- sum(ll)
  return(ll)
}

Convergedlogit  <- function(y, X, startMethod ="Nelder-Mead", init = 0){
  inits <- c(rep(init,ncol(X)))
  logitResults  <- optim(inits,ll.logit,control=list(fnscale=-1),y=y,X=X,hessian=F
                         ,method="Nelder-Mead")
  if (logitResults$convergence == 0){print("Algorithm converged")
                                     return(logitResults)}
  if (logitResults$convergence == 1){
    print("Fucking shit. Algorithm didn't converge")
    logitResults2  <- optim(logitResults$par,ll.logit,control=list(fnscale=-1),y=y,X=X, method="Nelder-Mead")
    if (logitResults2$convergence == 0){print("Converged")
                                        return(logitResults2)}
    
    if (logitResults2$convergence == 1){
      logitResults3  <- optim(logitResults2$par,ll.logit,control=list(fnscale=-1),y=y,X=X,method="BFGS")
      
      if (logitResults3$convergence == 0){print("Converged")  
                                          return(logitResults3)}
      if (logitResults3$convergence == 1){
        logitResults4  <- optim(logitResults3$par,ll.logit,control=list(fnscale=-1),y=y,X=X,method="BFGS")
        if (logitResults4$convergence == 0){print("Converged") 
                                            return(logitResults4)} 
        if (logitResults4$convergence == 1){print("I give up. algorithm didn't converge")} 
      }
    }
  }
}

ols <- function(X,y){
  beta  <- solve(t(X)%*%X)%*%t(X)%*%y #calculate beta
  e  <- y - X%*%beta #residuals
  n = length(y) 
  k = ncol(X)
  sigma.sq <- (t(e)%*%e)/(n-k)
  var  <- as.numeric(sigma.sq)*solve(t(X)%*%X)
  std.er  <- sqrt(diag(var))
  return(cbind(beta, std.er))
}


# TASK 4

## 4 Specification Issues 
set.seed(123)
n <- 1000000
X1 <- rnorm(n,0,1)
X2 <- rnorm(n,0,3)
beta0 <- 0.5
beta1 <- 1
beta2 <- 1.5
mu <- beta0 + beta1*X1 + beta2*X2
p <- 1/(1 + exp(-mu))
y <- rbinom(n,1,p)
# 1 
cor(X1, X2) #uncorrelated
X  <- cbind(1, X1, X2)
# 2
full.logit  <- Convergedlogit(y =y, X=cbind(1, X1, X2))
restricted.logit  <- Convergedlogit(y =y, X=cbind(1, X1))

full.ols  <- ols(y =y, X=cbind(1, X1, X2))
restricted.ols  <- ols(y =y, X=cbind(1, X1))
library(xtable)

results <- cbind(round(full.logit$par, 4), c(round(restricted.logit$par, 4), "NA"), round(full.ols[,1], 4), c(round(restricted.ols[,1], 4), "NA"))
rownames(results) <- c ("Constant","X1","X2")
colnames(results) <- c("Full Logit","Restricted Logit","Full OLS","Restricted OLS")
xtable(results)

# 3
# Logit

monte.logit(X=X, y=y)

monteCarloLogit <- function(y, X, startMethod ="Nelder-Mead", init = 0){
  id.sample <- sample(length(y),size=500)
  reg <- Convergedlogit(y=y[id.sample], X=X[id.sample, ],startMethod = "Nelder-Mead", init = init ) 
  b1.hat <- reg$par # get coefficient and standard error
  return(b1.hat)
}
monteLogitFull  <- replicate(1000, monteCarloLogit( y=y, X=X))
monteLogitRestricted  <- replicate(1000, monteCarloLogit( y=y, X=X[,1:2]))
mean(monteLogitFull[2,])
mean(monteLogitRestricted[2,])
#plot beta zeros LOGIT
plotBetas  <- function(mcDataFull, mcDataRestricted, beta = 0, model = "logit"){
  if (model == "logit"){
  if (beta == 0){
    betas0  <- melt(data.frame("1" = mcDataFull[1,], "2"=mcDataRestricted[1,]))
    meansb0 <- melt(data.frame("1"=mean(betas0[1:1000,2]), "2"=mean(betas0[1001:2000,2])))
    betaZeroPlot  <- ggplot(data =  NULL, aes(x = betas0[,2], colour = betas0[,1])) +
      geom_density() +geom_vline(data=NULL, aes(xintercept=meansb0[,2], colour = meansb0[,1]),
                                 linetype="dashed", size=1) +
      labs(title = "Plot of MC results for restricted and \n full logit models. Beta 0", x = "Beta Zero") + 
      scale_colour_discrete(name="Models",
                          labels=c("Full", "Restricted")) + theme_bw()
                          
    return(betaZeroPlot)
  }
  if (beta == 1){
    betas1  <- data.frame("1" = mcDataFull[2,], "2"=mcDataRestricted[2,])
    betas1  <- melt(betas1)
    meansb1 <- melt(data.frame("1"=mean(betas1[1:1000,2]), "2"=mean(betas1[1001:2000,2])))
    beta1Plot  <- ggplot(data =  NULL, aes(x = betas1[,2], colour = betas1[,1])) +
      geom_density() + geom_vline(data=NULL, aes(xintercept=meansb1[,2], colour = meansb1[,1]),
                                 linetype="dashed", size=1) + 
      labs(title = "Plot of MC results for restricted and \n full logit models. Beta 1", x= "Beta One") + 
      scale_colour_discrete(name="Models",
                            labels=c("Full", "Restricted")) + theme_bw()
    return(beta1Plot)
  }
  }
  if (model == "ols"){
    if (beta == 0){
      betas0  <- melt(data.frame("1" = mcDataFull[1,], "2"=mcDataRestricted[1,]))
      meansb0 <- melt(data.frame("1"=mean(betas0[1:1000,2]), "2"=mean(betas0[1001:2000,2])))
      betaZeroPlot  <- ggplot(data =  NULL, aes(x = betas0[,2], colour = betas0[,1])) +
        geom_density() +geom_vline(data=NULL, aes(xintercept=meansb0[,2], colour = meansb0[,1]),
                                   linetype="dashed", size=1) +
        labs(title = "Plot of MC results for restricted and \n full OLS models. Beta 0", x = "Beta Zero") + 
        scale_colour_discrete(name="Models",
                              labels=c("Full", "Restricted")) + theme_bw()
      
      return(betaZeroPlot)
    }
    if (beta == 2){
      betas1  <- data.frame("1" = mcDataFull[2,], "2"=mcDataRestricted[2,])
      betas1  <- melt(betas1)
      meansb1 <- melt(data.frame("1"=mean(betas1[1:1000,2]), "2"=mean(betas1[1001:2000,2])))
      beta1Plot  <- ggplot(data =  NULL, aes(x = betas1[,2], colour = betas1[,1])) +
        geom_density() + geom_vline(data=NULL, aes(xintercept=meansb1[,2], colour = meansb1[,1]),
                                    linetype="dashed", size=1) + 
        labs(title = "Plot of MC results for restricted and \n full OLS models. Beta 1", x= "Beta One") + 
        scale_colour_discrete(name="Models",
                              labels=c("Full", "Restricted")) + theme_bw()
      return(beta1Plot)
    }
  }
}

plotLogitB0  <- plotBetas(monteLogitFull, monteLogitRestricted, beta = 0, model = "logit"  )
plotLogitB1  <- plotBetas(monteLogitFull, monteLogitRestricted, beta = 1 , model = "logit")
#plot beta ones OLS
# ols
monteCarloOls <- function(y, X){
  id.sample <- sample(length(y),size=500)
  reg <- ols(y=y[id.sample], X=X[id.sample,]) 
  b1.hat <- reg[,1] # get coefficient and standard error
  return(b1.hat)
}
monteOlsFull <-  replicate(1000, monteCarloOls( y=y, X=X))
monteOlsRestricted <- replicate(1000, monteCarloOls( y=y, X=X[,1:2]))


plotOLSB0  <- plotBetas(monteOlsFull, monteOlsRestricted, beta = 0, model = "ols" )
plotOLSB1  <- plotBetas(monteOlsFull, monteOlsRestricted, beta = 1, , model = "ols" )
mean(monteLogitFull[1,])
mean(monteLogitRestricted[1,])
#plot results
pdf("MonteCarloPlots1.pdf", width = 13)
multiplot(plotLogitB0, plotOLSB0,plotLogitB1, plotOLSB1, cols=2)
dev.off()
# 6 
#bias should increase 
# Increase variance and mean
set.seed(123)
X1 <- rnorm(n,0,1)
X2 <- rnorm(n,2,6)
mu <- beta0 + beta1*X1 + beta2*X2
p <- 1/(1 + exp(-mu))
y <- rbinom(n,1,p)

# 1 
cor(X1, X2) #uncorrelated
X  <- cbind(1, X1, X2)
# 2
monteLogitFull  <- replicate(1000, monteCarloLogit(y=y, X=X))
monteLogitRestricted  <- replicate(1000, monteCarloLogit(y=y, X=X[,1:2]))

plotLogitB0IncreasedVarMean  <- plotBetas(monteLogitFull, monteLogitRestricted, beta = 0 ) + labs(title = "Plot of MC results for restricted and \n full logit models. Beta 0. \n Variance = 6, mean =2 ")
plotLogitB1IncreasedVarMean  <- plotBetas(monteLogitFull, monteLogitRestricted, beta = 1 )+ labs(title = "Plot of MC results for restricted and \n full logit models. Beta 1. \n Variance = 6, mean =2 ")
mean(monteLogitFull[2,])-mean(monteLogitRestricted[2,])
mean(monteLogitFull[1,])-mean(monteLogitRestricted[1,])

set.seed(123)
X2 <- rnorm(n,3,7)
mu <- beta0 + beta1*X1 + beta2*X2
p <- 1/(1 + exp(-mu))
y <- rbinom(n,1,p)
# 1 
cor(X1, X2) #uncorrelated
X  <- cbind(1, X1, X2)

# 2
# res  <- glm(y~X, data =NULL, family=binomial)
monteLogitFull  <- replicate(1000, monteCarloLogit(y=y, X=X))
monteLogitRestricted  <- replicate(1000, monteCarloLogit( y=y, X=X[,1:2]))

plotLogitB0IncreasedVarMean2  <- plotBetas(monteLogitFull, monteLogitRestricted, beta = 0 ) + labs(title = "Plot of MC results for restricted and \n full logit models. Beta 0. \n Variance = 8, mean =4 ")
plotLogitB1IncreasedVarMean2  <- plotBetas(monteLogitFull, monteLogitRestricted, beta = 1 )+ labs(title = "Plot of MC results for restricted and \n full logit models. Beta 1. \n Variance = 8, mean =4 ")



#plot results
pdf("MonteCarloPlots2.pdf", width = 13)
multiplot(plotLogitB0IncreasedVarMean, plotLogitB0IncreasedVarMean2,plotLogitB1IncreasedVarMean, plotLogitB1IncreasedVarMean2, cols=2)
dev.off()

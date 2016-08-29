#### HW 4
setwd("/Users/nastyashukhova/Dropbox/UniMannheim/PoliSci/SS_16/AQM/HW/HW4")
library(MASS)
library(ggplot2)
library(reshape)

# 1 logit, Probit, and the Linear Probability Model
load("fl.three.RData")

y.fl  <-  fl.three$onset
X.fl   <- cbind(1, fl.three$warl, fl.three$gdpenl, fl.three$lpopl1, fl.three$lmtnest, fl.three$ncontig,
                fl.three$Oil, fl.three$nwstate, fl.three$instab, fl.three$polity2l, fl.three$ethfrac, 
                fl.three$relfrac)

# LOGIT 

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
  logitResults  <- optim(inits,ll.logit,control=list(fnscale=-1),y=y,X=X,hessian=T
                          ,method="Nelder-Mead")
  if (logitResults$convergence == 0){print("Algorithm converged")
                                      return(logitResults)}
  if (logitResults$convergence == 1){
    print("Fucking shit. Algorithm didn't converge")
    logitResults2  <- optim(logitResults$par,ll.logit,control=list(fnscale=-1),y=y,X=X, hessian=T, method="Nelder-Mead")
    if (logitResults2$convergence == 0){print("Converged")
                                         return(logitResults2)}
    
    if (logitResults2$convergence == 1){
      logitResults3  <- optim(logitResults2$par,ll.logit,control=list(fnscale=-1),y=y,X=X,hessian=T, method="BFGS")
      
      if (logitResults3$convergence == 0){print("Converged")  
                                           return(logitResults3)}
      if (logitResults3$convergence == 1){
        logitResults4  <- optim(logitResults3$par,ll.logit,control=list(fnscale=-1),y=y,X=X,hessian=T, method="BFGS")
        if (logitResults4$convergence == 0){print("Converged") 
                                             return(logitResults4)} 
        if (logitResults4$convergence == 1){print("I give up. algorithm didn't converge")} 
      }
    }
  }
}
logitResults <-  Convergedlogit(X=X.fl, y=y.fl)
logitResults


#PROBIT
ll.probit <- function(theta, y, X){
  # theta consists merely of beta (dim is ncol(X))
  beta <- theta[1:ncol(X)]
  # linear predictor; make sure that X is stored as.matrix
  mu <- as.matrix(X) %*% beta
  # link function
  p <- pnorm(mu)
  # log-likelihood
  ll <- y * log(p) + (1 - y) * log(1 - p)
  # sum
  ll <- sum(ll)
  return(ll)
}
inits <- c(rep(0,ncol(X.fl)))
probitResults  <- optim(inits,ll.probit,control=list(fnscale=-1),y=y.fl,X=X.fl,hessian=TRUE
                        ,method="Nelder-Mead")
ConvergedProbit  <- function(y, X, startMethod ="Nelder-Mead", init = 0){
  inits <- c(rep(init,ncol(X)))
  probitResults  <- optim(inits,ll.probit,control=list(fnscale=-1),y=y,X=X,hessian=T
                         ,method="Nelder-Mead")
  if (probitResults$convergence == 0){print("Algorithm converged")
                                     return(probitResults)}
  if (probitResults$convergence == 1){
    print("Fucking shit. Algorithm didn't converge")
    probitResults2  <- optim(probitResults$par,ll.probit,control=list(fnscale=-1),y=y,X=X, hessian=T, method="Nelder-Mead")
    if (probitResults2$convergence == 0){print("Converged")
                                        return(probitResults2)}
    
    if (probitResults2$convergence == 1){
      probitResults3  <- optim(probitResults2$par,ll.probit,control=list(fnscale=-1),y=y,X=X,hessian=T, method="BFGS")
      
      if (probitResults3$convergence == 0){print("Converged")  
                                            return(probitResults3)}
      if (probitResults3$convergence == 1){
        probitResults4  <- optim(probitResults3$par,ll.probit,control=list(fnscale=-1),y=y,X=X,hessian=T,method="BFGS")
        if (probitResults4$convergence == 0){print("Converged") 
                                                    return(probitResults4)} 
        if (probitResults4$convergence == 1){print("I give up. algorithm didn't converge")} 
      }
    }
    }
  }
probitResults <- ConvergedProbit(y =y.fl, X = X.fl)
probitResults$par


### OLS 
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

olsResults  <- ols(X=X.fl, y=y.fl )

# 2 

#### logit ####
logit.coef <- logitResults$par
logit.variance <-  solve(-logitResults$hessian)
logit.se <- sqrt(diag(solve(-logitResults$hessian)))
sim.coef.logit <- mvrnorm(1000,logit.coef,logit.variance)
dim(sim.coef.logit)
## B. Calculate Expected Values
gdp.sim <- seq(min(fl.three$gdpenl), max(fl.three$gdpenl), length.out =100)
#no prior war, non-contiguous state, no oil, not new state, stability (instab = 0)
scenario <- cbind(1,0, gdp.sim, mean(fl.three$lpopl1), mean(fl.three$lmtnest), 1,
                  0, 0, 0, mean(fl.three$polity2l), mean(fl.three$ethfrac), 
                  mean(fl.three$relfrac))

logit.Xbeta <- sim.coef.logit %*% t(scenario)

# To get expected values for p, we need to plug in the Xbeta values into
# the link function to get simulatd probabilities
logitProbsSim <- (exp(logit.Xbeta))/(1 + exp(logit.Xbeta))
# Means and Quantiles
logitProbsMean <- apply(logitProbsSim,2,mean)
logitProbsQuantiles <- t(apply(logitProbsSim,2,quantile,prob=c(0.025,0.975)))
##Logit Plot
logitPlot  <- ggplot(aes(x=gdp.sim,y =logitProbsMean)) + geom_line(color = "black") +
  geom_smooth(aes(ymin=logitProbsQuantiles[,1], ymax = logitProbsQuantiles[2,]), stat="identity")+
  theme_bw()
logitPlot


#### PROBIT ####
probit.coef <- probitResults$par
probit.variance <-  solve(-probitResults$hessian)
probit.se <- sqrt(diag(solve(-probitResults$hessian)))
sim.coef.probit <- mvrnorm(1000,probit.coef,probit.variance)
dim(sim.coef.probit)
##  Calculate Expected Values
gdp.sim <- seq(min(fl.three$gdpenl), max(fl.three$gdpenl), length.out =100)
#no prior war, non-contiguous state, no oil, not new state, stability (instab = 0)
scenario <- cbind(1,0, gdp.sim, mean(fl.three$lpopl1), mean(fl.three$lmtnest), 1,
                  0, 0, 0, mean(fl.three$polity2l), mean(fl.three$ethfrac), 
                  mean(fl.three$relfrac))

probit.Xbeta <- sim.coef.probit %*% t(scenario)

probitProbsSim <- pnorm(probit.Xbeta)

# Means and Quantiles
probitProbsMean <- apply(probitProbsSim,2,mean)
probitProbsQuantiles <- t(apply(probitProbsSim,2,quantile,prob=c(0.025,0.975)))
dim(logitProbsSim)

##Logit Plot
probitPlot  <- ggplot(aes(x=gdp.sim,y =probitProbsMean), data = NULL) + geom_line(color = "black") +
  geom_smooth(aes(ymin=probitProbsQuantiles[,1], ymax = probitProbsQuantiles[,2]), stat="identity")+
  theme_bw()
probitPlot

# LM 
varcov  <- function(X,y){
  beta  <- ols(X, y)[,1] #calculate beta
  e  <- y - X%*%beta #residuals
  n = length(y) 
  k = ncol(X)
  sigma.sq <- (t(e)%*%e)/(n-k)
  var  <- as.numeric(sigma.sq)*solve(t(X)%*%X)
  return(var)
}

#### Linear Model ####
olsResults  <- ols(X=X.fl, y=y.fl )
lm.coef <- olsResults[,1]
lm.variance <-  varcov(X=X.fl, y=y.fl)
sim.coef.lm <- mvrnorm(1000,lm.coef,lm.variance)
dim(sim.coef.lm)
## B. Calculate Expected Values
gdp.sim <- seq(min(fl.three$gdpenl), max(fl.three$gdpenl), length.out =100)
#no prior war, non-contiguous state, no oil, not new state, stability (instab = 0)
scenario <- cbind(1,0, gdp.sim, mean(fl.three$lpopl1), mean(fl.three$lmtnest), 1,
                  0, 0, 0, mean(fl.three$polity2l), mean(fl.three$ethfrac), 
                  mean(fl.three$relfrac))
dim(scenario)
lmProbsSim <- sim.coef.lm %*% t(scenario)

# Means and Quantiles
lmProbsMean <- apply(lmProbsSim,2,mean)
lmProbsQuantiles <- t(apply(lmProbsSim,2,quantile,prob=c(0.025,0.975)))
dim(logitProbsSim)
##Logit Plot
lmPlot  <- ggplot(aes(x=gdp.sim,y =lmProbsMean), data = NULL) + geom_line(color = "black") +
  geom_smooth(aes(ymin=lmProbsQuantiles[,1], ymax = lmProbsQuantiles[,2]), stat="identity")+
  theme_bw()
lmPlot


####### 2 An Alternative Link Function 
# 1
# linear model 
# curve()
# logit 
curve(1/(1 + exp(-x)), -5, 5)
#probit
curve(pnorm, -5, 5)
#log-log 
curve(1-exp(-exp(x)), -5, 5)

curvePlot  <- ggplot(data = NULL, aes(x = seq(-5,5, length.out =100), y = 1/(1 + exp(-seq(-5,5, length.out =100))))) + geom_line()
curvePlot  <- curvePlot +geom_line(data = NULL, aes(x = seq(-5,5, length.out =100), y =pnorm(seq(-5,5, length.out =100))), color = "blue")
curvePlot  <- curvePlot + 
  geom_line(data = NULL, aes(x = seq(-5,5, length.out =100), y =1-exp(-exp(seq(-5,5, length.out =100)))), color = "red")


curve(pnorm,-5,5,col = "red", main = "Plot of different link functions" )#probit
curve(1/(1 + exp(-x)),-5, 5, col = "blue", add = TRUE)#logit
curve(1- exp(-exp(x)), -5,5, col = "green", add = TRUE)#cloglog
#we omitted the linear model as you wrote in the email

#add legend
legend(-5, 0.6, legend=c("Probit", "Logit", "Cloglog"),
       col=c("red", "blue", "green"), lty=1:2, cex=0.4)
?legend
# 2 log-log
ll.loglog <- function(theta, y, X) {
  # theta consists merely of beta (dim is ncol(X))
  beta <- theta[1:ncol(X)]
  # linear predictor; make sure that X is stored as.matrix
  mu <- X %*% beta
  # link function
  p <- 1-exp(-exp(mu))
  # log-likelihood
  ll <- sum(y * log(p) + (1 - y) * log(1 - p))
  return(ll)
}
logLog  <- function(X, y){
  inits <- c(rep(0, ncol(X)))
  logitResults  <- optim(inits
                         ,ll.loglog
                         ,control=list(fnscale=-1)
                         ,y=y,X=X
                         ,hessian=TRUE
                         ,method="CG"
  ) 
  return(logitResults)
}
logLogResults <- logLog(X=X.fl, y = y.fl)

loglog.coef <- logLogResults$par
loglog.variance <-  solve(-logLogResults$hessian)
loglog.se  <- sqrt(diag(solve(-logLogResults$hessian)))
sim.coef.loglog <- mvrnorm(1000,loglog.coef,loglog.variance)



## B. Calculate Expected Values
gdp.sim <- seq(min(fl.three$gdpenl), max(fl.three$gdpenl), length.out =100)
#no prior war, non-contiguous state, no oil, not new state, stability (instab = 0)
scenario <- cbind(1,0, gdp.sim, mean(fl.three$lpopl1), mean(fl.three$lmtnest), 1,
                  0, 0, 0, mean(fl.three$polity2l), mean(fl.three$ethfrac), 
                  mean(fl.three$relfrac))

loglog.Xbeta <- sim.coef.loglog %*% t(scenario)

# To get expected values for p, we need to plug in the Xbeta values into
# the link function to get simulatd probabilities
loglogProbsSim <- 1- exp(-exp(loglog.Xbeta))
# Means and Quantiles
loglogProbsMean <- apply(loglogProbsSim,2,mean)
loglogProbsQuantiles <- t(apply(loglogProbsSim,2,quantile,prob=c(0.025,0.975)))
dim(loglogProbsSim)

##loglog Plot
loglogPlot  <- ggplot(aes(x=gdp.sim,y =loglogProbsMean), data = NULL) + geom_line(color = "black") +
  geom_smooth(aes(ymin=loglogProbsQuantiles[,1], ymax = loglogProbsQuantiles[,2]), stat="identity")+
  theme_bw()
loglogPlot

summary(m1 <- glm(onset ~ warl + gdpenl + lpopl1 + lmtnest + ncontig + Oil + nwstate + instab + polity2l + ethfrac + relfrac, 
                  data = fl.three, family=binomial(link="logit")))
summary(m2 <- glm(onset ~ warl + gdpenl + lpopl1 + lmtnest + ncontig + Oil + nwstate + instab + polity2l + ethfrac + relfrac, 
                  data = fl.three, family=binomial(link="probit")))
summary(m3 <- lm(onset ~ warl + gdpenl + lpopl1 + lmtnest + ncontig + Oil + nwstate + instab + polity2l + ethfrac + relfrac, 
                 data = fl.three))
summary(m_c.ll <- glm(onset ~ warl + gdpenl + lpopl1 + lmtnest + ncontig + Oil + nwstate + instab + polity2l + ethfrac + relfrac, 
                      data = fl.three, family=binomial(link="cloglog")))
suppressMessages(library(stargazer))
stargazer(m1,m2,m2,m_c.ll, font.size = "footnotesize")



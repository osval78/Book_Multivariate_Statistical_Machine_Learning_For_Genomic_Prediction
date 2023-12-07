rm(list=ls(all=TRUE))
library(mvtnorm)
library(MASS)
source('TR_RR.R')
set.seed(10)
n =  100
p = 500
Var =  0.25*diag(p)+0.75
Tab =  data.frame()
betav =  rnorm(p,rpois(p,1),2)
plot(betav,xlab=expression(j),ylab=expression(beta[j]))
for(i in 1:100)
{
  X =  rmvnorm(n,rep(0,p),Var)
  dim(X)
  y = 5+X%*%betav + rnorm(n,0,0.5)
  Pos_tr = sample(1:n,n*0.80)
  y_tr = y[Pos_tr]; X_tr =  X[Pos_tr,]
  y_tst = y[-Pos_tr]; X_tst =  X[-Pos_tr,]
  #Training RR
  TR_RR = Tr_RR_f(y_tr,X_tr,K=5,KG=100,KR=1)
  TR_RR$lamb_o
  lambv =  TR_RR$lambv
  plot(log(TR_RR$lambv),TR_RR$iPEv_mean)
  #Prediction RR in testing data
  Pred_RR = Pred_RR_f(y_tst,X_tst,TR_RR)
  Pred_RR$MSEP
  
  #OLS
  Pred_ols = Pred_ols_f(y_tst,X_tst,y_tr,X_tr)
  Pred_ols
  
  Tab =  rbind(Tab,data.frame(Sim=i,MSEP_RR = Pred_RR$MSEP,
                              MSEP_ols = Pred_ols$MSEP))
  cat('i = ', i,'\n')}

Mean_v = colMeans(Tab)
(Mean_v[3]-Mean_v[2])/Mean_v[2]*100

mean(Tab$MSEP_RR<Tab$MSEP_ols)
Pos =  which(Tab$MSEP_RR<Tab$MSEP_ols)
mean((Tab$MSEP_ols[Pos]-Tab$MSEP_RR[Pos])/Tab$MSEP_RR[Pos])*100
mean((Tab$MSEP_RR[-Pos]-Tab$MSEP_ols[-Pos])/Tab$MSEP_ols[-Pos])*100

plot(Tab$MSEP_RR, Tab$MSEP_ols,
     col=ifelse(Tab$MSEP_RR<Tab$MSEP_ols,3,2),
     xlab='MSE Ridge Regression', ylab='MSE OLS')
abline(a=0,b=1)

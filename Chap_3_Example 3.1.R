
#################R code for example 3.1 #################################
#Syntetic data: 100 observation and 2 regressors with Gradient descendent method
rm(list=ls())
library(mvtnorm)
set.seed(1)
X =  cbind(1,rmvnorm(100,c(0,0),diag(2)))
betav = c(5,1,2.1)
y = X%*%betav+rnorm(100,0,1.2)
dat = data.frame(y=y,x1 = X[,2],x2=X[,3])
plot(dat)

alpha = 1e-2
tol      = 1e-8 
p =  2
betav_0 =  c(mean(y),rep(0,p))
tol.e = 1
Iter = 0
tX =  t(X)
while(tol<tol.e)
{
  Iter =  Iter + 1
  e = y-X%*%betav_0
  betav_t = betav_0 + alpha*tX%*%e
  tol.e = max(abs(betav_t-betav_0))
  betav_0 =  betav_t
}
betav_t
tol.e
Iter

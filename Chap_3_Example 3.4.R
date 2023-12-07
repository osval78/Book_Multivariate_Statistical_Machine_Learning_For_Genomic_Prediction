rm(list=ls())
library(BMTME)
data("WheatMadaToy")
dat_F = phenoMada
dim(dat_F)
dat_F$GID = as.character(dat_F$GID) 
G =  genoMada
eig_G =  eigen(G)
G_0.5 = eig_G$vectors%*%diag(sqrt(eig_G$values))%*%t(eig_G$vectors)
X =  G_0.5
y = dat_F$PH
n =  length(y)
library(glmnet)
#5FCV
set.seed(3)
K = 5
Tab =  data.frame()
set.seed(1) 
for(k in 1:100)
{
  Pos_tr = sample(1:n,n*0.8)
  y_tr = y[Pos_tr] ; X_tr =  X[Pos_tr,]; n_tr =  dim(X_tr)[1]
  y_tst = y[-Pos_tr]; X_tst = X[-Pos_tr,]
  #Partition for internal training the model
  Grpv_k = findInterval(cut(sample(1:n_tr,n_tr),breaks=5),1:n_tr)
  #RR
  A_RR = cv.glmnet(X_tr,y_tr,alpha=0,foldid=Grpv_k,type.measure='mse')
  yp_RR = predict(A_RR,newx=X_tst,s='lambda.min')
  #LR
  A_LR = cv.glmnet(X_tr,y_tr,alpha=1,foldid=Grpv_k,type.measure='mse')
  yp_LR = predict(A_LR,newx=X_tst,s='lambda.min')
  #OLS
  A_OLS = glmnet(X_tr,y_tr,alpha=1,lambda=0)
  yp_OLS = predict(A_OLS,newx=X_tst)
  
  
  Tab =  rbind(Tab,data.frame(PT=k,MSEP_RR =  mean((y_tst-yp_RR)^2),
                              MSEP_LR =  mean((y_tst-yp_LR)^2),
                              MSEP_OLS = mean((y_tst-yp_OLS)^2)))
  cat('k = ', k,'\n')
}
Tab

Mean_v = colMeans(Tab)
(Mean_v[3]-Mean_v[2])/Mean_v[2]*100
Mean_v

plot(Tab$MSEP_RR, Tab$MSEP_LR,ylim=c(0,1200), xlim=c(0,300),
     col=ifelse(Tab$MSEP_RR<Tab$MSEP_LR,3,2),
     xlab='MSE Ridge Regression', ylab='MSE Lasso Regression')
abline(a=0,b=1)


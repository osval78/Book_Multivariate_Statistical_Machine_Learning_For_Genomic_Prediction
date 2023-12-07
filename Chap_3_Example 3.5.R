#############R code for example 3.5###################################
load(file ='dat-E3.5.RData')
dat_F = dat$dat_F
dat_F =  dat_F[order(dat_F$Rep,dat_F$GID),]
head(dat_F)
G = dat$G
dat_F$y = dat_F$Height
ZL =  model.matrix(~0+GID,data=dat_F)
colnames(ZL)
Pos =  match(colnames(ZL),paste('GID',colnames(G),sep=''))
max(abs(diff(Pos)))
y =  dat_F$y
ei =  eigen(G)
X  = ZL%*%ei$vectors%*%diag(sqrt(ei$values))%*%t(ei$vectors)
n =  length(y)
library(glmnet)
#5FCV
set.seed(1)
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
  A_RR = cv.glmnet(X_tr,y_tr,family='binomial',
                   alpha=0,foldid=Grpv_k,type.measure='class')
  yp_RR = as.numeric(predict(A_RR,newx=X_tst,s='lambda.min',type='class'))
  #LR
  A_LR = cv.glmnet(X_tr,y_tr,family='binomial',
                   alpha=1,foldid=Grpv_k,type.measure='class')
  yp_LR = as.numeric(predict(A_LR,newx=X_tst,s='lambda.min',type='class'))
  #SLR
  A_SLR = glmnet(X_tr,y_tr,family='binomial',alpha=0,lambda=0)
  yp_SLR = as.numeric(predict(A_SLR,newx=X_tst,type='class'))
  
  
  Tab =  rbind(Tab,data.frame(PT=k,PCCC_RR =  1-mean(y_tst!=yp_RR),
                              PCCC_LR =  1-mean(y_tst!=yp_LR),
                              PCCC_SLR = 1-mean(y_tst!=yp_SLR)))
  cat('k = ', k,'\n')
}
Tab

Tr_RR_f<-function(y_tr,X_tr,K=5,KG=100,KR=1)
{
  n_tr =  dim(X_tr)[1]
  X_tr_s = scale2_f(X_tr)
  mu      = mean(y_tr)
  y_tr  = y_tr-mu
  #Inner CV
  lambv = lamb_f(X_tr_s,K=KG,li=1e-7,ls=1-1e-7)
  iPEv_mean = 0 
  for(ir in 1:KR)
  {
    iPE_mat = matrix(0,nr=length(lambv),nc=K)
    MSEP_mat=iPE_mat
    Grpv = findInterval(cut(sample(1:n_tr,n_tr),breaks=K),1:n_tr) 
    for(i in 1:K)
    {
      Pos_itr = which(Grpv!=i)
      X_itr  = X_tr[Pos_itr,]
      y_itr   = y_tr[Pos_itr]; 
      y_itst   = y_tr[-Pos_itr]; 
      dat_itr = data.frame(y=y_itr,X=X_itr)
      iPEv  = lambv
      n_itr =  dim(X_itr)[1]
      betav_itr = A_RR_f_a_V(y_itr,X_itr,lambv)
      yp_mat     = 0  + scale2_f(X_tr[-Pos_itr,])%*%betav_itr
      iPEv       = colMeans((matrix(y_itst,nr=length(y_itst),
                                    nc=length(lambv))-yp_mat)**2)
      iPE_mat[,i] =iPEv
    }
    iPEv_mean = (ir-1)/ir*iPEv_mean + rowMeans(iPE_mat)/ir
  }  
  
  plot(log(lambv),iPEv_mean)#,xlim=c(3,5))
  Pos_o = which.min(iPEv_mean)
  lamb_o = lambv[Pos_o]
  A_RR  = RR_f(y=y_tr+mu,X=X_tr,lamb=lamb_o,Intercept = TRUE,Std=TRUE)
  betav_s = A_RR$betav_s
  betav_o = A_RR$betav_o
  list(lamb_o = lamb_o, betav_s =  betav_s,
       betav_o = betav_o,
       lambv=lambv,iPEv_mean=iPEv_mean,X_tr=X_tr,y_tr =y_tr,Grpv=Grpv)
}

Pred_RR_f<-function(y_tst,X_tst,TR_RR)
{
  #betava_RR = TR_RR$betava_RR
  y_tr = TR_RR$y_tr; X_tr = TR_RR$X_tr
  X    = rbind(X_tr,X_tst)
  betav_s = TR_RR$betav_s
  betav_o = TR_RR$betav_o
  plot(y_tst,yp_tst); abline(a=0,b=1)
  MSEP_RR =  mean((y_tst-yp_tst)**2)
  list(MSEP=MSEP_RR, betav_s = betav_s,betav_o = betav_o)
}
Pred_ols_f<-function(y_tst,X_tst,y_tr,X_tr)
{
  #OLS
  p  =   dim(X_tr)[2]
  A  = lm_f(y_tr,X_tr)
  sdv   =  apply(X_tr,2,sd2_f)
  A_inv = diag(c(1,1/sdv))
  A_inv[1,1] = 1
  A_inv[1,2:(p+1)] = -apply(X_tr,2,mean)/sdv
  betav_o = A_inv%*%A$betav_s
  
  yp_tst_ols = cbind(1,X_tst)%*%A$betav_o
  plot(y_tst,yp_tst_ols); abline(a=0,b=1)
  MSEP_ols =  mean((y_tst-yp_tst_ols)**2)
  list(MSEP=MSEP_ols,betav_s = A$betav_s)
}

#rm(list=ls(all=TRUE))
sd2_f<-function(x)
{
  (mean((x-mean(x))**2))^0.5
}
scale2_f<-function(X)
{
  scale(X,center=TRUE,scale = apply(X,2,sd2_f))
}

#RR fit for a given lambda
#internal standarized: zij = (xij-\bar{x}_j)/sqrt(ss_j)
#s2_j =  mean (xij-\bar{x}_j)^2
RR_f<-function(y,X,lamb = 0, Intercept=TRUE,Std=TRUE)
{
  p     =  dim(X)[2]; n = dim(X)[1]
  if(Std == TRUE)
  {
    sdv   =  apply(X,2,sd2_f)
    A_inv = diag(c(1,1/sdv))
    A_inv[1,1] = 1
    A_inv[1,2:(p+1)] = -apply(X,2,mean)/sdv
    X_s  =  scale2_f(X)
    Xa = X_s
    svd_X = svd(Xa); d1   = svd_X$d
    Gama1 = svd_X$v
    U = svd_X$u
    betav =  c(mean(y),Gama1%*%(((1/(d1**2+lamb))*d1)*(t(U)%*%(y-mean(y)))))
    betav_o = A_inv%*%betav# betav in orginal scale
    list(betav_s = betav,betav_o=betav_o)
  }
  
  else
  {
    Xa = scale2_f(X); p =  dim(X)[2]
    svd_X = svd(Xa)
    d   = svd_X$d
    Ind =  which(d>0)
    d1 = d[Ind]
    Gama1 = svd_X$v[,Ind]; U = svd_X$u[,Ind]
    pa    = dim(Gama1)[2]
    betav =  Gama1%*%((1/(d1**2+lamb)*d1)*(t(U)%*%(y)))#-mean(y)
    betav = betav#/c(1,sdv)# betav in orginal scale
    betav
  }
}

RR_f_V = Vectorize(RR_f,'lamb')

A_RR_f_a <-function(y_itr,X_itr,lamb)
{
  A = RR_f(y=y_itr,X=X_itr,lamb=lamb,Std=FALSE,Intercept = FALSE)
  A
}
A_RR_f_a_V<-Vectorize(A_RR_f_a,'lamb')

lamb_f<-function(Xac,K=100,li=0.001,ls=0.999)
{
  n =  dim(Xac)[1]
  R2v = seq(li,ls,length=K)
  lambv = (1-R2v)/R2v*sum(diag(Xac%*%t(Xac)))/n
  lambv = exp(seq(min(log(lambv)),max(log(lambv)),length=K))
  lambv
}

lm_f<-function(y,X)
{
  p =  dim(X)[2]
  sdv   =  apply(X,2,sd2_f)
  A_inv = diag(c(1,1/sdv))
  A_inv[1,1] = 1
  A_inv[1,2:(p+1)] = -apply(X,2,mean)/sdv
  X  =  scale2_f(X); mu =  mean(y); 
  svd_X = svd(X)
  d1   = svd_X$d; Gama1 = svd_X$v; U = svd_X$u
  betav_s =  c(mu,Gama1%*%((1/d1)*(t(U)%*%(y))))
  betav_o = A_inv%*%betav_s# betav in orginal scale
  list(betav_s = betav_s,betav_o=betav_o)
}

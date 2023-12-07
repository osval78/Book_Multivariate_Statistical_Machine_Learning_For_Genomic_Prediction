library(epiR)
source("RR.R")
Tr_RR_f<-function(y_tr,X_tr,K=5,KG=100,KR=1)
{
  n_tr =  dim(X_tr)[1]
  X_tr_s = scale2_f(X_tr)
  mu      = mean(y_tr)
  y_tr  = y_tr-mu
  #Inner CV
  lambv = lamb_f(X_tr,K=KG,li=1e-7,ls=1-1e-7)
  iPEv_mean = 0 
  for(ir in 1:KR)
  {
    iPE_mat = matrix(0,nr=length(lambv),nc=K)
    MSEP_mat=iPE_mat
    Grpv = findInterval(cut(sample(1:n_tr,n_tr),breaks=K),1:n_tr) 
    for(i in 1:K)
    {
      Pos_itr = which(Grpv!=i)
      X_itr  = X_tr_s[Pos_itr,]
      y_itr   = y_tr[Pos_itr]; 
      y_itst   = y_tr[-Pos_itr]; 
      dat_itr = data.frame(y=y_itr,X=X_itr)
      n_itr =  dim(X_itr)[1]
      betav_itr = A_RR_f_a_V(y_itr,X_itr,lambv)
      yp_mat     = 0  + (X_tr_s[-Pos_itr,])%*%betav_itr
      iPEv       = colMeans((matrix(y_itst,nr=length(y_itst),
                                    nc=length(lambv))-yp_mat)**2)
      iPE_mat[,i] =iPEv
    }
    iPEv_mean = (ir-1)/ir*iPEv_mean + rowMeans(iPE_mat)/ir
    
  }  
  
  # 
  plot(log(lambv),iPEv_mean)
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
  yp_tst  =  c(cbind(1,X_tst)%*%betav_o)
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
  betav_s =  A$betav_s
  betav_o = A_inv%*%betav_s
  yp_tst_ols = c(cbind(1,X_tst)%*%betav_o)
  MSEP_ols =  mean((y_tst-yp_tst_ols)**2)
  list(MSEP=MSEP_ols,betav_s = betav_s)
}


  

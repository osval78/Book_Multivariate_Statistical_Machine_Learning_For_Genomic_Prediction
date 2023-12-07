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

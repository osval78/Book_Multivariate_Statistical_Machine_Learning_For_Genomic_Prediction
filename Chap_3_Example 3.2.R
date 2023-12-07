rm(list=ls())
library(BMTME)
data("WheatMadaToy")
dat_F = phenoMada
dim(dat_F)
ls()
dat_F$GID = as.character(dat_F$GID) 
G =  genoMada
eig_G =  eigen(G)
G_0.5 = eig_G$vectors%*%diag(sqrt(eig_G$values))%*%t(eig_G$vectors)
X =  G_0.5
y = dat_F$PH
n =  length(y)

source('TR_RR.R')
#5FCV
set.seed(3)
K = 5
Tab=data.frame()
Grpv = findInterval(cut(sample(1:n,n),breaks=K),1:n) 
for(i in 1:K)
{
  Pos_tr = which(Grpv!=i)
  y_tr = y[Pos_tr]
  X_tr =  X[Pos_tr,]
  TR_RR = Tr_RR_f(y_tr,X_tr,K=5,KG=100,KR=1)
  lambv =  TR_RR$lambv
  #Tst
  y_tst = y[-Pos_tr]; X_tst = X[-Pos_tr,]
  #RR
  Pred_RR = Pred_RR_f(y_tst,X_tst,TR_RR)
  #OLS
  Pred_ols = Pred_ols_f(y_tst,X_tst,y_tr,X_tr)
  Tab =  rbind(Tab,data.frame(Sim=i,MSEP_RR = Pred_RR$MSEP,
                              MSEP_ols = Pred_ols$MSEP))
  cat('i = ', i,'\n')
}

colMeans(Tab)
mean(Tab$MSEP_RR<Tab$MSEP_ols)
Pos = which(Tab$MSEP_RR<Tab$MSEP_ols)
mean((Tab$MSEP_ols[Pos]-Tab$MSEP_RR[Pos])/Tab$MSEP_RR[Pos]*100)
mean((Tab$MSEP_ols[-Pos]-Tab$MSEP_RR[-Pos])/Tab$MSEP_RR[-Pos]*100)
plot(Tab$MSEP_RR,Tab$MSEP_RR_gl);abline(a=0,b=1)


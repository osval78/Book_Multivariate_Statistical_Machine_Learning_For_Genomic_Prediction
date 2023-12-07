
Data_RCBD=read.table("Example_RCBD.csv", header =T, sep = ",")
Data_RCBD

library(lme4)
Data_RCBD$Genotype=as.factor(Data_RCBD$Genotype)
Data_RCBD$Block=as.factor(Data_RCBD$Block)

Fitted=lmer(Yield~ Genotype + (1 | Block), Data_RCBD)
Fitted
####Extracting design matrix
X=Fitted@pp$X
X=X[!duplicated(X), ]
X
####Extracting the beta coeffcients of fixed effects
Beta=Fitted@beta
####EStimating the BLUEs of genotypes
BLUEs_Gen=X%*%Beta
BLUEs_Gen
#str(Fitted)

#####BLUP of genotypes
Fitted2=lmer(Yield~ (1|Genotype) + (1 | Block), Data_RCBD)
Fitted2
#####Fixed effect=Intercept######
Intercept=fixef(Fitted2)
str(Intercept)
#####Random effects of genoytpes
U_ref=c(ranef(Fitted2)$Genotype)
U_ref

########BLUP of Genotypes#####
BLUP_Gen2=Intercept+U_ref$`(Intercept)`
BLUP_Gen2
cor(c(BLUEs_Gen),BLUP_Gen2)

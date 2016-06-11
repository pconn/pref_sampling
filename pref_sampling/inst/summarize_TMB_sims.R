#Summarize TMB pref sampling sim results

load('Pref_sims_result.Rda')
Dat=Result
load('Pref_sims_result2.Rda')
Dat=rbind(Dat,Result)
load('Pref_sims_result3.Rda')
Dat=rbind(Dat,Result)
load('Pref_sims_result4.Rda')
Dat=rbind(Dat,Result)

#remove 24 results with NAs occurring in fixed effect REs
#Dat=Dat[-which(Dat[,"NA.fixed"]>0),]
#remove missinig sim records (original Results rows > # sims conducted)
Dat=Dat[-which(Dat[,1]==0),]
Dat=data.frame(Dat)
for(i in 1:4)Dat[,i]=as.factor(as.integer(Dat[,i]))

Dat$Bias=(Dat[,"N.est"]-Dat[,"N.true"])/Dat[,"N.true"]

library(doBy)
summaryBy(Bias~Gen.mod+Samp.intens+Covs+Est.mod,data=Dat,FUN='mean')
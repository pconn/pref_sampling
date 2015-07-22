#test basic spatial regression models in TMB

#  1) ICAR: Y ~ Pois( exp(X*beta + Eta) )
#  2) Thin plate spline formulation of above
#  3) Geostatistical model

#funtions and libraries for data generating model
source("./pref_sampling/R/sim_data_abundance.R")
source("./pref_sampling/R/util_funcs.R")
source("./pref_sampling/R/spat_pred.R")
library(spsurvey)
library(rgeos)
library(ggplot2)
library(RColorBrewer)

#stuff for estimation model (TMB)
library( RandomFields )
library( TMB )
library( INLA )

COMPILE=FALSE
SIM=T
prop.sampled=0.05
n.transects=90


if(COMPILE){
  Tmb.dir="./TMB_version/inst/executables/"
  TmbFile = paste0(Tmb.dir,"spatial_only.cpp")  #includes shared spatial model specification
  compile(TmbFile )
  TmbFile2 = paste0(Tmb.dir,"spatial_only_skaug.cpp")  #includes shared spatial model specification
  compile(TmbFile2 )
}

set.seed(223456)
n.sims=100 #number of simulations at each design point
S=900
Adj=rect_adj(sqrt(S),sqrt(S))
Q=-Adj
diag(Q)=apply(Adj,2,'sum')+0.01  #the +1 is just for debugging to make it proper
Q=Matrix(Q)

#initialize Results data frame
Result=matrix(0,n.sims,11)  #number of simulations * number of generating models * number of estimation models
colnames(Result)=c("Gen.mod","Samp.intens","Covs","Est.mod","Sim.no","N.est","N.true","SE","MARE","NA.fixed","N.est2")


if(SIM){   #true abundance is independent of data generating model
  for(isim in 1:n.sims){
    #simulate covariates, abundance
    Grid=sim_data_generic(S=S,n.covs=3,tau.epsilon=10000)
    Grid$cov1.quad=Grid$cov1^2
    Grid$cov2.quad=Grid$cov2^2
    Grid$cov3.quad=Grid$cov3^2
    Grid$cov12=Grid$cov1*Grid$cov2
    Grid$cov13=Grid$cov1*Grid$cov3
    Cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
    save(Grid,file=Cur.file)
  }
}

if(SIM){
  for(isim in 1:n.sims){
      cat('isim ',isim,'\n')
      Cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
      load(Cur.file)
      Effort=sim_effort(Data=Grid,S=S,type="balanced",prop.sampled=prop.sampled,n.points=n.transects)    
      Cur.file=paste("./Sim_data/Effort_bal_high",isim,".Rda",sep='')
      save(Effort,file=Cur.file)
  }
}

set.seed(423456)
for(isim in 1:n.sims){
  counter=0
  Cur.file=paste("./Sim_data/Effort_bal_high",isim,".Rda",sep='')
  load(Cur.file) #load Effort 
  Cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
  load(Cur.file)
  
  
  #reformat data to consist of R (binary sampled/not sampled) and Y (count)
  Y=Offset=rep(NA,S)
  Y[Effort$Mapping]=Effort$Counts
  Offset[Effort$Mapping]=prop.sampled
  Area.adjust=rep(1,S)   
  
  X = model.matrix(~1,data=Grid@data) 
  Data=list("n_cells"=S,"n_transects"=length(Effort$Counts),"ncol_X"=ncol(X),"Y_i"=Y, "Log_area_i"=log(Area.adjust),"Prop_sampled_i"=Offset,"X"=X,"Q"=Q)
  
  
  #### Unconditional spec
  
  TmbFile="./TMB_version/inst/executables/spatial_only"
  dyn.load(dynlib(TmbFile))
  
  Map = list()  #specify fixed parameter values
  
  Params = list("logtau_Eta"=log(1),"log_rho"=-5,"Beta"=0,"Eta"=rep(0,S))
  Random=c("Eta")  #specify which params are random effects
  #Map[["beta_pref"]] = factor( NA )
  
  obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, inner.control=list(maxit=10000) )
  obj$control <- c( obj$control, list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=1000) )
  obj$env$inner.control <- c(obj$env$inner.control, list("step.tol"=1e-8, "tol10"=1e-6, "grad.tol"=1e-8) )
  
  # First run
  Init = obj$fn( obj$par )
  Initial_gradient = obj$gr( obj$par )
  
  # Bounds
  Upper = rep(Inf, length(obj$par) )
  Lower = rep(-Inf, length(obj$par) )
  
  Upper[grep("logtau_Eta",names(obj$par))] = 8
  Lower[grep("logtau_Eta",names(obj$par))] = -8
  
  Start=obj$env$last.par.best
  if(length(Random)>0)Start=Start[-c(obj$env$random)]
  
  # Run model
  # opt = nlminb(start=Start, objective=obj$fn, gradient=obj$gr, upper=Upper, lower=Lower, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
  # opt[["final_gradient"]] = obj$gr( opt$par )
  # opt[["AIC"]] = 2*opt$objective + 2*length(opt$par)
  # #opt[["BIC"]] = 2*opt$objective + length(opt$par) * log(n.transects)
  # 
  # # Diagnostics
  # na=0
  # Report=sdreport(obj,bias.correct=TRUE)
  # Random=summary(Report,"random")
  # if(sum(is.na(Random[,2]))>0)na=1
  # 
  # Fixed=summary(Report,"fixed")
  # 
  # if(sum(is.na(Fixed[,2]))>0)na=1
  # 
  # Report2 = obj$report()
  # 
  # n.hat=Report$unbiased$value
  # n.se=Report$sd
  # n.true=sum(Grid[["N"]])
  
  
  ###########  alternative formulation using Hodges & Reich Eqn 3
  
  TmbFile="./TMB_version/inst/executables/spatial_only_skaug"
  dyn.load(dynlib(TmbFile))
  
  Map = list()  #specify fixed parameter values
  
  Params = list("logtau_Eta"=log(1),"Beta"=0,"Eta"=rep(0,S))
  Random=c("Eta")  #specify which params are random effects
  #Map[["beta_pref"]] = factor( NA )
  
  #Data2=Data[1:7]
  #Data2$W=Matrix(Adj)
  obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, inner.control=list(maxit=1000) )
  obj$control <- c( obj$control, list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100) )
  obj$env$inner.control <- c(obj$env$inner.control, list("step.tol"=1e-8, "tol10"=1e-6, "grad.tol"=1e-8) )
  
  # First run
  Init = obj$fn( obj$par )
  Initial_gradient = obj$gr( obj$par )
  
  # Bounds
  Upper = rep(Inf, length(obj$par) )
  Lower = rep(-Inf, length(obj$par) )
  
  #Upper[grep("logtau_Eta",names(obj$par))] = 8
  #Lower[grep("logtau_Eta",names(obj$par))] = -8
  
  Start=obj$env$last.par.best
  if(length(Random)>0)Start=Start[-c(obj$env$random)]
  
  # Run model
  opt = nlminb(start=Start, objective=obj$fn, gradient=obj$gr, upper=Upper, lower=Lower, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
  opt[["final_gradient"]] = obj$gr( opt$par )
  opt[["AIC"]] = 2*opt$objective + 2*length(opt$par)
  #opt[["BIC"]] = 2*opt$objective + length(opt$par) * log(n.transects)
  
  # Diagnostics
  na=0
  Report=sdreport(obj,bias.correct=TRUE)
  Random=summary(Report,"random")
  if(sum(is.na(Random[,2]))>0)na=1
  
  Fixed=summary(Report,"fixed")
  
  if(sum(is.na(Fixed[,2]))>0)na=1
  
  Report2 = obj$report()
  
  n.hat2=Report$unbiased$value
  n.se2=Report$sd
  n.true=sum(Grid[["N"]])
  
  
  
  #Result[counter,]=c(imod,isamp,icov,iest,isim,n.hat,n.true,n.se,sum(abs(Report2$Ypred_i-Grid[["N"]]))/n.true,na)
  
  #run in INLA g
  Adj.inla=Adj
  diag(Adj.inla)=1
  Dat.inla=data.frame(Y=Y,Offset=Offset,id=c(1:S))
  #Which.na=which(is.na(Y))
  #Dat.inla$Y=c(Y[-Which.na],Y[Which.na])
  #Adj2=rbind(Adj[-Which.na,],Adj[Which.na,])
  Out.inla=inla(formula=Y~1+f(id,model="besag",hyper=c(1,0.01),graph=Adj),family="poisson",data=Dat.inla,control.predictor=list(compute=T,link=1),E=Offset)
  
  n.inla.mean=sum(Out.inla$summary.fitted[,"mean"])
  n.inla.median=sum(Out.inla$summary.fitted[,"0.5quant"])
  
  Result[isim,]=c(imod,isamp,icov,iest,isim,n.hat2,n.true,n.se,sum(abs(Report2$Ypred_i-Grid[["N"]]))/n.true,na,n.inla.mean)
  
}
prop.sampled
summary((Result[,"N.est"]-Result[,"N.true"])/Result[,"N.true"])

#make some plots
# plot_N_map(N=Y,Grid,highlight=NULL,cell.width=1,leg.title="Abundance") #True abundance
# 
# plot_N_map(N=Report2$Ypred_i,Grid,highlight=NULL,cell.width=1,leg.title="Abundance") #TMB abundance
# plot_N_map(N=Grid[["N"]],Grid,highlight=NULL,cell.width=1,leg.title="Abundance") #True abundance
# plot_N_map(N=Out.inla$summary.fitted[,"mean"],Grid,highlight=NULL,cell.width=1,leg.title="Abundance") #INLA abundance
# plot_N_map(N=Report2$Ypred_i-Grid[["N"]],Grid,highlight=NULL,cell.width=1,leg.title="Abundance") #residual
# summary(Report2$Ypred_i-Grid[["N"]])
#       
    
    #save(Result,file="Pref_sims_result4.Rda")    





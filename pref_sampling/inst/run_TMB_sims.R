# run_generic_sims.R
# script to run generic spatio-temporal count data simulations
# using preferential sampling model 
# Y|R=1 ~ Pois( exp(Xy*betay + Nu + Delta + Epsilon) )
# R ~ Bernoulli( exp(Xr*betar + Gamma + Delta) )

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

COMPILE=TRUE
SIM=FALSE

if(COMPILE){
  Tmb.dir="./TMB_version/inst/executables/"
  TmbFile = paste0(Tmb.dir,"pref_samp_v3.cpp")  #includes shared spatial model specification
  compile(TmbFile )
  TmbFile = paste0(Tmb.dir,"pref_samp_linear.cpp")  #linear predictor for R_i a linear function of linear predictor for Y_i
  compile(TmbFile )
}

set.seed(123456)
n.sims=125 #number of simulations at each design point
S=900
Adj=rect_adj(sqrt(S),sqrt(S))
Q=-Adj
diag(Q)=apply(Adj,2,'sum')+0.01  #the +1 is just for debugging to make it proper
Q=Matrix(Q)

#initialize Results data frame
Result=matrix(0,n.sims*8*10,10)  #number of simulations * number of generating models * number of estimation models
colnames(Result)=c("Gen.mod","Samp.intens","Covs","Est.mod","Sim.no","N.est","N.true","SE","MARE","NA.fixed")


if(SIM){   #true abundance is independent of data generating model
  for(isim in 1:n.sims){
    #simulate covariates, abundance
    Grid=sim_data_generic(S=S,n.covs=3,tau.epsilon=10)
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
  for(isamp in 1:2){  #two different levels of sampling intensity
    if(isamp==1){
      prop.sampled=0.5
      n.transects=180
    }
    else{
      prop.sampled=0.1
      n.transects=90
    }
    
    for(isim in 1:n.sims){
      cat('isim ',isim,'\n')
      Cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
      load(Cur.file)
      Effort=sim_effort(Data=Grid,S=S,type="random",prop.sampled=prop.sampled,n.points=n.transects)
      Cur.file=paste("./Sim_data/Effort_random",isim,"samp",isamp,".Rda",sep='')
      save(Effort,file=Cur.file)
      Effort=sim_effort(Data=Grid,S=S,type="preferential",prop.sampled=prop.sampled,n.points=n.transects,pref.mult=2)    
      Cur.file=paste("./Sim_data/Effort_pref2",isim,"samp",isamp,".Rda",sep='')
      save(Effort,file=Cur.file)
      Effort=sim_effort(Data=Grid,S=S,type="preferential",prop.sampled=prop.sampled,n.points=n.transects,pref.mult=10)    
      Cur.file=paste("./Sim_data/Effort_pref10",isim,"samp",isamp,".Rda",sep='')
      save(Effort,file=Cur.file)
      Effort=sim_effort(Data=Grid,S=S,type="balanced",prop.sampled=prop.sampled,n.points=n.transects)    
      Cur.file=paste("./Sim_data/Effort_bal",isim,"samp",isamp,".Rda",sep='')
      save(Effort,file=Cur.file)
    }
  }
}

set.seed(123456)
counter=0
for(isamp in 1:2){  #sampling intensity
  if(isamp==1){
    prop.sampled=0.5
    n.transects=180
  }
  else{
    prop.sampled=0.1
    n.transects=90
  }
  for(imod in 1:4){  #data generating model
    for(isim in 1:n.sims){
      cat(paste("samp ",isamp," mod ",imod," sim ",isim,'\n'))
      #load true abundance, covariates
      Cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
      load(Cur.file) #load abundance,covariate grid
      
      #load effort data
      if(imod==1){
        Cur.file=paste("./Sim_data/Effort_random",isim,"samp",isamp,".Rda",sep='')
        load(Cur.file) #load Effort 
      }
      if(imod==2){
        Cur.file=paste("./Sim_data/Effort_pref2",isim,"samp",isamp,".Rda",sep='')
        load(Cur.file) #load Effort 
      }    
      if(imod==3){
        Cur.file=paste("./Sim_data/Effort_pref10",isim,"samp",isamp,".Rda",sep='')
        load(Cur.file) #load Effort 
      } 
      if(imod==4){
        Cur.file=paste("./Sim_data/Effort_bal",isim,"samp",isamp,".Rda",sep='')
        load(Cur.file) #load Effort 
      } 

      #reformat data to consist of R (binary sampled/not sampled) and Y (count)
      R=Y=Offset=rep(NA,S)
      R[Effort$Mapping]=1
      R[which(is.na(R))]=0
      Y[Effort$Mapping]=Effort$Counts
      Offset[Effort$Mapping]=prop.sampled
      Area.adjust=rep(1,S)   
      
      for(icov in 1:2){
        #formulate design matrices
        if(icov==1){
          Xr = model.matrix(~1,data=Grid@data) 
          Xy = model.matrix(~1,data=Grid@data) 
        }
        else{
          Xr = model.matrix(~1,data=Grid@data)
          Xy = model.matrix(~cov1+cov2+cov1.quad+cov2.quad+cov12,data=Grid@data) 
        }
        Data=list("n_cells"=S,"n_transects"=length(Effort$Counts),"ncol_Xr"=ncol(Xr), "ncol_Xy"=ncol(Xy), "R_i"=R, "Y_i"=Y, "Log_area_i"=log(Area.adjust),"Prop_sampled_i"=Offset,"Xr"=Xr, "Xy"=Xy, "Q1"=Q, "Q2"=Q,"Q3"=Q)
        
        for(iest in c(1,2,4,5)){
          counter=counter+1
          if(iest==1)TmbFile="./TMB_version/inst/executables/pref_samp_v3"
          else TmbFile="./TMB_version/inst/executables/pref_samp_linear"
          dyn.load(dynlib(TmbFile))
                   
          # Parameters
          #if(Version=="pref_samp_v3") Params = list("logtau_Eta"=log(1),"logtau_Gamma"=log(1),"logtau_Delta"=log(1),"logsigma_Epsilon"=log(1),"betar"=rep(-1,ncol(Xr)),"betay"=rep(0,ncol(Xy)),"Eta"=rep(0,S),"Gamma"=rep(0,S),"Delta"=rep(0,S),"Epsilon"=rep(0,S))
          Map = list()  #specify fixed parameter values
          if(iest==1){ #full model
            Params = list("logtau_Eta"=log(1),"logtau_Gamma"=log(1),"logtau_Delta"=log(1),"betar"=rep(-1,ncol(Xr)),"betay"=rep(0,ncol(Xy)),"Eta"=rep(0,S),"Gamma"=rep(0,S),"Delta"=rep(0,S))
            Random=c("Eta","Gamma","Delta")  #specify which params are random effects
          }
          if(iest==2){ #spatial model with linear preferential sampling
            Params = list("logtau_Eta"=log(1),"logtau_Gamma"=log(1),"beta_pref"=0,"betar"=rep(-1,ncol(Xr)),"betay"=rep(0,ncol(Xy)),"Eta"=rep(0,S),"Gamma"=rep(0,S))
            Random=c("Eta","Gamma")  #specify which params are random effects
          }
          #if(iest==3){ #nonspatial model with linear preferential sampling  ###this one is overparamterized I think
          #  Params = list("logtau_Eta"=log(1),"logtau_Gamma"=log(1),"beta_pref"=0,"betar"=rep(-1,ncol(Xr)),"betay"=rep(0,ncol(Xy)),"Eta"=rep(0,S),"Gamma"=rep(0,S))
          #  Random=NULL  #specify which params are random effects
          #  Map[["Gamma"]] = factor( rep(NA, length(Params[["Gamma"]])) )  #fix parameter values
          #  Map[["logtau_Gamma"]] = factor( NA )
          #  Map[["Eta"]] = factor( rep(NA, length(Params[["Eta"]])) )
          #  Map[["logtau_Eta"]] = factor( NA ) 
          #}
          if(iest==4){ #spatial model without preferential sampling
            Params = list("logtau_Eta"=log(1),"logtau_Gamma"=log(1),"beta_pref"=0,"betar"=rep(-1,ncol(Xr)),"betay"=rep(0,ncol(Xy)),"Eta"=rep(0,S),"Gamma"=rep(0,S))
            Random=c("Eta","Gamma")  #specify which params are random effects
            Map[["beta_pref"]] = factor( NA )
          }          
          if(iest==5){ #nonspatial model without preferential sampling
            Params = list("logtau_Eta"=log(1),"logtau_Gamma"=log(1),"beta_pref"=0,"betar"=rep(-1,ncol(Xr)),"betay"=rep(0,ncol(Xy)),"Eta"=rep(0,S),"Gamma"=rep(0,S))
            Random=NULL  #specify which params are random effects
            Map[["beta_pref"]] = factor( NA )
            Map[["Gamma"]] = factor( rep(NA, length(Params[["Gamma"]])) )  #fix parameter values
            Map[["logtau_Gamma"]] = factor( NA )
            Map[["Eta"]] = factor( rep(NA, length(Params[["Gamma"]])) )
            Map[["logtau_Eta"]] = factor( NA )            
          }           
          
          if(iest<6){   #(iest != 3){  #this one is analytically not identifiable when no covs; probably not even when there are covariates... (Hessian for Bernoulli fixed effects not pos def)
            obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, inner.control=list(maxit=1000) )
            obj$control <- c( obj$control, list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100) )
            obj$env$inner.control <- c(obj$env$inner.control, list("step.tol"=1e-8, "tol10"=1e-6, "grad.tol"=1e-8) )
            
            # First run
            Init = obj$fn( obj$par )
            Initial_gradient = obj$gr( obj$par )
            
            # Bounds
            Upper = rep(Inf, length(obj$par) )
            Lower = rep(-Inf, length(obj$par) )
            
            Upper[grep("logtau_Gamma",names(obj$par))] = 8
            Lower[grep("logtau_Gamma",names(obj$par))] = -8
            Upper[grep("logtau_Eta",names(obj$par))] = 8
            Lower[grep("logtau_Eta",names(obj$par))] = -8
            Upper[grep("logtau_Delta",names(obj$par))] = 8
            Lower[grep("logtau_Delta",names(obj$par))] = -8
                        
            Start=obj$env$last.par.best
            if(length(Random)>0)Start=Start[-c(obj$env$random)]
            
            # Run model
            opt = nlminb(start=Start, objective=obj$fn, gradient=obj$gr, upper=Upper, lower=Lower, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
            opt[["final_gradient"]] = obj$gr( opt$par )
            opt[["AIC"]] = 2*opt$objective + 2*length(opt$par)
            #opt[["BIC"]] = 2*opt$objective + length(opt$par) * log(n.transects)
            
            # Diagnostics
            na=0
            if(length(Random)>0){
              Report=sdreport(obj,bias.correct=TRUE)
              Random=summary(Report,"random")
              if(sum(is.na(Random[,2]))>0)na=1
            }
            else{
              Report=sdreport(obj)
            }
            Fixed=summary(Report,"fixed")
            
             if(sum(is.na(Fixed[,2]))>0)na=1
            
            Report2 = obj$report()
            
            n.hat=Report$value
            n.se=Report$sd
            n.true=sum(Grid[["N"]])
            
            Result[counter,]=c(imod,isamp,icov,iest,isim,n.hat,n.true,n.se,sum(abs(Report2$Ypred_i-Grid[["N"]]))/n.true,na)
            
            
          }
        } 
      }
    }
    save(Result,file="Pref_sims_result.Rda")    
  }
}

    
    # Data
  
    #if(Version=="pref_samp_v2") Params = list("log_kappa"=rep(log(1),3), "logtau_Nu"=log(1), "logtau_Gamma"=log(1), "logtau_Delta"=log(1), "betar"=rep(0,ncol(Xr)), "betay"=rep(0,ncol(Xy)), "Nu_input"=rep(0,Data$n_knots), "Gamma_input"=rep(0,Data$n_knots), "Delta_input"=rep(0,Data$n_knots) )
    
    
      
#       # Initialization
#       obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, inner.control=list(maxit=1000) )
#       obj$control <- c( obj$control, list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100) )
#       obj$env$inner.control <- c(obj$env$inner.control, list("step.tol"=1e-8, "tol10"=1e-6, "grad.tol"=1e-8) )
#       
#       # First run
#       Init = obj$fn( obj$par )
#       Initial_gradient = obj$gr( obj$par )
#       
#       # Bounds
#       Upper = rep(Inf, length(obj$par) )
#       Lower = rep(-Inf, length(obj$par) )
#       Start=obj$env$last.par.best[-c(obj$env$random)]
#         
#       # Run model
#       opt = nlminb(start=Start, objective=obj$fn, gradient=obj$gr, upper=Upper, lower=Lower, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
#       opt[["final_gradient"]] = obj$gr( opt$par )
#       opt[["AIC"]] = 2*opt$objective + 2*length(opt$par)
#       #opt[["BIC"]] = 2*opt$objective + length(opt$par) * log(n.transects)
#       
#       # Diagnostics
#       Report=sdreport(obj,bias.correct=TRUE)
#       Fixed=summary(Report,"fixed")
#       Random=summary(Report,"random")
# 
#       Report2 = obj$report()
# 
#       n.hat=Report$value
#       n.se=Report$sd
#       #plot maps
#       #plot_N_map(N=Grid[["N"]],Grid,highlight=NULL,cell.width=1,leg.title="Abundance") #True abundance
#       #plot_N_map(N=Report2$Ypred_i,Grid,highlight=NULL,cell.width=1,leg.title="Abundance") #True abundance
# 
#       
#       
#       #Results[modelI,i,c("Sigma_Nu","Sigma_Gamma","Sigma_Delta",paste0("Range_",c("Nu","Gamma","Delta")))] = unlist( Report[c("Sigma_Nu","Sigma_Gamma","Sigma_Delta","Range")] )
#       #if(modelI==2) Results[modelI,i,c("Sigma_Delta","Range_Delta")] = NA
#       
#       # Calculate error in prediction of Ypred
#       #Ypred_error = sum(Report[["Ypred_i"]] - Ypred) / sum(Ypred)
#       #Results[modelI,i,"Ypred_error"] = Ypred_error
#     #}






# Z -- true abundance
# Y -- sampled abundance
# R -- predicted probability of sampling
# delta -- affects density
# eta -- affects sampling probability

library( RandomFields )
library( TMB )
library( INLA )
library(PrefSampling)

#library( ThorsonUtilities )

Version = "simple_v6"
# v2 -- added option for changing distribution for random effects, and implemented ICAR distribution
# v3 -- added Bernoulli model for sampling locations, and removed the ICAR functionality (by commenting it out, to save speed during building the object)
# v4 -- added simpler less-smooth SPDE option (may be faster!), and option to turn off delta or eta
# v5 -- Implements Devin's finite population correction approach to posterior prediction; proportion area sampled additional input vector  
# v6 -- Fixes link function for R_i; now logit() instead of multinomial-logit()

# Compile
#TmbFile = "C:/Users/paul.conn/git/pref_sampling/pref_sampling/src/PrefSampling"
#compile(paste0(TmbFile,".cpp")) 


# Settings
grid_dim = c("x"=25, "y"=25)
n_samp = 50
n_cells = grid_dim[1]*grid_dim[2]
prop_sampled=0.5
Prop_sampled=rep(prop_sampled,n_samp)
SpatialScale = sqrt(prod(grid_dim))/5  # Range ~ 2*Scale
SD_eta = SD_x = SD_delta = 1
beta0 = 2
Use_REML = FALSE   # 
Spatial_sim_model = c("GP_gaussian", "ICAR")[1]
Spatial_model = c("SPDE_GMRF", "ICAR")[1]
Alpha = 1  # Smoothness for GMRF, 1 or 2 (1 is faster)
RandomSeed = 123456
n_sim = 500

# Configurations
b_set = c(0,1,3)            # Impact of delta on sampling intensity
EM_set = c("fix_b","est_b")

# Configurations

# Source ICAR precision matrix function
#Q = rect_adj( x=grid_dim['x'], y=grid_dim['y'] )
#Q = -1*Q
#diag(Q) = -1 * colSums(Q)

# Results
Results = data.frame(matrix(NA,length(b_set)*length(EM_set)*n_sim,12))
colnames(Results)=c("sim","b","EM","b_est","SEb","N_true","N_est","N_SE","Convergence","LogLik","BiasBeta11","BiasBeta12")

# Loop through replicates
counter=1
for(i in 1:n_sim){
  for(SimI in 1:length(b_set)){
    for(EstI in 1:length(EM_set)){
      cat(paste("Simulation ",i," SimI ",SimI," EstI ",EstI,"\n"))
      # Settings
      b = b_set[SimI]
      EM = EM_set[EstI]
      set.seed( RandomSeed + i )
      
      #note: betax set to zero in last commit by Jim (v4)
      betax = betax_prob = runif(1,-.5,.5)       # Impact of x on density
      # did this so places that have covariate values with higher abundance are more likely to be sampled
      
      # Spatial model
      loc_s = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
      model_delta <- RMgauss(var=SD_delta^2, scale=SpatialScale)
      model_x <- RMgauss(var=SD_x^2, scale=SpatialScale)
      model_eta <- RMgauss(var=SD_eta^2, scale=SpatialScale)
      
      # Realization from GRF
      # delta_s -- spatial variation that affects both density and sampling intensity
      # x_s -- spatial variation in 
      if( Spatial_sim_model=="GP_gaussian"){
        delta_s = RFsimulate(model=model_delta, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
        x_s = RFsimulate(model=model_x, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
        eta_s = RFsimulate(model=model_eta, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
      }
      if( Spatial_sim_model=="ICAR"){
        delta_s = rrw( Q )[,1]
        x_s = rrw( Q )[,1]
        eta_s = rrw( Q )[,1]
      }
      
      # Total abundance
      Ztrue_s = exp( beta0 + betax*x_s + delta_s )
      Ztrue_s = rpois(n_cells,Ztrue_s)     #added to v5, 11/9/2015
      
      # Samping intensity
      R_s = exp( betax_prob*x_s + eta_s + b*delta_s )
      R_s = R_s / sum(R_s)
      
      # Process for locating samples
      s_i = sample(1:prod(grid_dim), size=n_samp, replace=FALSE, prob=R_s)
      y_s = ifelse(1:prod(grid_dim) %in% s_i, 1, 0)
      
      # Counting process
      #c_i = rpois( n=n_samp, lambda=Ztrue_s[s_i])
      c_i = rbinom(n=n_samp,Ztrue_s[s_i],prop_sampled)   #changed to binom for v5, 11/9
      
      # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
      mesh = inla.mesh.create( loc_s )
      
      # Plot stuff
      #par( mfrow=c(1,2) )
      #plot( x=Ztrue_s[s_i], y=c_i)
      #plot( x=Ztrue_s, y=R_s)
      
      # Loop 
      # Options
      Options_vec = c( 'Prior'=switch(Spatial_model,"ICAR"=1,"SPDE_GMRF"=0), 'Alpha'=Alpha, 'IncludeDelta'=1, 'IncludeEta'=1,"SE"=1)
      
      # Data
      spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
      #Data = list( "Options_vec"=Options_vec, "c_i"=c_i, "P_i"=Prop_sampled,"A_s"=rep(1,n_cells),"s_i"=s_i-1, "X_sj"=cbind(1,x_s), "y_s"=y_s, "X_sk"=cbind(1,x_s),"X_sb"=matrix(1,n_cells,1), "spde"=spde)
      Data = list( "Options_vec"=Options_vec, "c_i"=c_i, "P_i"=Prop_sampled,"A_s"=rep(1,n_cells),"s_i"=s_i-1, "X_sj"=cbind(1,x_s), "y_s"=y_s, "X_sk"=cbind(1,x_s),"X_sb"=matrix(1,n_cells,1), "spde"=spde)
      
      
      # Parameters
      # Intercept is needed for beta_j (delta -- abundance) but not beta_k (eta -- sampling intensity)
      if( Options_vec['Prior']==0 ) etainput_s = deltainput_s = rep(0,mesh$n)
      if( Options_vec['Prior']==1 ) etainput_s = deltainput_s = rep(0,prod(grid_dim))
      Params = list("beta_j"=rep(0,ncol(Data$X_sj)), "beta_k"=rep(0,ncol(Data$X_sk)), "beta_b"=0, "logtau_z"=rep(0,2), "logkappa_z"=rep(0,2), "deltainput_s"=deltainput_s, "etainput_s"=etainput_s)
      Params$beta_j[1]=beta0 + runif(1,-0.1,0.1) #set mean expected abundance close to truth for faster optimization
      if(ncol(Data$X_sj)>1)Params$beta_j[2]=betax + runif(1,-0.1,0.1)
      Params$beta_k[1]=log(n_samp/(n_cells-n_samp)) #just proportion sampled
      if(ncol(Data$X_sk)>1)Params$beta_k[2]=betax_prob + runif(1,-0.1,0.1)
  
      # Random
      if(Version %in% c("simple_v2","simple_v1")) Random = c( "nuinput_s" )
      if(Version %in% c("simple_v6","simple_v5","simple_v4","simple_v3")) Random = c( "deltainput_s", "etainput_s" )
      if(Use_REML==TRUE) Random = c(Random,"beta_j","beta_k")
      if(Use_REML==TRUE & Version%in%c("simple_v6","simple_v5","simple_v4","simple_v3")) Random = c( Random, "beta_b" )
      
      # Fix parameters
      Map = list()
      # Fix common spatial-scale
      Map[["logkappa_z"]] = factor( rep(1,length(Params[["logkappa_z"]])) )
      # Elminate intercept and marginal-SD if using ICAR model
      if( Options_vec['Prior']==1 ){
        Params[["logtau"]] = log(1)
        Map[["logtau"]] = factor(NA)
        Params[["beta_j"]] = log(1)
        Map[["beta_j"]] = factor(NA)
      }
      # Eliminate linkeage of density and sampling intensity
      if( EM=="fix_b" ){
        Map[["beta_b"]] = factor(NA)
      }
      
      # Make object
      #compile( paste0(Version,".cpp") )
      #dyn.load( dynlib(TmbFile) )
      Start_time = Sys.time()
      Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE, DLL="PrefSampling")
      Obj$fn( Obj$par )
      
      # Run
      #Lower = -Inf
      #Upper = Inf
      Lower = -50  #getting a few cases where -Inf,Inf bounds result in nlminb failure (NaN gradient)
      Upper = 50
      Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
      Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
      Report = Obj$report()
      
      # Potentially fix random fields with zero sample or population variance
      if( any(Report$MargSD_z<0.001) ){
        Which = which(Report$MargSD_z<0.001)
        Map[["logtau_z"]] = factor( ifelse(1:2==Which,NA,1:2) )
        if(length(Which)==2){
          Map[["logkappa_z"]] = factor( c(NA,NA) )
        }
        if( any(Which==1) ){
          Map[["deltainput_s"]] = factor( rep(NA,length(Params[["deltainput_s"]])) )
          Params[["deltainput_s"]][] = 0
          Map[["beta_b"]] = factor(NA)
          Params[["beta_b"]] = 0
        }
        if( any(Which==2) ){
          Map[["etainput_s"]] = factor( rep(NA,length(Params[["etainput_s"]])) )
          Params[["etainput_s"]][] = 0
        }
        Data$Options_vec[Which+2] = 0
        # Re-run
        if( length(Which)!=2 ) Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE,DLL="PrefSampling")
        if( length(Which)==2 ) Obj = MakeADFun( data=Data, parameters=Params, random=NULL, map=Map, silent=TRUE,DLL="PrefSampling")
        Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
        Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
      }
      Converge=Opt$convergence
      
      # SD
      if(Converge==0){
        Report = Obj$report()
        if( all(c("etainput_s","deltainput_s")%in%names(Map)) ){
          SD = sdreport( Obj, bias.correct=FALSE )
          SD$unbiased$value = c("total_abundance"=Report$total_abundance)
        }else{
          SD = sdreport( Obj, bias.correct=TRUE )  
        }
        Opt[["run_time"]] = Sys.time()-Start_time
      }
      
      # Record results
      if(EstI==2){
        which_beta_b = which(colnames(SD$value)=="beta_b")
        Results[counter,1:9]=c(i,SimI,EstI,Report$beta_b,SD$sd[which(names(SD$value)=="beta_b")],sum(Ztrue_s), SD$unbiased$value['total_abundance'],SD$sd[which(names(SD$value)=="total_abundance")],Converge)
      }
      if(EstI==1){
        Results[counter,1:9]=c(i,SimI,EstI,0,NA,sum(Ztrue_s), SD$unbiased$value['total_abundance'],SD$sd[which(names(SD$value)=="total_abundance")],Converge)
      }
      Results[counter,"LogLik"]=-Opt$objective
      Results[counter,"BiasBeta11"]=Opt$par[1]-beta0
      Results[counter,"BiasBeta12"]=Opt$par[2]-betax
      
      #dyn.unload( dynlib(TmbFile) )
      
      # output
      save( Results, file="Results.RData")
      
      counter=counter+1
      # Plot stuff
      #par( mfrow=c(1,3) )
      #plot( x=Ztrue_s[s_i], y=c_i)
      #plot( x=Ztrue_s, y=R_s)
      #plot( x=Report$Z_s, y=Report$R_s)
    }
  }     
}

####################
# Read results
####################
TmbFile = "C:/Users/paul.conn/git/pref_sampling/" #TMB_version/inst/executables"
setwd(TmbFile)
load( "Results.RData")
WhichDone = which(Results[,"Convergence"]==0)
Results=Results[WhichDone,]
Results[,"RelBias_N"]=(Results[,"N_est"]-Results[,"N_true"])/Results[,"N_true"]
True.b = Results[,"b"]
True.b[which(True.b==1)]=0
True.b[which(True.b==2)]=1
True.b[which(True.b==3)]=3
Results[,"b"]=True.b
Which.est = which(Results[,"EM"]>1)
Results[,"Bias_b"]=NA
Results[Which.est,"Bias_b"]=(Results[Which.est,"b_est"]-Results[Which.est,"b"])

#
pdf("B_est.pdf")
  Which_plot = which(Results[,"EM"]>1)
  plot(Results[Which_plot,"b_est"],xlab="Simulation",ylab=expression(hat(b)))
dev.off()


#For figures, limit to results for which |b estimate| < 10 
Results[,"OnBoundary"]=rep(0,nrow(Results))
Results[which(abs(Results[,"b_est"])>10),"OnBoundary"]=1

Results_plot=Results[Results[,"OnBoundary"]==0,]
Results_plot[,"Est.model"]=rep("Independent",nrow(Results_plot))
Results_plot[Results_plot[,"EM"]==2,"Est.model"]="Joint"

library(doBy)
sd <- function(x)sqrt(var(x))
Bias_table = summaryBy(Bias_b+RelBias_N~b+EM,data=Results_plot,FUN=c(median,mean,sd))

#plot bias_N by b_est for EstMod=2
Cur_data=Results_plot[Results_plot[,"EM"]==2,]
pdf("Bias_bN.pdf")
plot(Cur_data[,"Bias_b"],Cur_data[,"RelBias_N"],xlab="Bias (b)",ylab="Relative bias (N)")
dev.off()



#plot bias as function of b, estimation method
#first, rearrange relative bias in form for ggplot
Bias.df=Results_plot
Bias.df[,"B"]=Bias.df[,"b"]
Bias.df[,"Bias"]=Bias.df[,"RelBias_N"]


library(ggplot2)
#plot proportion relative bias
bias.plot = ggplot(Bias.df,aes(factor(Est.model),Bias))+geom_boxplot()+facet_grid(~b) #,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
#bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
#bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + labs(x = "Estimation model", y="Proportion relative bias")
bias.plot
pdf("bias.pdf")
bias.plot
dev.off()

#bias of b parameter
BiasB.df = Bias.df[which(!is.na(Bias.df[,"Bias_b"])),]
bias.plot = ggplot(BiasB.df,aes(factor(b),Bias_b))+geom_boxplot() #+facet_grid(~b) #,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
#bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
#bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + labs(x = "True b parameter", y="Absolute bias")
bias.plot
pdf("biasB.pdf")
bias.plot
dev.off()

#produce plot of bias of species-habitat relationship parameters
bias.plot = ggplot(Bias.df,aes(factor(Est.model),BiasBeta11))+geom_boxplot()+facet_grid(~b) #,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
#bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
#bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
bias.plot
pdf("Bias_beta0.pdf")
bias.plot
dev.off()

#produce plot of bias of species-habitat relationship parameters
bias.plot = ggplot(Bias.df,aes(factor(Est.model),BiasBeta12))+geom_boxplot()+facet_grid(~b) #,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20)) #+ coord_cartesian(ylim = c(-1., 6.0))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
#bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
#bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
bias.plot
pdf("Bias_beta1.pdf")
bias.plot
dev.off()


Results[,"Converge"]=Results[,"OnBoundary"]+Results[,"Convergence"]
Results[which(Results[,"Converge"]==2),"Converge"]=1
Converge_table = summaryBy(Converge~b+EM,data=Results,FUN=sum)



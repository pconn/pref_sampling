
# Z -- true abundance
# Y -- sampled abundance
# R -- predicted probability of sampling
# delta -- affects density
# eta -- affects sampling probability

library( RandomFields )
library( TMB )
library( INLA )
library( ThorsonUtilities )

Version = "simple_v4"
# v2 -- added option for changing distribution for random effects, and implemented ICAR distribution
# v3 -- added Bernoulli model for sampling locations, and removed the ICAR functionality (by commenting it out, to save speed during building the object)
# v4 -- added simpler less-smooth SPDE option (may be faster!), and option to turn off delta or eta

# Compile
TmbFile = "C:/Users/James.Thorson/Desktop/Project_git/pref_sampling/TMB_version/inst/executables"
setwd( TmbFile )
compile( paste0(Version,".cpp") )

# SOurce
source( paste0(TmbFile,"/../../examples/rect_adj.R") )
source( paste0(TmbFile,"/../../examples/rrw.R") )

# Settings
grid_dim = c("x"=30, "y"=30)
n_samp = 50
SpatialScale = sqrt(prod(grid_dim))/5  # Range ~ 2*Scale
SD_eta = SD_x = SD_delta = 1
beta0 = 2
betax = 0        # Impact of x on density
betax_prob = 0   # Impact of x on sampling intensity
Use_REML = FALSE   # 
Spatial_sim_model = c("GP_gaussian", "ICAR")[1]
Spatial_model = c("SPDE_GMRF", "ICAR")[1]
Alpha = 1  # Smothness for GMRF, 1 or 2 (1 is faster)
RandomSeed = ceiling(runif(1,min=1,max=1e6))
n_sim = 500

# Configurations
b_set = c(0,5)            # Impact of delta on sampling intensity
EM_set = c("fix_b","est_b")

# Configurations

# Source ICAR precision matrix function
Q = rect_adj( x=grid_dim['x'], y=grid_dim['y'] )
Q = -1*Q
diag(Q) = -1 * colSums(Q)

# Results
Results = array(NA, dim=c(length(b_set),length(EM_set),n_sim,3), dimnames=list(paste0("b=",b_set),EM_set,NULL,c("b","sum_Ztrue","sum_Zpred")) )

# Loop through replicates
for(i in 1:n_sim){
for(SimI in 1:length(b_set)){
for(EstI in 1:length(EM_set)){
 
  # Settings
  b = b_set[SimI]
  EM = EM_set[EstI]
  set.seed( RandomSeed + i )
  
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

  # Samping intensity
  R_s = exp( betax_prob*x_s + eta_s + b*delta_s )
  R_s = R_s / sum(R_s)
  
  # Process for locating samples
  s_i = sample(1:prod(grid_dim), size=n_samp, replace=FALSE, prob=R_s)
  y_s = ifelse(1:prod(grid_dim) %in% s_i, 1, 0)
  
  # Counting process
  c_i = rpois( n=n_samp, lambda=Ztrue_s[s_i])

  # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
  mesh = inla.mesh.create( loc_s )
  
  # Plot stuff
  par( mfrow=c(1,2) )
  plot( x=Ztrue_s[s_i], y=c_i)
  plot( x=Ztrue_s, y=R_s)

  # Loop 
  # Options
  Options_vec = c( 'Prior'=switch(Spatial_model,"ICAR"=1,"SPDE_GMRF"=0), 'Alpha'=Alpha, 'IncludeDelta'=1, 'IncludeEta'=1 )

  # Data
  spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
  if(Version%in%c("simple_v1")) Data = list( "c_i"=c_i, "s_i"=s_i-1, "X_sj"=Zy, "spde"=spde)
  if(Version%in%c("simple_v2")) Data = list( "Options_vec"=Options_vec, "c_i"=c_i, "s_i"=s_i-1, "X_sj"=cbind(1,x_s), "spde"=spde, "Q_ICAR"=Matrix(Q))
  if(Version%in%c("simple_v4","simple_v3")) Data = list( "Options_vec"=Options_vec, "c_i"=c_i, "s_i"=s_i-1, "X_sj"=cbind(1,x_s), "y_s"=y_s, "X_sk"=cbind(x_s), "spde"=spde)

  # Parameters
  # Intercept is needed for beta_j (delta -- abundance) but not beta_k (eta -- sampling intensity)
  if( Options_vec['Prior']==0 ) etainput_s = deltainput_s = rep(0,mesh$n)
  if( Options_vec['Prior']==1 ) etainput_s = deltainput_s = rep(0,prod(grid_dim))
  if(Version %in% c("simple_v2","simple_v1")) Params = list("beta_j"=rep(0,ncol(Data$X_sj)), "logtau"=log(1), "logkappa"=log(1), "nuinput_s"=nuinput_s )
  if(Version %in% c("simple_v4","simple_v3")) Params = list("beta_j"=rep(0,ncol(Data$X_sj)), "beta_k"=rep(0,ncol(Data$X_sk)), "b"=0, "logtau_z"=rep(0,2), "logkappa_z"=rep(0,2), "deltainput_s"=deltainput_s, "etainput_s"=etainput_s )

  # Random
  if(Version %in% c("simple_v2","simple_v1")) Random = c( "nuinput_s" )
  if(Version %in% c("simple_v4","simple_v3")) Random = c( "deltainput_s", "etainput_s" )
  if(Use_REML==TRUE) Random = c(Random,"beta_j","beta_k")
  if(Use_REML==TRUE & Version%in%c("simple_v4","simple_v3")) Random = c( Random, "b" )

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
    Map[["b"]] = factor(NA)
  }
  
  # Make object
  #compile( paste0(Version,".cpp") )
  dyn.load( dynlib(Version) )
  Start_time = Sys.time()
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE)
  Obj$fn( Obj$par )

  # Run
  Lower = -Inf
  Upper = Inf
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
      Map[["b"]] = factor(NA)
      Params[["b"]] = 0
    }
    if( any(Which==2) ){
      Map[["etainput_s"]] = factor( rep(NA,length(Params[["etainput_s"]])) )
      Params[["etainput_s"]][] = 0
    }
    Data$Options_vec[Which+2] = 0
    # Re-run
    if( length(Which)!=2 ) Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE)
    if( length(Which)==2 ) Obj = MakeADFun( data=Data, parameters=Params, random=NULL, map=Map, silent=TRUE)
    Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, maxit=1000))         #
    Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  }

  # SD
  Report = Obj$report()
  if( all(c("etainput_s","deltainput_s")%in%names(Map)) ){
    SD = sdreport( Obj, bias.correct=FALSE )
    SD$unbiased$value = c("total_abundance"=Report$total_abundance)
  }else{
    SD = sdreport( Obj, bias.correct=TRUE )  
  }
  Opt[["run_time"]] = Sys.time()-Start_time
  
  # Record                  # "betay","Sigma_Nu","Range_Nu",
  Results[SimI,EstI,i,c("b","sum_Ztrue","sum_Zpred")] = c( Report$b, sum(Ztrue_s), SD$unbiased$value['total_abundance'] )
  dyn.unload( dynlib(Version) )
  
  # output
  write.table( aperm(Results,c(3,4,1,2)), file="Results.txt")            # aperm(Results,c(2,3,1))
  save( Results, file="Results.RData")

  # Plot stuff
  par( mfrow=c(1,3) )
  plot( x=Ztrue_s[s_i], y=c_i)
  plot( x=Ztrue_s, y=R_s)
  plot( x=Report$Z_s, y=Report$R_s)

}}}

####################
# Read results
####################
TmbFile = "C:/Users/James.Thorson/Desktop/Project_git/pref_sampling/TMB_version/inst/executables"
Stats = function(vec){ c("mean"=mean(vec,na.rm=TRUE),"median"=median(vec,na.rm=TRUE),"SE"=sqrt(var(vec,na.rm=TRUE)/sum(!is.na(vec))))}
setwd(TmbFile)
load( "Results.RData")
WhichDone = which( !is.na(Results[1,1,,'sum_Ztrue']))
length(WhichDone)

# Print to screen
Results[,,WhichDone,]
apply(Results[,,WhichDone,], MARGIN=c(1:2,4), FUN=Stats)
Stats( Results[,'sum_Zpred',drop=FALSE]/Results[,'sum_Ztrue',drop=FALSE] )
apply(Results[,,,'sum_Zpred',drop=FALSE]/Results[,,,'sum_Ztrue',drop=FALSE], MARGIN=c(1:2), FUN=Stats)

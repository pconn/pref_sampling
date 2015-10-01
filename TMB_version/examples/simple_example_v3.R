
# Z -- true abundance
# Y -- sampled abundance
# R -- predicted probability of sampling
# delta -- affects density
# eta -- affects sampling probability

library( RandomFields )
library( TMB )
library( INLA )

Version = "simple_v3"
# v2 -- added option for changing distribution for random effects, and implemented ICAR distribution
# v3 -- added Bernoulli model for sampling locations, and removed the ICAR functionality (by commenting it out, to save speed during building the object)

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
betax = 1        # Impact of x on density
betax_prob = 1   # Impact of x on sampling intensity
b = 0            # Impact of delta on sampling intensity
Use_REML = TRUE
Spatial_sim_model = c("GP_gaussian", "ICAR")[1]
Spatial_model = c("SPDE_GMRF", "ICAR")[1]

# Source ICAR precision matrix function
Q = rect_adj( x=grid_dim['x'], y=grid_dim['y'] )
Q = -1*Q
diag(Q) = -1 * colSums(Q)
# Sigma = solve(Q)   # Appears to be singular

# Calculate C (pg. 185 Cressie and Wikle)
#C = rect_adj( x=grid_dim['x'], y=grid_dim['y'] )

# Results
Results = array(NA, dim=c(500,3), dimnames=list(NULL,c("b","sum_Ztrue","sum_Zpred")) )

# Loop through replicates
for(i in 1:dim(Results)[1]){
 
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
  Ztrue_s = exp( beta0 + betax*x_s )

  # Samping intensity
  R_s = exp( betax_prob + betax_prob*x_s + eta_s + b*delta_s )
  R_s = R_s / sum(R_s)
  
  # Process for locating samples
  s_i = sample(1:prod(grid_dim), size=n_samp, replace=FALSE, prob=R_s)
  y_s = ifelse(1:prod(grid_dim) %in% s_i, 1, 0)
  
  # Counting process
  c_i = rpois( n=n_samp, lambda=Ztrue_s[s_i])

  # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
  mesh = inla.mesh.create( loc_s )
  plot( x=Ztrue_s[s_i], y=c_i)

  # Loop 
  # Options
  Options_vec = c( 'Prior'=switch(Spatial_model,"ICAR"=1,"SPDE_GMRF"=0) )

  # Data
  spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
  if(Version%in%c("simple_v1")) Data = list( "c_i"=c_i, "s_i"=s_i-1, "X_sj"=Zy, "spde"=spde)
  if(Version%in%c("simple_v2")) Data = list( "Options_vec"=Options_vec, "c_i"=c_i, "s_i"=s_i-1, "X_sj"=cbind(1,x_s), "spde"=spde, "Q_ICAR"=Matrix(Q))
  if(Version%in%c("simple_v3")) Data = list( "Options_vec"=Options_vec, "c_i"=c_i, "s_i"=s_i-1, "X_sj"=cbind(1,x_s), "y_s"=y_s, "X_sk"=cbind(x_s), "spde"=spde)

  # Parameters
  # Intercept is needed for beta_j (delta -- abundance) but not beta_k (eta -- sampling intensity)
  if( Options_vec['Prior']==0 ) etainput_s = deltainput_s = rep(0,mesh$n)
  if( Options_vec['Prior']==1 ) etainput_s = deltainput_s = rep(0,prod(grid_dim))
  if(Version %in% c("simple_v2","simple_v1")) Params = list("beta_j"=rep(0,ncol(Data$X_sj)), "logtau"=log(1), "logkappa"=log(1), "nuinput_s"=nuinput_s )
  if(Version %in% c("simple_v3")) Params = list("beta_j"=rep(0,ncol(Data$X_sj)), "beta_k"=rep(0,ncol(Data$X_sk)), "b"=0, "logtau_z"=rep(0,2), "logkappa_z"=rep(0,2), "deltainput_s"=deltainput_s, "etainput_s"=etainput_s )

  # Random
  if(Version %in% c("simple_v2","simple_v1")) Random = c( "nuinput_s" )
  if(Version %in% c("simple_v3")) Random = c( "deltainput_s", "etainput_s" )
  if(Use_REML==TRUE) Random = c(Random,"beta_j","beta_k")

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
  
  # Make object
  #compile( paste0(Version,".cpp") )
  dyn.load( dynlib(Version) )
  Start_time = Sys.time()
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE)
  Obj$fn( Obj$par )

  # Run
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list(trace=1, maxit=1000))         #
  Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  Opt[["run_time"]] = Sys.time()-Start_time

  # SD
  SD = sdreport( Obj, bias.correct=TRUE)
  Report = Obj$report()
  
  # Record                  # "betay","Sigma_Nu","Range_Nu",
  Results[i,c("b","sum_Ztrue","sum_Zpred")] = c( Report$b, sum(Ztrue_s), SD$unbiased$value['total_abundance'] )
  dyn.unload( dynlib(Version) )
  
  # output
  write.table( Results, file="Results.txt")            # aperm(Results,c(2,3,1))
  save( Results, file="Results.RData")
}

####################
# Read results
####################
TmbFile = "C:/Users/James.Thorson/Desktop/Project_git/pref_sampling/TMB_version/inst/executables"
setwd(TmbFile)
load( "Results.RData")
WhichDone = which( !is.na(Results[,'sum_Ztrue']))

# Print to screen
Results[WhichDone,]
Stats = function(vec){ c("mean"=mean(vec,na.rm=TRUE),"median"=median(vec,na.rm=TRUE),"SE"=sqrt(var(vec,na.rm=TRUE)/sum(!is.na(vec))))}
apply(Results, MARGIN=c(2), FUN=Stats)
#apply(, MARGIN=1, FUN=Stats)
Stats( Results[,'sum_Zpred',drop=FALSE]/Results[,'sum_Ztrue',drop=FALSE] )
length(WhichDone)

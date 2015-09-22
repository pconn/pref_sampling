
# Z -- true abundance
# Y -- sampled abundance

library( RandomFields )
library( TMB )
library( INLA )

Version = "simple_v2"
# v2 -- added option for changing distribution for random effects, and implemented ICAR distribution

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
SD_Nu = 1
betay = 2
Use_REML = FALSE
Spatial_model_set = c("SPDE_GMRF", "ICAR")[2]
Spatial_sim_model = c("GP_gaussian", "ICAR")[2]

# Source ICAR precision matrix function
Q = rect_adj( x=grid_dim['x'], y=grid_dim['y'] )
Q = -1*Q
diag(Q) = -1 * colSums(Q)
# Sigma = solve(Q)   # Appears to be singular

# Calculate C (pg. 185 Cressie and Wikle)
#C = rect_adj( x=grid_dim['x'], y=grid_dim['y'] )

# Results
Results = array(NA, dim=c(length(Spatial_model_set),500,5), dimnames=list(Spatial_model_set,NULL,c("betay","Sigma_Nu","Range_Nu","sum_Ztrue","sum_Zpred")) )

# Loop through replicates
for(i in 1:dim(Results)[2]){
  # Spatial model
  loc_s = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
  model_Nu <- RMgauss(var=SD_Nu^2, scale=SpatialScale)

  # Realization from GRF
  if( Spatial_sim_model=="GP_gaussian") Nu = RFsimulate(model=model_Nu, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
  if( Spatial_sim_model=="ICAR") Nu = rrw( Q )[,1]

  # Linear predictors
  Zy = cbind( rep(1,prod(grid_dim)) )

  # Total abundance
  Ztrue_s = exp(Zy*betay + Nu)

  # Counting process
  s_i = sample(1:prod(grid_dim), size=n_samp, replace=TRUE)
  c_i = rpois( n=n_samp, lambda=Ztrue_s[s_i])
  plot( x=Ztrue_s[s_i], y=c_i)

  # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
  mesh = inla.mesh.create( loc_s )

  # Loop 
  for(model_i in 1:length(Spatial_model_set)){

    # Options
    Spatial_model = Spatial_model_set[model_i]
    Options_vec = c( 'Prior'=switch(Spatial_model,"ICAR"=1,"SPDE_GMRF"=0) )
  
    # Data
    spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
    if(Version%in%c("simple_v1")) Data = list( "c_i"=c_i, "s_i"=s_i-1, "X_sj"=Zy, "spde"=spde)
    if(Version%in%c("simple_v2")) Data = list( "Options_vec"=Options_vec, "c_i"=c_i, "s_i"=s_i-1, "X_sj"=Zy, "spde"=spde, "Q_ICAR"=Matrix(Q))
  
    # Parameters
    if( Options_vec['Prior']==0 ) nuinput_s = rep(0,mesh$n)
    if( Options_vec['Prior']==1 ) nuinput_s = rep(0,prod(grid_dim))
    if(Version %in% c("simple_v2","simple_v1")) Params = list("beta_j"=rep(0,ncol(Data$X_sj)), "logtau"=log(1), "logkappa"=log(1), "nuinput_s"=nuinput_s )
  
    # Random
    Random = c( "nuinput_s" )
    if( Use_REML==TRUE) Random = c(Random,"beta_j")
  
    # Fix parameters
    Map = list()
    if( Options_vec['Prior']==1 ){
      Params[["logtau"]] = log(1)
      Map[["logtau"]] = factor(NA)
      Params[["beta_j"]] = log(1)
      Map[["beta_j"]] = factor(NA)
    }
    
    # Make object
    dyn.load( dynlib(Version) )
    Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE)
    Obj$env$silent=TRUE
  
    # Run
    Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list(trace=1))         #
    Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  
    # SD
    SD = sdreport( Obj, bias.correct=TRUE)
    Report = Obj$report()
    
    # Record
    Results[Spatial_model,i,c("betay","Sigma_Nu","Range_Nu","sum_Ztrue","sum_Zpred")] = c( NA, Report$MargSD, Report$Range, sum(Ztrue_s), SD$unbiased$value['total_abundance'] )
    if( Use_REML==TRUE ) Results[Spatial_model,i,'betay'] = SD$unbiased$value['beta_j']
  }
  
  # output
  write.table( aperm(Results,c(2,3,1)), file="Results.txt")
  save( Results, file="Results.RData")
}

####################
# Read results
####################
TmbFile = "C:/Users/James.Thorson/Desktop/Project_git/pref_sampling/TMB_version/inst/executables"
setwd(TmbFile)
load( "Results.RData")
WhichDone = which( !is.na(Results[1,,'sum_Ztrue']))

# Print to screen
Results[,WhichDone,]
Stats = function(vec){ c("mean"=mean(vec,na.rm=TRUE),"median"=median(vec,na.rm=TRUE),"SE"=sqrt(var(vec,na.rm=TRUE)/sum(!is.na(vec))))}
apply(Results, MARGIN=c(1,3), FUN=Stats)
apply(Results[,,'sum_Zpred',drop=FALSE]/Results[,,'sum_Ztrue',drop=FALSE], MARGIN=1, FUN=Stats)
length(WhichDone)

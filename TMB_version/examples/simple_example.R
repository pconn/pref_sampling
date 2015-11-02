
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

# Settings
grid_dim = c("x"=30, "y"=30)
n_samp = 50
SpatialScale = sqrt(prod(grid_dim))/5  # Range ~ 2*Scale
SD_Nu = 1
betay = 2
Use_REML = TRUE
Spatial_model = c("ICAR","SPDE_GMRF")[1]

# Derived
Options_vec = c( 'Prior'=switch(Spatial_model,"ICAR"=1,"SPDE_GMRF"=0) )

# Source ICAR precision matrix function
source( paste0(TmbFile,"/../../examples/rect_adj.R") )
Q = rect_adj( x=grid_dim['x'], y=grid_dim['y'] )
Q = -1*Q
diag(Q) = -1 * colSums(Q)
# Sigma = solve(Q)   # Appears to be singular

# Results
Results = array(NA, dim=c(500,5), dimnames=list(NULL,c("betay","Sigma_Nu","Range_Nu","sum_Ztrue","sum_Zpred")) )

# Lop through replicates
for(i in 1:nrow(Results)){
  # Spatial model
  loc_s = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
  model_Nu <- RMgauss(var=SD_Nu^2, scale=SpatialScale)

  # Realization from GRF
  Nu = RFsimulate(model=model_Nu, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]

  # Linear predictors
  Zy = cbind( rep(1,prod(grid_dim)) )

  # Total abundance
  Ztrue_s = exp(Zy*betay + Nu)

  # Counting process
  s_i = sample(1:prod(grid_dim), size=n_samp, replace=TRUE)
  c_i = rpois( n=n_samp, lambda=Ztrue_s[s_i])
  plot( x=Ztrue[WhichSample], y=c_i)

  # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
  mesh = inla.mesh.create( loc_s )

  # Data
  spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
  if(Version%in%c("simple_v1")) Data = list( "c_i"=c_i, "s_i"=s_i-1, "X_sj"=Zy, "spde"=spde)

  # Parameters
  if(Version=="simple_v1") Params = list("beta_j"=rep(0,ncol(Data$X_sj)), "logtau"=log(1), "logkappa"=log(1), "nuinput_s"=rep(0,mesh$n) )

  # Random
  Random = c( "nuinput_s" )
  if( Use_REML==TRUE) Random = c(Random,"beta_j")

  # Make object
  dyn.load( dynlib(Version) )
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, silent=TRUE)
  Obj$env$silent=TRUE

  # Run
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list(trace=1))         #

  # SD
  SD = sdreport( Obj, bias.correct=TRUE)
  Report = Obj$report()
  
  # Record
  Results[i,c("betay","Sigma_Nu","Range_Nu","sum_Ztrue","sum_Zpred")] = c( SD$unbiased$value['beta_j'], Report$MargSD, Report$Range, sum(Ztrue_s), SD$unbiased$value['total_abundance'] )

  # output
  write.table( Results, file="Results.txt")
}

# Read results
TmbFile = "C:/Users/James.Thorson/Desktop/Project_git/pref_sampling/TMB_version/inst/executables"
setwd(TmbFile)
Results = read.table( "Results.txt")
WhichDone = which( !is.na(Results[,'sum_Ztrue']))

# Print to screen
Results[WhichDone,]
Stats = function(vec){ c("mean"=mean(vec,na.rm=TRUE),"median"=median(vec,na.rm=TRUE),"SE"=sqrt(var(vec,na.rm=TRUE)/sum(!is.na(vec))))}
apply(Results, MARGIN=2, FUN=Stats)
Stats(Results[,'sum_Zpred']/Results[,'sum_Ztrue'])
length(WhichDone)



#  Y|R=1 ~ Pois( exp(Xy*betay + Nu + Delta) )
#  R ~ Bernoulli( exp(Xr*betar + Gamma + Delta) )
library( RandomFields )
library( TMB )
library( INLA )

# Compile
Version = "pref_samp_v1"
  #TmbFile = system.file("executables", package="SpatialDFA")
  TmbFile = "C:/Users/James.Thorson/Desktop/Project_git/preferential_sampling/inst/executables"
setwd( TmbFile )
compile( paste0(Version,".cpp") )

# Parameters
n_stations = 200
SpatialScale = 0.25
SD_Nu = 0.5
SD_Gamma = 0.5
SD_Delta = 1
betay = 2
betar = 0

###################
# Simulate loop
###################
Results = array(NA, dim=c(100,4), dimnames=list(NULL,c("Sigma_Nu","Sigma_Gamma","Sigma_Delta","Range")) )

for(i in 1:nrow(Results)){
  ###################
  # Simulate data
  ###################

  #Sim_Fn( n_stations, SD_Nu, SD_Gamma, SD_Gamma, SD_Delta, SpatialScale, betay, betar )
  #Return = list("Xr"=Xr, "Xy"=Xy, "R"=R, "Y"=Y, "Loc"=Loc, "Nu"=Nu, "Gamma"=Gamma, "Delta"=Delta)

  # Spatial model
  Loc = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )
  model_Nu <- RMgauss(var=SD_Nu^2, scale=SpatialScale)
  model_Gamma <- RMgauss(var=SD_Gamma^2, scale=SpatialScale)
  model_Delta <- RMgauss(var=SD_Delta^2, scale=SpatialScale)

  # Realization from GRF
  Nu = RFsimulate(model=model_Nu, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
  Gamma = RFsimulate(model=model_Gamma, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
  Delta = RFsimulate(model=model_Delta, x=Loc[,'x'], y=Loc[,'y'])@data[,1]

  # Linear predictors
  Xr = cbind( rep(1,n_stations) )
  Xy = cbind( rep(1,n_stations) )

  # Sampling process
  Rpred = plogis( Xr%*%betar + Gamma + Delta)
  R = rbinom( n=n_stations, size=1, prob=Rpred )

  # Counting process
  Ypred = exp(Xy*betay + Nu + Delta)
  Y = rpois( n=n_stations, lambda=ifelse(R==1, Ypred, NA))


  ###################
  # Fit model
  ###################

  # Build SPDE object using INLA
  mesh = inla.mesh.create( Loc, plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=list(min.angle=26) )

  # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
  spde = inla.spde2.matern(mesh, alpha=2)

  # Data
  if(Version=="pref_samp_v1") Data = list("n_stations"=n_stations, "n_knots"=mesh$n, "ncol_Xr"=ncol(Xr), "ncol_Xy"=ncol(Xy), "R_i"=R, "Y_i"=Y, "Xr"=Xr, "Xy"=Xy, "G0"=spde$param.inla$M0, "G1"=spde$param.inla$M1, "G2"=spde$param.inla$M2 )

  # Parameters
  if(Version=="pref_samp_v1") Params = list("log_kappa"=log(1), "logtau_Nu"=log(1), "logtau_Gamma"=log(1), "logtau_Delta"=log(1), "betar"=rep(0,ncol(Xr)), "betay"=rep(0,ncol(Xy)), "Nu_input"=rep(0,Data$n_knots), "Gamma_input"=rep(0,Data$n_knots), "Delta_input"=rep(0,Data$n_knots) )

  # Random
  Random = c( "Nu_input", "Gamma_input", "Delta_input" )
  #if(Use_REML==TRUE) ...

  # Fixed values
  Map = list()

  # Load DLL
  setwd( TmbFile )
  dyn.load( dynlib(Version) )                                                         # log_tau=0.0,

  # Initialization
  obj <- MakeADFun(data=Data, parameters=Params, random=Random, map=Map, hessian=FALSE, inner.control=list(maxit=1000) )
  obj$control <- c( obj$control, list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100) )
  obj$env$inner.control <- c(obj$env$inner.control, list("step.tol"=1e-8, "tol10"=1e-6, "grad.tol"=1e-8) )

  # First run
  Init = obj$fn( obj$par )
  Initial_gradient = obj$gr( obj$par )

  # Bounds
  Upper = rep(Inf, length(obj$par) )
  Lower = rep(-Inf, length(obj$par) )

  # Run model
  opt = nlminb(start=obj$env$last.par.best[-c(obj$env$random)], objective=obj$fn, gradient=obj$gr, upper=Upper, lower=Lower, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )
  opt[["final_gradient"]] = obj$gr( opt$par )
  opt[["AIC"]] = 2*opt$objective + 2*length(opt$par)
  opt[["BIC"]] = 2*opt$objective + length(opt$par) * log(n_stations)

  # Diagnostics
  Report = obj$report()
  Results[i,] = unlist( Report[c("Sigma_Nu","Sigma_Gamma","Sigma_Delta","Range")] )
}


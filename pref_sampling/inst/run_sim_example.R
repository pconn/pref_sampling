### run simple preferential sampling example.  Doing in this in two steps since there seems to be
# a conflict between TMB loop and rgeos ('gCentroid') functionality
# need to declare an ADREPORT for Z_s in TMB for this
library(PrefSampling)

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
RandomSeed = ceiling(runif(1,min=1,max=1e6))

# Configurations
b_set = c(0,1,3)            # Impact of delta on sampling intensity
EM_set = c("fix_b","est_b")

# Configurations

set.seed(2222222)

Counts=Est=vector("list",3)


# Settings
betax = runif(1,-.5,.5)       # Impact of x on density
betax_prob = runif(1,-.5,.5)   # Impact of x on sampling intensity
#betax=0
#betax_prob=0

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
 
#loop over sampling, estimation
for(SimI in 1:length(b_set)){
  b = b_set[SimI]
  R_s = exp( betax_prob*x_s + eta_s + b*delta_s )
  R_s = R_s / sum(R_s)
  
  # Process for locating samples
  s_i = sample(1:prod(grid_dim), size=n_samp, replace=FALSE, prob=R_s)
  y_s = ifelse(1:prod(grid_dim) %in% s_i, 1, 0)
  
  # Counting process
  #c_i = rpois( n=n_samp, lambda=Ztrue_s[s_i])
  c_i = rbinom(n=n_samp,Ztrue_s[s_i],prop_sampled)   #changed to binom for v5, 11/9
  
  Cur.count=rep(NA,prod(grid_dim))
  Cur.count[s_i]=c_i
  Counts[[SimI]]=Cur.count
  
  for(EstI in 1:1){
    EM = EM_set[EstI]
    
    # Create the SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
    mesh = inla.mesh.create( loc_s )
    
    # Loop 
    # Options
    Options_vec = c( 'Prior'=switch(Spatial_model,"ICAR"=1,"SPDE_GMRF"=0), 'Alpha'=Alpha, 'IncludeDelta'=1, 'IncludeEta'=1, 'OutputSE'=1)
    
    # Data
    spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
    Data = list( "Options_vec"=Options_vec, "c_i"=c_i, "P_i"=Prop_sampled,"A_s"=rep(1,n_cells),"s_i"=s_i-1, "X_sj"=cbind(1,x_s), "y_s"=y_s, "X_sk"=cbind(x_s), "spde"=spde, "X_sb"=matrix(1,nrow=n_cells))
    
    # Parameters
    # Intercept is needed for beta_j (delta -- abundance) but not beta_k (eta -- sampling intensity)
    if( Options_vec['Prior']==0 ) etainput_s = deltainput_s = rep(0,mesh$n)
    if( Options_vec['Prior']==1 ) etainput_s = deltainput_s = rep(0,prod(grid_dim))
    Params = list("beta_j"=rep(0,ncol(Data$X_sj)), "beta_k"=rep(0,ncol(Data$X_sk)), "beta_b"=0, "logtau_z"=rep(0,2), "logkappa_z"=rep(0,2), "deltainput_s"=deltainput_s, "etainput_s"=etainput_s)
    Params$beta_j[1]=median(log(Ztrue_s+0.01)) #set mean expected abundance close to truth for faster optimization
    
    # Random
    Random = c( "deltainput_s", "etainput_s" )
    if(Use_REML==TRUE) Random = c(Random,"beta_j","beta_k")
    if(Use_REML==TRUE) Random = c( Random, "b" )
    
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
    #dyn.load( dynlib(Version) )
    Start_time = Sys.time()
    Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, silent=TRUE, DLL="PrefSampling")
    Obj$fn( Obj$par )
    
    # Run
    Lower = -50
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
        Map[["b"]] = factor(NA)
        Params[["b"]] = 0
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
    
    # SD
    #Report = Obj$report()
    SD=sdreport(Obj,bias.correct=TRUE)
    Est[[SimI]]=SD$unbiased$value[6:630]
    #if( all(c("etainput_s","deltainput_s")%in%names(Map)) ){
    #  Est[[EstI]]=Report$Z_s
    #}else{
    #  SD = sdreport( Obj, bias.correct=TRUE )  
    #  Est[[EstI]]=SD$unbiased$Z_s
    #}
  }     
}

Sim=list("N.1"=Est[[1]],"N.2"=Est[[2]],"N.3"=Est[[3]],"C.1"=Counts[[1]],"C.2"=Counts[[2]],"C.3"=Counts[[3]],"N.true"=Ztrue_s,"Delta"=delta_s,"Cov"=x_s)
save(Sim,file="sim_plot_data.Rdata")


// Space time 
#include <TMB.hpp>
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type plogis(Type x){
  return 1.0 / (1.0 + exp(-x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Indices
  DATA_INTEGER(n_cells);       // Total number of areal survey units in study area
  DATA_INTEGER(ncol_Xr);       // Total number of observations (i)
  DATA_INTEGER(ncol_Xy);       // Total number of observations (i)

  // Data
  DATA_FACTOR( R_i );         // Count for observation
  DATA_VECTOR( Y_i );       	// Species for observation
  DATA_VECTOR( Log_area_i );        // A vector giving the (standardized) area of suitable habitat associated with each cell (length = n_cells)
  DATA_VECTOR( Prop_sampled_i);  // Proportion of suitable habitat in each cell that is sampled (length = n_cells) 

  // Covariates
  DATA_MATRIX( Xr );  //design matrix for Bernoulli sampling model
  DATA_MATRIX( Xy );  //design matrix for Poisson count model
  
  // Spatial precision matrix (Q where Sigma^{-1} = \tau Q)  Note: not all of these are used
  DATA_SPARSE_MATRIX(Q1);
  DATA_SPARSE_MATRIX(Q2);
  DATA_SPARSE_MATRIX(Q3);

 
  // Fixed effects
  PARAMETER(logtau_Eta);           // log-inverse SD of Eta
  PARAMETER(logtau_Gamma);        // log-inverse SD of Gamma
  PARAMETER(beta_pref);    //preferential sampling effect
  //PARAMETER(logsigma_Epsilon);   
  PARAMETER_VECTOR(betar);
  PARAMETER_VECTOR(betay);

  // Random effects
  PARAMETER_VECTOR(Eta);  // Poisson count spatial random effect
  PARAMETER_VECTOR(Gamma);   // Bernoulli site selection spatial random effect
  //PARAMETER_VECTOR(Epsilon);  //residual error random effect on linear predictor for count model
  
  

  // global stuff
  using namespace density; //gives access to multivariate distributions in density.cpp
  Type jnll = 0;
  
  // Derived quantities related to GMRF.  NOTE: No residual white noise assumed
  Eigen::SparseMatrix<Type> Prec;   //precision matrix for GMRF

  //Probability of random fields
  Prec = exp(logtau_Eta)*Q1;
  //Type sigma=exp(logsigma_Epsilon);
  jnll += GMRF(Prec)(Eta);  //no residual error in here yet
  Prec = exp(logtau_Gamma)*Q2;
  jnll += GMRF(Prec)(Gamma);  //no residual error in here yet
  //for(int n=0; n<n_cells; n++){
  //  jnll -= dnorm(Epsilon(n),Type(0.0),sigma,1);
  //}  
  
  //Type pi = 3.141592;
  //vector<Type> Range(3);
  //for(int r=0; r<Range.size(); r++) Range(r) = sqrt(8) / exp( log_kappa(r) );
  //Type Sigma_Nu = 1 / sqrt(4*pi*exp(2*logtau_Nu)*exp(2*log_kappa(0)));
  //Type Sigma_Gamma = 1 / sqrt(4*pi*exp(2*logtau_Gamma)*exp(2*log_kappa(1)));
  //Type Sigma_Delta = 1 / sqrt(4*pi*exp(2*logtau_Delta)*exp(2*log_kappa(2)));

  // Probability of random fields
  //Prec =  = exp(log_kappa(0))*G0 + Type(2.0)*exp(2.0*log_kappa(0))*G1 + G2;
  //jnll += GMRF(Q)(Nu_input);
  //Q = exp(4.0*log_kappa(1))*G0 + Type(2.0)*exp(2.0*log_kappa(1))*G1 + G2;
  //jnll += GMRF(Q)(Gamma_input);
  //Q = exp(4.0*log_kappa(2))*G0 + Type(2.0)*exp(2.0*log_kappa(2))*G1 + G2;
  //jnll += GMRF(Q)(Delta_input);

  // Derived quantities
  vector<Type> Rpred_i(n_cells);
  vector<Type> Ypred_i(n_cells);
  Type tot_abund=0;
  
  // Probability of observations
  for(int i=0; i<n_cells; i++){
    // probability of counts
    Ypred_i(i) = Log_area_i(i)+Eta(i); // + Delta(i); // + Epsilon(i);
    for(int c=0; c<ncol_Xy; c++) Ypred_i(i) += Xy(i,c) * betay(c) ;
    Rpred_i(i) = Gamma(i) + beta_pref*Ypred_i(i);
    Ypred_i(i) = exp(Ypred_i(i));
    if( !isNA(Y_i(i)) ) jnll -= dpois( Y_i(i), Prop_sampled_i(i)*Ypred_i(i), true );
    tot_abund += Ypred_i(i);
    // probability of sampling
     for(int c=0; c<ncol_Xr; c++) Rpred_i(i) += Xr(i,c) * betar(c);
    Rpred_i(i) = plogis(Rpred_i(i));
    if(R_i(i)==0) jnll -= log( 1-Rpred_i(i) );
    if(R_i(i)==1) jnll -= log( Rpred_i(i) );
  }
  //std::cout<<jnll;
  

  // Spatial field summaries
  REPORT( Rpred_i );
  REPORT( Ypred_i );
  ADREPORT( tot_abund);

  return jnll;
}

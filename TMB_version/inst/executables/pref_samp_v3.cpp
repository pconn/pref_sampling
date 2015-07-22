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
  
  // ICAR structure
  DATA_SPARSE_MATRIX(Q);
 
  // Fixed effects
  PARAMETER(logtau_Eta);           // log-inverse SD of Eta
  PARAMETER(logtau_Gamma);        // log-inverse SD of Gamma
  PARAMETER(logtau_Delta);        // log-inverse SD of Delta
  //PARAMETER(logsigma_Epsilon);   
  PARAMETER_VECTOR(betar);
  PARAMETER_VECTOR(betay);

  // Random effects
  PARAMETER_VECTOR(Eta);  // Poisson count spatial random effect
  PARAMETER_VECTOR(Gamma);   // Bernoulli site selection spatial random effect
  PARAMETER_VECTOR(Delta);  // shared spatial random effect
  //PARAMETER_VECTOR(Epsilon);  //residual error random effect on linear predictor for count model
  
  

  // global stuff
  using namespace density; //gives access to multivariate distributions in density.cpp
  Type jnll = 0;
  
  
  //random effects priors
  Type tau_eta = exp(logtau_Eta);
  vector<Type> Tmp = Q*Eta;
  jnll -= Type(0.5)*(n_cells*logtau_Eta-tau_eta*(Eta*Tmp).sum());
  
  Type tau_gamma = exp(logtau_Gamma);
  Tmp = Q*Gamma;
  jnll -= Type(0.5)*(n_cells*logtau_Gamma-tau_gamma*(Gamma*Tmp).sum());
  
  Type tau_delta = exp(logtau_Delta);
  Tmp = Q*Delta;
  jnll -= Type(0.5)*(n_cells*logtau_Delta-tau_delta*(Delta*Tmp).sum());
  
  //prior on tau
  jnll -= dgamma(tau_eta,Type(1.0),Type(0.01), true ); 
  jnll -= dgamma(tau_gamma,Type(1.0),Type(0.01), true ); 
  jnll -= dgamma(tau_delta,Type(1.0),Type(0.01), true ); 

  // Derived quantities
  vector<Type> Rpred_i(n_cells);
  vector<Type> Ypred_i(n_cells);
  Type tot_abund=0;
  
  // Probability of observations
  for(int i=0; i<n_cells; i++){
    // probability of sampling
    Rpred_i(i) = 0 + Gamma(i) + Delta(i);
    for(int c=0; c<ncol_Xr; c++) Rpred_i(i) += Xr(i,c) * betar(c);
    Rpred_i(i) = plogis(Rpred_i(i));
    if(R_i(i)==0) jnll -= log( 1-Rpred_i(i) );
    if(R_i(i)==1) jnll -= log( Rpred_i(i) );
    // probability of counts
    Ypred_i(i) = Log_area_i(i)+Eta(i) + Delta(i); // + Epsilon(i);
    for(int c=0; c<ncol_Xy; c++) Ypred_i(i) += Xy(i,c) * betay(c) ;
    Ypred_i(i) = exp(Ypred_i(i));
    if( !isNA(Y_i(i)) ) jnll -= dpois( Y_i(i), Prop_sampled_i(i)*Ypred_i(i), true );
    tot_abund += Ypred_i(i);
  }
  //std::cout<<jnll;
  

  // Spatial field summaries
  REPORT( Rpred_i );
  REPORT( Ypred_i );
  ADREPORT( tot_abund);

  return jnll;
}

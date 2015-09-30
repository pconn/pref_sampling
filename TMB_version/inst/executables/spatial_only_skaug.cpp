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
  DATA_INTEGER(ncol_X);       // Number of covariates

  // Data
  DATA_VECTOR( Y_i );       	// Counts
  DATA_VECTOR( Log_area_i );        // A vector giving the (standardized) area of suitable habitat associated with each cell (length = n_cells)
  DATA_VECTOR( Prop_sampled_i);  // Proportion of suitable habitat in each cell that is sampled (length = n_cells) 

  // Covariates
  DATA_MATRIX(X);  //design matrix for Poisson count model
  
  // ICAR structure
  DATA_SPARSE_MATRIX(Q);
 
  // Fixed effects
  PARAMETER(logtau_Eta);           // log-inverse SD of Eta

  //PARAMETER_VECTOR(Beta); //(run without since ICAR rank deficient by 1)

  // Random effects
  PARAMETER_VECTOR(Eta);  // random effects  

  // global stuff
  using namespace density; //gives access to multivariate distributions in density.cpp
  Type jnll = 0;
  
   //
  //Probability of random fields
  // ICAR prior
  Type tau = exp(logtau_Eta);
  vector<Type> Tmp = Q*Eta;
  jnll -= Type(0.5)*((n_cells-1)*logtau_Eta-tau*(Eta*Tmp).sum());

  //vector<Type> Eta = Tmp/tau;       // Doing the scaling as a post-step. 100% OK!!

  //prior on tau
  jnll -= dgamma(tau,Type(1.0),Type(0.01), true ); 

  // Derived quantities
  vector<Type> Ypred_i(n_cells);
  Type tot_abund=0;
  
  // Probability of observations
  for(int i=0; i<n_cells; i++){
    // probability of counts
    Ypred_i(i) = Log_area_i(i)+Eta(i);
    //for(int c=0; c<ncol_X; c++) Ypred_i(i) += X(i,c) * Beta(c) ; //run without intercept temporarily
    Ypred_i(i) = exp(Ypred_i(i));
    if( !isNA(Y_i(i)) ) jnll -= dpois( Y_i(i), Prop_sampled_i(i)*Ypred_i(i), true );
    tot_abund += Ypred_i(i);
  }
  //std::cout<<jnll;
  

  // Spatial field summaries
  REPORT( Ypred_i );
  ADREPORT( tot_abund);

  return jnll;
}

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
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;
  
  // options vec
  DATA_FACTOR( Options_vec );
  // Slot 0: prior on random effects (0=SPDE_GMRF; 1=ICAR_GMRF)
  
  // Data
  DATA_VECTOR( c_i );       	// Response (count) for each observation
  DATA_FACTOR( s_i ); // Site for each sample
  DATA_MATRIX( X_sj );  //design matrix for Bernoulli sampling model
  
  // Aniso objects
  DATA_STRUCT(spde,spde_t);
  DATA_SPARSE_MATRIX(Q_ICAR);
 
  // Parameters 
  PARAMETER_VECTOR(beta_j);  
  PARAMETER(logtau);      
  PARAMETER(logkappa)
  
  // Random effects
  PARAMETER_VECTOR( nuinput_s );
  
  // derived sizes
  int n_i = c_i.size();
  int n_j = X_sj.row(0).size();
  int n_s = X_sj.col(0).size();
  int n_k = nuinput_s.size();
  
  // global stuff
  Type jnll = 0;
  Type MargSD = 1 / sqrt(4*M_PI) / exp(logtau) / exp(logkappa);
  Type Range = sqrt(8) / exp( logkappa );
    
  //random effects priors
  if( Options_vec(0)==0 ){
    SparseMatrix<Type> Q = Q_spde(spde, exp(logkappa));
    jnll += GMRF(Q)(nuinput_s);
  }
  if( Options_vec(0)==1 ){
    Type tau = exp(logkappa);
    vector<Type> Tmp = Q_ICAR*nuinput_s;
    jnll -= Type(0.5) * ( (n_k-Type(1.0) ) * log(tau) - tau * (nuinput_s*Tmp).sum() );
  }
  
  // Predicted densities
  vector<Type> Z_s(n_s);
  for(int s=0; s<n_s; s++){
    Z_s(s) = exp( ( X_sj.row(s).array() * beta_j ).sum() + nuinput_s(s)/exp(logtau) );
  }
  Type total_abundance = Z_s.sum();
  
  // Probability of counts
  for(int i=0; i<n_i; i++){
    if(!isNA(c_i(i))) jnll -= dpois( c_i(i), Z_s(s_i(i)), true );
  }

  // Reporting
  REPORT( Z_s );
  REPORT( total_abundance );
  REPORT( Range );
  REPORT( MargSD );
  REPORT( beta_j );
  
  // Bias correction output
  ADREPORT( beta_j );
  ADREPORT( total_abundance);

  return jnll;
}

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
Type dbern(Type x, Type prob, int give_log=1){
  Type logres;
  if( x==0 ) logres = log( 1-prob );
  if( x==1 ) logres = log( prob );
  if(give_log) return logres; else return exp(logres);
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
  DATA_VECTOR( c_i );       	// Response (count) for each observation i
  DATA_FACTOR( s_i ); // Site for each sample
  DATA_MATRIX( X_sj );  //design matrix for Observation model
  DATA_VECTOR( y_s );       	// Response (0:not surveyed, 1:surveyed) for each site
  DATA_MATRIX( X_sk );  //design matrix for Bernoulli sampling model
  
  // Aniso objects
  DATA_STRUCT(spde,spde_t);
  //DATA_SPARSE_MATRIX(Q_ICAR);
 
  // Parameters 
  PARAMETER_VECTOR(beta_j);  
  PARAMETER_VECTOR(beta_k);  
  PARAMETER( b );
  PARAMETER_VECTOR(logtau_z);      
  PARAMETER_VECTOR(logkappa_z)
  
  // Random effects
  PARAMETER_VECTOR( deltainput_s );
  PARAMETER_VECTOR( etainput_s );
  
  // derived sizes
  int n_i = c_i.size();
  int n_j = X_sj.row(0).size();
  int n_k = X_sk.row(0).size();
  int n_s = X_sj.col(0).size();
  int n_z = logkappa_z.size();  

  // global stuff
  Type jnll = 0;
  vector<Type> MargSD_z(n_z);
  vector<Type> Range_z(n_z);
  for( int z=0; z<n_z; z++){
    MargSD_z(z) = 1 / sqrt(4*M_PI) / exp(logtau_z(z)) / exp(logkappa_z(z));
    Range_z(z) = sqrt(8) / exp( logkappa_z(z) );
  }
    
  //random effects priors
  if( Options_vec(0)==0 ){
    SparseMatrix<Type> Q;
    Q = Q_spde(spde, exp(logkappa_z(0)));
    jnll += GMRF(Q)(deltainput_s);
    Q = Q_spde(spde, exp(logkappa_z(1)));
    jnll += GMRF(Q)(etainput_s);
  }
  if( Options_vec(0)==1 ){
    //vector<Type> tau_z(n_z);
    //vector<Type> Tmp(n_s);
    //tau_z = exp(logkappa_z);
    //Tmp = Q_ICAR * deltainput_s;
    //jnll -= Type(0.5) * ( (n_s-Type(1.0) ) * log(tau_z(0)) - tau_z(0) * (deltainput_s*Tmp).sum() );
    //Tmp = Q_ICAR * etainput_s;
    //jnll -= Type(0.5) * ( (n_s-Type(1.0) ) * log(tau_z(1)) - tau_z(1) * (etainput_s*Tmp).sum() );
  }
  
  // Predicted densities
  vector<Type> Z_s(n_s);
  vector<Type> debug_j(n_j);
  vector<Type> linpred_s = X_sj * beta_j;
  for(int s=0; s<n_s; s++){
    Z_s(s) = exp( linpred_s(s) + deltainput_s(s)/exp(logtau_z(0)) );
  }
  Type total_abundance = Z_s.sum();
  
  // Probability of counts
  for(int i=0; i<n_i; i++){
    if( !isNA(c_i(i)) ) jnll -= dpois( c_i(i), Z_s(s_i(i)), true );
  }

  // Predicted probability of sampling
  vector<Type> R_s(n_s);
  linpred_s = X_sk * beta_k;
  for(int s=0; s<n_s; s++){
    R_s(s) = exp( linpred_s(s) + b*deltainput_s(s)/exp(logtau_z(0)) + etainput_s(s)/exp(logtau_z(1)) );
  }
  R_s = R_s / R_s.sum();
  
  // Probability of sampling locations
  for(int s=0; s<n_s; s++){
    jnll -= dbern( y_s(s), R_s(s), true );
  }

  // Reporting
  REPORT( Z_s );
  REPORT( R_s );
  REPORT( total_abundance );
  REPORT( Range_z );
  REPORT( MargSD_z );
  REPORT( beta_j );
  REPORT( beta_k );
  REPORT( b );
  
  // Bias correction output
  ADREPORT( beta_j );
  ADREPORT( total_abundance);

  return jnll;
}

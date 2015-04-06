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
  DATA_INTEGER(n_stations);       // Total number of observations (i)
  DATA_INTEGER(n_knots);       // Total number of observations (i)
  DATA_INTEGER(ncol_Xr);       // Total number of observations (i)
  DATA_INTEGER(ncol_Xy);       // Total number of observations (i)

  // Data
  DATA_FACTOR( R_i );         // Count for observation
  DATA_VECTOR( Y_i );       	// Species for observation

  // Covariates
  DATA_MATRIX( Xr );
  DATA_MATRIX( Xy );
  
  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER_VECTOR(log_kappa);         // Controls range of spatial variation
  PARAMETER(logtau_Nu);        // log-inverse SD of Epsilon
  PARAMETER(logtau_Gamma);        // log-inverse SD of Epsilon
  PARAMETER(logtau_Delta);        // log-inverse SD of Epsilon
  PARAMETER_VECTOR(betar);
  PARAMETER_VECTOR(betay);

  // Random effects
  PARAMETER_VECTOR(Nu_input);  // Spatial process variation
  PARAMETER_VECTOR(Gamma_input);   // Spatial variation in carrying capacity
  PARAMETER_VECTOR(Delta_input);

  // global stuff
  using namespace density;
  Type jnll = 0;
  
  // Derived quantities related to GMRF
  Type pi = 3.141592;
  vector<Type> Range(3);
  for(int r=0; r<Range.size(); r++) Range(r) = sqrt(8) / exp( log_kappa(r) );
  Type Sigma_Nu = 1 / sqrt(4*pi*exp(2*logtau_Nu)*exp(2*log_kappa(0)));
  Type Sigma_Gamma = 1 / sqrt(4*pi*exp(2*logtau_Gamma)*exp(2*log_kappa(1)));
  Type Sigma_Delta = 1 / sqrt(4*pi*exp(2*logtau_Delta)*exp(2*log_kappa(2)));

  // Probability of random fields
  Eigen::SparseMatrix<Type> Q;
  Q = exp(4.0*log_kappa(0))*G0 + Type(2.0)*exp(2.0*log_kappa(0))*G1 + G2;
  jnll += GMRF(Q)(Nu_input);
  Q = exp(4.0*log_kappa(1))*G0 + Type(2.0)*exp(2.0*log_kappa(1))*G1 + G2;
  jnll += GMRF(Q)(Gamma_input);
  Q = exp(4.0*log_kappa(2))*G0 + Type(2.0)*exp(2.0*log_kappa(2))*G1 + G2;
  jnll += GMRF(Q)(Delta_input);

  // Transform random fields
  vector<Type> Nu(n_knots);
  vector<Type> Gamma(n_knots);
  vector<Type> Delta(n_knots);
  for(int n=0; n<n_knots; n++){
    Nu(n) = Nu_input(n) / exp(logtau_Nu);
    Gamma(n) = Gamma_input(n) / exp(logtau_Gamma);
    Delta(n) = Delta_input(n) / exp(logtau_Delta);
  }
  
  // Derived quantities
  vector<Type> Rpred_i(n_stations);
  vector<Type> Ypred_i(n_stations);
  
  // Probability of observations
  for(int i=0; i<n_stations; i++){
    // probability of sampling
    Rpred_i(i) = Gamma(i) + Delta(i);
    for(int c=0; c<ncol_Xr; c++) Rpred_i(i) += Xr(i,c) * betar(c);
    Rpred_i(i) = plogis(Rpred_i(i));
    if(R_i(i)==0) jnll -= log( 1-Rpred_i(i) );
    if(R_i(i)==1) jnll -= log( Rpred_i(i) );
    // probability of counts
    Ypred_i(i) = Nu(i) + Delta(i);
    for(int c=0; c<ncol_Xy; c++) Ypred_i(i) += Xy(i,c) * betay(c) ;
    Ypred_i(i) = exp(Ypred_i(i));
    if( !isNA(Y_i(i)) ) jnll -= dpois( Y_i(i), Ypred_i(i), true );
  }

  // Spatial field summaries
  REPORT( Range );
  REPORT( Sigma_Nu );
  REPORT( Sigma_Gamma );
  REPORT( Sigma_Delta );
  REPORT( Rpred_i );
  REPORT( Ypred_i );

  return jnll;
}

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
  
  // Spatial precision matrix (Q where Sigma^{-1} = \tau Q)
  DATA_SPARSE_MATRIX(Q);
 
  // Fixed effects
  PARAMETER(logtau_Eta);           // log-inverse SD of Eta
  PARAMETER(log_rho);

  PARAMETER_VECTOR(Beta);

  // Random effects
  PARAMETER_VECTOR(Eta);  // Poisson count spatial random effect  

  // global stuff
  using namespace density; //gives access to multivariate distributions in density.cpp
  Type jnll = 0;
  
  // Derived quantities related to GMRF.  NOTE: No residual white noise assumed
  Eigen::SparseMatrix<Type> Prec;   //precision matrix for GMRF
  //matrix<Type> Prec(n_cells,n_cells);
  //Eigen::SparseMatrix<Type> Q_rho;
  
  //Q_rho=Q;
  Type rho=exp(log_rho);
  Type tau=exp(logtau_Eta);
  //
  //Probability of random fields
  Prec=Q;
  for(int i=0;i<n_cells;i++)Prec.coeffRef(i,i)=Prec.coeffRef(i,i)+rho;
  Prec = exp(logtau_Eta)*Prec; //Q_rho;


  //Type sigma=exp(logsigma_Epsilon);
  jnll += GMRF(Prec)(Eta);  //no residual error in here yet
  
  //prior for tau
  jnll -= dgamma(tau,Type(1.0),Type(0.01), true ); 
  for(int i=0;i<ncol_X;i++)jnll -= dnorm(Beta[i],Type(0.0),Type(100.0),true);
  
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
  vector<Type> Ypred_i(n_cells);
  Type tot_abund=0;
  
  // Probability of observations
  for(int i=0; i<n_cells; i++){
    // probability of counts
    Ypred_i(i) = Log_area_i(i)+Eta(i);
    for(int c=0; c<ncol_X; c++) Ypred_i(i) += X(i,c) * Beta(c) ;
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

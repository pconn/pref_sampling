// Space time 
#include <TMB.hpp>
//#include <atomic_math.hpp>  //for D_lgamma


/** Precision matrix for the anisotropic case, eqn (20) in Lindgren et al. (2011) */    
namespace R_inla_generalized {
using namespace Eigen;
using namespace tmbutils;
using namespace R_inla;

template<class Type>
  SparseMatrix<Type> Q_spde_generalized(spde_t<Type> spde, Type kappa, int alpha=2){
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  	
  if( alpha==1 ) return kappa_pow2*spde.M0 + spde.M1;
  if( alpha==2 ) return kappa_pow4*spde.M0 + Type(2.0)*kappa_pow2*spde.M1 + spde.M2;
}

template<class Type>
  SparseMatrix<Type> Q_spde_generalized(spde_aniso_t<Type> spde, Type kappa, matrix<Type> H, int alpha=2){

  int i;
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  
  int n_s = spde.n_s;
  int n_tri = spde.n_tri;
  vector<Type> Tri_Area = spde.Tri_Area;
  matrix<Type> E0 = spde.E0;
  matrix<Type> E1 = spde.E1;
  matrix<Type> E2 = spde.E2;
  matrix<int> TV = spde.TV;
  SparseMatrix<Type> G0 = spde.G0;
  SparseMatrix<Type> G0_inv = spde.G0_inv;
	  	  
  //Type H_trace = H(0,0)+H(1,1);
  //Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);
  SparseMatrix<Type> G1_aniso(n_s,n_s); 
  SparseMatrix<Type> G2_aniso(n_s,n_s); 
  // Calculate adjugate of H
  matrix<Type> adj_H(2,2);
  adj_H(0,0) = H(1,1);
  adj_H(0,1) = -1 * H(0,1);
  adj_H(1,0) = -1 * H(1,0);
  adj_H(1,1) = H(0,0);
  // Calculate new SPDE matrices

  // Calculate G1 - pt. 1
  array<Type> Gtmp(n_tri,3,3);
  for(i=0; i<n_tri; i++){    
    // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.    
    Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
    Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
  }
  // Calculate G1 - pt. 2
  for(i=0; i<n_tri; i++){
    G1_aniso.coeffRef(TV(i,1),TV(i,0)) = G1_aniso.coeffRef(TV(i,1),TV(i,0)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,0),TV(i,1)) = G1_aniso.coeffRef(TV(i,0),TV(i,1)) + (Gtmp(i,0,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,1)) = G1_aniso.coeffRef(TV(i,2),TV(i,1)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,1),TV(i,2)) = G1_aniso.coeffRef(TV(i,1),TV(i,2)) + (Gtmp(i,1,2));  
    G1_aniso.coeffRef(TV(i,2),TV(i,0)) = G1_aniso.coeffRef(TV(i,2),TV(i,0)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,2)) = G1_aniso.coeffRef(TV(i,0),TV(i,2)) + (Gtmp(i,0,2));  
    G1_aniso.coeffRef(TV(i,0),TV(i,0)) = G1_aniso.coeffRef(TV(i,0),TV(i,0)) + (Gtmp(i,0,0));  
    G1_aniso.coeffRef(TV(i,1),TV(i,1)) = G1_aniso.coeffRef(TV(i,1),TV(i,1)) + (Gtmp(i,1,1));  
    G1_aniso.coeffRef(TV(i,2),TV(i,2)) = G1_aniso.coeffRef(TV(i,2),TV(i,2)) + (Gtmp(i,2,2));  
  }
  G2_aniso = G1_aniso * G0_inv * G1_aniso; 

  if( alpha==1 ) return kappa_pow2*G0 + G1_aniso;
  if( alpha==2 ) return kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1_aniso + G2_aniso;
}
} // end namespace R_inla

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
  using namespace R_inla_generalized;
  using namespace Eigen;
  using namespace density;
  
  // options vec
  DATA_FACTOR( Options_vec );
  // Slot 0: prior on random effects (0=SPDE_GMRF; 1=ICAR_GMRF)
  // Slot 1: Alpha
  // Slot 2: Include hyperdistribution for delta
  // Slot 3: Include hyperdistribution for eta
  // Slot 4: Output Z_s in ADREPORT (0=no, 1=yes); needed for unbiased map estimates using bias.correct
  // Slot 5: Indicator for whether to put a prior on preferential sampling parameters
  // Data
  DATA_VECTOR( c_i );       	// Response (count) for each observation i
  DATA_VECTOR( P_i );      // Proportion of survey unit that is surveyed for each observation i
  DATA_VECTOR( A_s );      // Relative area of each survey unit s (can set = 1.0 if all the same size)
  DATA_FACTOR( s_i ); // Site for each sample
  DATA_MATRIX( X_sj );  //design matrix for Observation model
  DATA_VECTOR( y_s );       	// Response (0:not surveyed, 1:surveyed) for each site
  DATA_MATRIX( X_sk );  //design matrix for Bernoulli sampling model
  DATA_MATRIX( X_sb );  //design matrix for preferential sampling model
  
  // Aniso objects
  DATA_STRUCT(spde, spde_t);
  //DATA_SPARSE_MATRIX(Q_ICAR);
 
  // Parameters 
  PARAMETER_VECTOR(beta_j);              // Covariates affecting density
  PARAMETER_VECTOR(beta_k);              // Coefficients governing probability of sampling
  PARAMETER_VECTOR(beta_b );             // Preferential sampling parameters, linking density submodel to sampling probability submodel
  PARAMETER_VECTOR(logtau_z);      
  PARAMETER_VECTOR(logkappa_z);
  
  // Random effects
  PARAMETER_VECTOR( deltainput_s );         // Spatial variation in  density Z_s, and sampling probability R_s conditional on preferential sampling parameters
  PARAMETER_VECTOR( etainput_s );           // Spatial variation in sampling probability R_s in excess of preferential sampling predictions

  // derived sizes
  int n_i = c_i.size();
  int n_j = X_sj.row(0).size();
  int n_k = X_sk.row(0).size();
  int n_b = X_sb.row(0).size();
  int n_s = X_sj.col(0).size();
  int n_z = logkappa_z.size();  

  // global stuff
  vector<Type> jnll_comp(4);
  jnll_comp.setZero();
  vector<Type> MargSD_z(n_z);
  vector<Type> Range_z(n_z);
  for( int z=0; z<n_z; z++){
    MargSD_z(z) = 1 / sqrt(4*M_PI) / exp(logtau_z(z)) / exp(logkappa_z(z));
    Range_z(z) = sqrt(8) / exp( logkappa_z(z) );
  }
  vector<Type> Thin_i(n_i);  
  vector<Type> Pc_i(n_i);
  vector<Type> ThinU_s(n_s); 
  for( int s=0;s<n_s;s++){
    ThinU_s(s)=A_s(s);
  }
  for( int i=0;i<n_i;i++){
    Pc_i(i)=1-P_i(i);
    Thin_i(i)=P_i(i)*A_s(s_i(i));
    ThinU_s(s_i(i))=ThinU_s(s_i(i))*Pc_i(i);
  }
  Type c_sum = c_i.sum();
  
  // Transform random effects
  vector<Type> delta_s( deltainput_s.size() );
  vector<Type> eta_s( etainput_s.size() );
  delta_s = deltainput_s / exp(logtau_z(0));
  eta_s = etainput_s / exp(logtau_z(1));

  //random effects priors
  if( Options_vec(0)==0 ){
    SparseMatrix<Type> Q;
    Q = Q_spde_generalized(spde, exp(logkappa_z(0)), Options_vec(1));
    if(Options_vec(2)==1) jnll_comp(2) = GMRF(Q)(deltainput_s);
    Q = Q_spde_generalized(spde, exp(logkappa_z(1)), Options_vec(1));
    if(Options_vec(3)==1) jnll_comp(3) = GMRF(Q)(etainput_s);
  }
  if( Options_vec(0)==1 ){
    //vector<Type> tau_z(n_z);
    //vector<Type> Tmp(n_s);
    //tau_z = exp(logkappa_z);
    //Tmp = Q_ICAR * deltainput_s;
    //if(Options_vec(2)==1) jnll -= Type(0.5) * ( (n_s-Type(1.0) ) * log(tau_z(0)) - tau_z(0) * (deltainput_s*Tmp).sum() );
    //Tmp = Q_ICAR * etainput_s;
    //if(Options_vec(3)==1) jnll -= Type(0.5) * ( (n_s-Type(1.0) ) * log(tau_z(1)) - tau_z(1) * (etainput_s*Tmp).sum() );
  }
  
  // Predicted densities
  vector<Type> Z_s(n_s);
  vector<Type> E_count(n_i);
  vector<Type> Unsampled_s(n_s);
  vector<Type> debug_j(n_j);
  vector<Type> linpredZ_s = X_sj * beta_j;
  for(int s=0; s<n_s; s++){
    Z_s(s) = exp( linpredZ_s(s) + delta_s(s) );
  }

  // Probability of counts
  for(int i=0; i<n_i; i++){
    E_count(i)=Z_s(s_i(i))*Thin_i(i);
    if( !isNA(c_i(i)) ) jnll_comp(0) -= dpois( c_i(i), E_count(i), true );
  }
  //std::cout<<Thin_i<<"\n\n";

  // Posterior predictions of abundance in unsampled areas
  for(int s=0; s<n_s; s++){
    Unsampled_s(s)=Z_s(s)*ThinU_s(s);
  }

  //Type total_abundance = Z_s.sum();
  Type total_abundance = Unsampled_s.sum() + c_sum;
  
  
  // Predicted probability of sampling
  vector<Type> R_s(n_s);
  vector<Type> linpredR_s = X_sk * beta_k;
  vector<Type> b_s = X_sb * beta_b;
  for(int s=0; s<n_s; s++){
    R_s(s) = 1 / (1 + exp(-linpredR_s(s) - b_s(s)*delta_s(s) - eta_s(s) ));
  }

  
  // Probability of sampling locations
  for(int s=0; s<n_s; s++){
    if( !isNA(y_s(s)) ) jnll_comp(1) -= dbern( y_s(s), R_s(s), true );
  }
  if(Options_vec(5)==1){
    for(int ib=0;ib<n_b;ib++){
      jnll_comp(1)+=0.5*beta_b(ib)*beta_b(ib);  // N(0,1) prior
    }
  }

  // Total objective
  Type jnll = jnll_comp.sum();

  // Reporting
  REPORT( Z_s );
  REPORT( R_s );
  REPORT( total_abundance );
  REPORT( Range_z );
  REPORT( MargSD_z );
  REPORT( beta_j );
  REPORT( beta_k );
  REPORT( beta_b );
  REPORT( delta_s );
  REPORT( eta_s );
  REPORT( linpredZ_s );
  REPORT( linpredR_s );
  REPORT( Unsampled_s );
  REPORT( jnll_comp );
  REPORT( jnll );
  
  // Bias correction output
  ADREPORT( beta_j );
  ADREPORT( total_abundance);
  if(Options_vec(4)==1){
    ADREPORT(beta_k);
    ADREPORT( beta_b ); 
    ADREPORT(Z_s);
  }
  return jnll;
}

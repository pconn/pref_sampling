#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec armaNorm(int n){
  NumericVector x = rnorm(n,0,1);
  arma::vec out(x.begin(), x.size(), false);
  return out;
}

// [[Rcpp::export]]
arma::vec GCN(const arma::mat& V, const arma::vec& v){
  arma::mat out = solve(V,v) + solve(chol(V), armaNorm(v.n_elem));
  return out;
}

// [[Rcpp::export]]
arma::vec GCN2(const arma::mat& V, const arma::vec& v){
  arma::mat out = solve(V,v) + solve(chol(V), arma::randn<arma::vec>(v.n_elem));
  return out;
}


// [[Rcpp::export]]
List PCspatregMCMC(const arma::vec& y, const arma::mat& X, const arma::mat& K, const arma::uvec& sampleIndex,
const arma::vec& betaMean, const arma::mat& betaSig, const double& phiScale, const int& burn, const int& iter) {
  int na = K.n_cols;
  int nsite = X.n_rows;
  int nb = X.n_cols;
  int ns = y.n_elem;
  arma::uvec idx = sampleIndex-1;
  arma::mat Xsmp = X.rows(idx);
  arma::mat Ksmp = K.rows(idx);
  
  arma::mat betaStor(iter, nb, fill::zeros);
  arma::mat alphaStor(iter, na, fill::zeros);
  arma::vec phiStor(iter, fill::zeros);
  arma::vec sigStor(iter, fill::zeros);
  arma::mat predStor(iter, nsite);
  
  // Beta items
  arma::mat betaPrec = betaSig.i();
  arma::mat V_beta_inv(nb, nb);
  arma::vec v_beta(nb);
  // Alpha items
  arma::mat I(na,na, fill::eye);
  arma::mat V_alpha_inv(na, na);
  arma::vec v_alpha(na);
  // phi items
  arma::vec Kbar(ns);
  double V_phi_inv;
  double v_phi;
  
  // Initial values
  arma::vec beta(nb);
  beta = solve(Xsmp.t() * Xsmp, Xsmp.t()*y);
  double tau = (ns-1)/as_scalar((y - Xsmp*beta).t()*(y - Xsmp*beta));
  double phi = 0;
  arma::vec alpha(na, fill::zeros);
  
  for(int i=0; i<iter+burn; i++){
    
    // update beta
    V_beta_inv = tau*Xsmp.t()*Xsmp + betaPrec;
    v_beta = tau * Xsmp.t() * (y - phi*Ksmp*alpha) + betaPrec*betaMean;
    beta = GCN(V_beta_inv, v_beta);
    if(i>=burn){betaStor.row(i-burn) = beta.t();}
    
    // update alpha
    V_alpha_inv = phi*phi*tau*Ksmp.t()*Ksmp + I; 
    v_alpha = phi*tau*Ksmp.t()*(y-Xsmp*beta);
    alpha = GCN(V_alpha_inv, v_alpha);
    if(i>=burn){alphaStor.row(i-burn) = phi*alpha.t();}
    
    // update phi
    Kbar = Ksmp*alpha;
    V_phi_inv = tau*as_scalar(Kbar.t()*Kbar) + 1/(phiScale*phiScale);
    v_phi = tau*dot(Kbar, y-Xsmp*beta);
    phi = as<double>(rnorm(1, v_phi/V_phi_inv, 1/sqrt(V_phi_inv)));
    if(i>=burn) phiStor(i-burn) = phi; 
    
    // update tau
    tau = as<double>(rgamma(1, ns/2, as_scalar((y-Xsmp*beta-phi*Ksmp*alpha).t()*(y-Xsmp*beta-phi*Ksmp*alpha))/2));
    if(i>=burn) sigStor(i-burn) = 1/sqrt(tau);
    
    // make prediction
    if(i>=burn){
      predStor.row(i-burn) = (X*beta + phi*K*alpha).t();
    }
    
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("beta") = betaStor, 
    Rcpp::Named("alpha") = alphaStor, 
    Rcpp::Named("phi")=phiStor,
    Rcpp::Named("sigma")=sigStor,
    Rcpp::Named("pred")=predStor
    );  
}

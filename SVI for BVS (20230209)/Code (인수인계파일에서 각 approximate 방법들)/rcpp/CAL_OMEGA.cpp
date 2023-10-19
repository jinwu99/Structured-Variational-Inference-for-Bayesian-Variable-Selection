#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix CAL_OMEGA(NumericMatrix SAMPLE_MAT,NumericVector weights, int P){
  IntegerVector cov_mem(P);
  NumericMatrix Omega(P,P);
  int total_samples = SAMPLE_MAT.nrow();
  
  for(int s=0; s<total_samples; s++){
    cov_mem = SAMPLE_MAT(s,_);
    double weight = weights(s);
    
    for(int j1=0; j1<P; j1++){
      int k1 = cov_mem(j1);
      if(k1==0) continue;
      
      Omega(j1,j1) += weight;
      
      for(int j2=(j1+1); j2<P; j2++){
        int k2 = cov_mem(j2);
        if(k2==0) continue;
        Omega(j1,j2) += weight;
        Omega(j2,j1) += weight;
      }
    }
  }
  
  return(Omega);
}

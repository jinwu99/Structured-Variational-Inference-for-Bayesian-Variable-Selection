#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix MEAN_FIELD_ALGORITHM(int Q, NumericVector psi, NumericMatrix PSI, NumericMatrix gam_matrix){
  for(int p1=0; p1<Q; p1++){
    // for star(0)
    NumericVector star(2);
    star(0) = 0;
    
    // for star(1)
    double term = psi(p1);
    
    for(int p2=0; p2<Q; p2++){
      if(p1==p2) continue;
      
      if(p2<Q){
        double gam_p2 = gam_matrix(p2,1);
        double tmp_p1p2 = 2*PSI(p1,p2);
        term += gam_p2*tmp_p1p2;
      }else continue;
    } // end p2
    star(1) = term;
    star = exp(star - max(star));
    star = star/sum(star);
    
    term = star(0); gam_matrix(p1,0) = term;
    term = star(1); gam_matrix(p1,1) = term;
  } // end p1
  return(gam_matrix);
}
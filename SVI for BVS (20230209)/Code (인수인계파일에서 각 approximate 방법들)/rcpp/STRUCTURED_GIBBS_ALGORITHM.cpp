#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace arma;

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

NumericVector EVALUATE_PSI_DIFF(NumericMatrix PSI_star, IntegerVector cov_mem, 
                                int P, int j1){
  ////////// for evaluating portion of j1 of PSI_star /////////
  NumericVector star(2);
  star(0) = 0;
  
  double term = PSI_star(j1,j1);
  for(int p2=0; p2<P; p2++){
    int k2 = cov_mem[p2];
    if(p2==j1) continue;
    if(k2==1){
      term += 2*PSI_star(j1,p2);
    }
  }
  star(1) = term;
  return(star);
}

// [[Rcpp::export]]
NumericMatrix STRUCTURED_GIBBS_ALGORITHM(NumericMatrix PSI_star,
                                         NumericVector rand_cov, IntegerVector cov_mem,
                                         int S, int P){
  NumericMatrix storage(S,P);
  NumericVector weights(S);
  
  // Then for recording trajectories :
  for(int s=0; s<S; s++){ // for each gibbs iter
    
    for(int j=0; j<P; j++){ // for sampling from each j1th cov
      int j1 = rand_cov(j);
      
      NumericVector star = EVALUATE_PSI_DIFF(PSI_star, cov_mem, P, j1);
      star = exp(star - max(star));
      star = star/sum(star);
      
      irowvec c(2);
      rmultinom(1,star.begin(),2,c.begin());
      int max_idx = c.index_max();
      cov_mem[j1] = max_idx;
    }
    
    storage(s,_) = cov_mem;
    std::random_shuffle(rand_cov.begin(), rand_cov.end(), randWrapper);
  }
  
  return(storage);
}

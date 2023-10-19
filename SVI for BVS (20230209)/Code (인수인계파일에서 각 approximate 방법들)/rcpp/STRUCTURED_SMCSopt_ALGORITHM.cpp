#include <Rcpp.h>
#include <random>
#include <RcppEigen.h>
#include <RcppThread.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppEigen)]]

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

double EVALUATE_PSI(MatrixXd& PSI_star, VectorXd& cov_mem, int P){
  double term = 0;
  ////////// for evaluating PSI_star /////////
  for(int p1=0; p1<P; p1++){
    int k1 = cov_mem(p1);
    if(k1==0) continue;
    
    double tmp = PSI_star(p1,p1);
    term += tmp;
    
    for(int p2=0; p2<P; p2++){
      int k2 = cov_mem(p2);
      if( (p2<=p1) || (k2==0)) continue;
      
      double tmp = PSI_star(p1,p2);
      term += 2*tmp;
    }
  } // end p1
  return(term);
}

double EVALUATE_PSI_DIFF(MatrixXd& PSI_star, VectorXd& cov_mem, 
                         int P, int j1){
  ////////// for evaluating portion of j1 of PSI_star /////////
  double term = PSI_star(j1,j1);
  for(int p2=0; p2<P; p2++){
    int k2 = cov_mem(p2);
    if(p2==j1) continue;
    if(k2==1){
      term += 2*PSI_star(j1,p2);
    }
  }
  return(term);
}

// [[Rcpp::export]]
List STRUCTURED_SMCSopt_ALGORITHM(std::vector<int>& rand_cov, MatrixXd& sample_mat_init,
                                     NumericVector NWeights,
                                     double tt, int particles, int sweeps, int P, double eta,
                                     MatrixXd& PSI_star_T, MatrixXd& PSI_star_O){
  // RcppEigen :
  // rand_cov, sample_mat_init, sample_mat_rs, PSI_star, PSI_star_bf, cov_mem
  
  std::random_device rd;
  std::mt19937 gen(rd());
  
  int nrow_weights = tt; // maximum iter
  
  NumericVector lweights_cml = log(NWeights);
  NumericVector log_weight_diff(particles);
  NumericVector log_weight_T(particles);
  NumericVector log_weight_t(particles);
  NumericMatrix Prop_weights(nrow_weights,particles);
  
  // for bisection
  double alpha_t = 1/tt; // assume only the first iterate is fixed
  double rho = alpha_t;
  NumericVector log_weight(particles);
  NumericVector weight(particles);
  int max_iter = 100;
  // int end = 0;
  
  MatrixXd SAMPLE_MAT(particles*nrow_weights,P);
  MatrixXd sample_mat_tmp(particles,P);
  VectorXd cov_mem(P);
  
  MatrixXd PSI_diff = PSI_star_T-PSI_star_O;
  MatrixXd PSI_star_t(P,P);
  
  int T = 0;
  int R = 0;
  for(int t=0; t<nrow_weights; t++){ // for each temperature
    // Current Weight
    for(int s=0; s<particles; s++){
      // for log(f(t)/f(t+1))[si(t+1)]
      cov_mem = sample_mat_init.row(s);
      SAMPLE_MAT.row(t*particles+s) = cov_mem;
      log_weight_diff(s) = EVALUATE_PSI(PSI_diff, cov_mem, P);
      // log_weight_T(s) = EVALUATE_PSI(PSI_star_T, cov_mem, P);
    }
    
    // Bisection
    int iter = 0;
    alpha_t = 1/tt;
    double l = 0;
    double u = 1+alpha_t;
    while(iter<max_iter){
      // obtain eta_beta_t
      log_weight = alpha_t*log_weight_diff;
      weight = exp(log_weight - max(log_weight));
      weight = weight/sum(weight);
      double ESSt = 1/sum(pow(weight,2));
      double eta_alpha_t = ESSt/particles;
      
      if(eta_alpha_t<eta){
        u = alpha_t;
        alpha_t = (alpha_t+l)/2;
      }else if(eta_alpha_t>eta){
        l = alpha_t;
        alpha_t = (alpha_t+u)/2;
      }
      
      if( (abs(u-l)<1e-5) || (l>(1-rho) || eta_alpha_t==eta ) ){
        break;
      }
      
      iter++;
    }
    alpha_t = std::min(alpha_t,(1-rho));
    rho += alpha_t;
    if(rho>=1){
      rho = 1;
      alpha_t = 0;
    }

    PSI_star_t = rho*PSI_diff + PSI_star_O;
    
    for(int s=0; s<particles; s++){
      // cov_mem = sample_mat_init.row(s);
      lweights_cml(s) += alpha_t*log_weight_diff(s);
      // log_weight_t(s) = EVALUATE_PSI(PSI_star_t, cov_mem, P);
    }
    
    // cal Weights for resampling
    NumericVector Weights = exp(lweights_cml - max(lweights_cml));
    NWeights = Weights/sum(Weights);
    lweights_cml = log(NWeights);

    for(int s=0; s<particles; s++){
      // cov_mem = sample_mat_tmp.row(s);
      cov_mem = sample_mat_init.row(s);
      // SAMPLE_MAT.row(t*particles+s) = cov_mem;
      log_weight_t(s) = EVALUATE_PSI(PSI_star_t, cov_mem, P);
      log_weight_T(s) = EVALUATE_PSI(PSI_star_T, cov_mem, P);
      // Prop_weights(t,s) = log_weight_T(s)-log_weight_t(s);
      Prop_weights(t,s) = log_weight_T(s)-log_weight_t(s)+lweights_cml(s);
    }
    
    // resample only number of num_resample particles
    double ESS = 1/sum(pow(NWeights,2));
    IntegerVector sampled_idx = sample(particles,particles,true,NWeights)-1;
    for(int s=0; s<particles; s++){
      sample_mat_tmp.row(s) = sample_mat_init.row(sampled_idx(s));
    }
    if(ESS<(particles/2)){
      R++;
      sample_mat_init = sample_mat_tmp;
    }
    
    // for(int s=0; s<particles; s++){
    //   cov_mem = sample_mat_tmp.row(s);
    //   // cov_mem = sample_mat_init.row(s);
    //   SAMPLE_MAT.row(t*particles+s) = cov_mem;
    //   log_weight_t(s) = EVALUATE_PSI(PSI_star_t, cov_mem, P);
    //   log_weight_T(s) = EVALUATE_PSI(PSI_star_T, cov_mem, P);
    //   Prop_weights(t,s) = log_weight_T(s)-log_weight_t(s);
    //   // Prop_weights(t,s) = log_weight_T(s)-log_weight_t(s)+lweights_cml(s);
    // }
    
    // Gibbs transition
    // using parallelFor
    RcppThread::parallelFor(0,particles,[&PSI_star_t, &sample_mat_init,
                            &rand_cov, &sweeps, &P, &gen] (int s) {
                              VectorXd cov_mem_s = sample_mat_init.row(s);
                              
                              std::vector<int> rand_cov_s = rand_cov;
                              std::random_shuffle(rand_cov_s.begin(), rand_cov_s.end(), randWrapper);
                              
                              for(int ss=0; ss<sweeps; ss++){
                                // sweeping
                                for(int j=0; j<P; j++){
                                  int j1 = rand_cov_s[j];
                                  
                                  double star = EVALUATE_PSI_DIFF(PSI_star_t, cov_mem_s, P, j1);
                                  if(star>=0){
                                    star = 1/(exp(-star)+1);
                                  }else{
                                    star = exp(star)/(exp(star)+1);
                                  }
                                  
                                  std::binomial_distribution<int> d(1,star);
                                  int max_idx = d(gen);
                                  cov_mem_s(j1) = max_idx;
                                } // end comb j
                              }
                              // store sweeped sample
                              sample_mat_init.row(s) = cov_mem_s;
                            });
    
    T += 1;
    // if(end==1) break;
    // T += 1;
    std::random_shuffle(rand_cov.begin(), rand_cov.end(), randWrapper);
    
    if(rho>=1){
      break;
    }
  } // end for T
  
  List final = List::create(sample_mat_init, NWeights, SAMPLE_MAT, Prop_weights, T, R);
  // should return sample_mat and all summed weight
  return(final);
}
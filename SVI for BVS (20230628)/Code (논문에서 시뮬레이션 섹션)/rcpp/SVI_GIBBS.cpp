#include <Rcpp.h>
#include <random>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

void nzero_zero_idx(const VectorXd& cov_mem, VectorXd& nzero_idx, VectorXd& zero_idx, const int& P){
  int nz = 0;
  int z = 0;
  for(int j=0; j<P; j++){
    if(cov_mem(j)==1){
      nzero_idx(nz) = j;
      nz++;
    }else{
      zero_idx(z) = j;
      z++;
    }
  }
}

MatrixXd submatrix(const MatrixXd& B, 
                   const VectorXd& cov_mem, const VectorXd& nzero_idx, const VectorXd& zero_idx,
                   const int& nzero_num, const int& zero_num, const int& P){
  MatrixXd B_r_tmp(P,nzero_num);
  MatrixXd B_r(nzero_num,nzero_num);
  // extract col
  for(int j=0; j<nzero_num; j++){
    B_r_tmp.col(j) = B.col(nzero_idx(j));
  }
  // extract row
  for(int j=0; j<nzero_num; j++){
    B_r.row(j) = B_r_tmp.row(nzero_idx(j));
  }
  return(B_r);
}

VectorXd subvector(const VectorXd& C, const VectorXd& nzero_idx, const int& nzero_num){
  VectorXd C_r(nzero_num);
  for(int j=0; j<nzero_num; j++){
    C_r(j) = C(nzero_idx(j));
  }
  return(C_r);
}

MatrixXd retrieve_mat(const MatrixXd& B_r, const VectorXd& inv_V,
                      const VectorXd& nzero_idx, const VectorXd& zero_idx,
                      const int& nzero_num, const int& zero_num, const int& P){
  MatrixXd B_t(P,P);
  B_t.setZero();
  for(int j1=0; j1<nzero_num; j1++){
    int org_idx1 = nzero_idx(j1);
    B_t(org_idx1,org_idx1) = B_r(j1,j1);
    for(int j2=0; j2<nzero_num; j2++){
      if(j2<=j1) continue;
      int org_idx2 = nzero_idx(j2);
      B_t(org_idx1,org_idx2) = B_r(j1,j2);
      B_t(org_idx2,org_idx1) = B_r(j2,j1);
    }
  }
  for(int j=0; j<zero_num; j++){
    B_t(zero_idx(j),zero_idx(j)) = 1/inv_V(zero_idx(j));
  }
  return(B_t);
}

MatrixXd retrieve_mat2(const VectorXd& inv_V, const int& P){
  MatrixXd V(P,P);
  V.setZero();
  for(int j=0; j<P; j++){
    V(j,j) = 1/inv_V(j);
  }
  return(V);
}

VectorXd expand_vec(const VectorXd& C_r, const VectorXd& nzero_idx, const int& nzero_num, const int& P){
  VectorXd C_t(P);
  C_t.setZero();
  for(int j=0; j<nzero_num; j++){
    int org_idx = nzero_idx(j);
    C_t(org_idx) = C_r(j);
  }
  return(C_t);
}

void log_prob_comp(std::vector<int>& rand_cov, VectorXd& cov_mem, int& nzero_num, int& zero_num,
                   VectorXd& C, MatrixXd& B, VectorXd& inv_V,
                   int& j, int& P,
                   double& log_det_B, double& CTBiC_r,
                   VectorXd& Theta_mu, MatrixXd& B_inv){
  
  nzero_num = cov_mem.sum();
  zero_num = P-nzero_num;
  if(nzero_num>0){
    // figuring indexes
    VectorXd nzero_idx(nzero_num);
    VectorXd zero_idx(zero_num);
    nzero_zero_idx(cov_mem,nzero_idx,zero_idx,P);
    
    // extract only-need components
    MatrixXd B_r = submatrix(B,cov_mem,nzero_idx,zero_idx,nzero_num,zero_num,P);
    VectorXd C_r = subvector(C,nzero_idx,nzero_num);
    LLT<MatrixXd> llt(B_r);
    
    // log_determinant : log(det(Br)) + sum(log(inv_Vo))
    double log_det_B_r_d = 2*MatrixXd(llt.matrixL()).diagonal().array().log().sum();
    double log_det_inv_V = 0;
    for(int j1=0; j1<zero_num; j1++)
      log_det_inv_V += log(inv_V(zero_idx(j1)));
    log_det_B = log_det_B_r_d + log_det_inv_V;
    
    // CTBiC
    VectorXd Theta_mu_r = llt.solve(C_r);
    MatrixXd CTBiC = C_r.adjoint()*Theta_mu_r;
    CTBiC_r = CTBiC(0,0);
    if(j==(P-1)){
      Theta_mu = expand_vec(Theta_mu_r,nzero_idx,nzero_num,P);
      MatrixXd B_r_inv = B_r.inverse();
      B_inv = retrieve_mat(B_r_inv,inv_V,nzero_idx,zero_idx,nzero_num,zero_num,P);
    }
  }else{
    // all the elements of cov_mem are 0
    // only inv_V is involved
    log_det_B = inv_V.array().log().sum();
    CTBiC_r = 0;
    if(j==(P-1)){
      Theta_mu.setZero();
      B_inv = retrieve_mat2(inv_V,P);
    }
  }
}

// [[Rcpp::export]]
List SVI_GIBBS(std::vector<int>& rand_cov, MatrixXd& B, VectorXd& C, VectorXd& inv_V,
                 double E_logit, double is2_a, double yTy,
                 VectorXd& cov_mem, int S){
  // Note that B = WTW + Inv_V
  std::random_device rd;
  std::mt19937 gen(rd());
  
  int P = B.cols();
  // Storage
  MatrixXd sample_mat(S,P);
  MatrixXd Theta_mu_lt(S,P);
  List B_inv_lt(0);
  VectorXd inv_sigma2_lt(S);
  VectorXd is2_b(S);
  // for log_prob
  double log_det_B_c = 0;
  double log_det_B_n = 0;
  double CTBiC_c = 0;
  double CTBiC_n = 0;
  int nzero_num = 0;
  int zero_num = 0;
  int nzero_num_n = 0;
  int zero_num_n = 0;
  // for Theta_mu,cov
  VectorXd Theta_mu_c(P);
  VectorXd Theta_mu_n(P);
  MatrixXd B_inv_c(P,P);
  MatrixXd B_inv_n(P,P);
  
  VectorXd log_Q_p(S);
  
  for(int s=0; s<S; s++){
    for(int j=0; j<P; j++){
      int current = cov_mem(rand_cov[j]);
      
      //////////////////// Current //////////////////////
      nzero_num = cov_mem.sum();
      zero_num = P-nzero_num;
      log_prob_comp(rand_cov,cov_mem,nzero_num,zero_num,
                    C,B,inv_V,
                    j,P,
                    log_det_B_c,CTBiC_c,
                    Theta_mu_c,B_inv_c);
      
      //////////////////// Swapped //////////////////////
      VectorXd cov_mem_n = cov_mem;
      cov_mem_n(rand_cov[j]) = abs(cov_mem(rand_cov[j])-1);
      nzero_num_n = cov_mem_n.sum();
      zero_num_n = P-nzero_num_n;
      log_prob_comp(rand_cov,cov_mem_n,nzero_num_n,zero_num_n,
                    C,B,inv_V,
                    j,P,
                    log_det_B_n,CTBiC_n,
                    Theta_mu_n,B_inv_n);
      
      NumericVector star(2);
      if(current==0){
        star(0) = -0.5*log_det_B_c + E_logit*nzero_num - is2_a*log(0.5*(yTy-CTBiC_c));
        star(1) = -0.5*log_det_B_n + E_logit*nzero_num_n - is2_a*log(0.5*(yTy-CTBiC_n));
      }else{ // current==1
        star(1) = -0.5*log_det_B_c + E_logit*nzero_num - is2_a*log(0.5*(yTy-CTBiC_c));
        star(0) = -0.5*log_det_B_n + E_logit*nzero_num_n - is2_a*log(0.5*(yTy-CTBiC_n));
      }
      NumericVector star_norm = exp(star - max(star));
      star_norm = star_norm/sum(star_norm);
      std::binomial_distribution<int> d(1,star_norm(1));
      int sample = d(gen);
      cov_mem(rand_cov[j]) = sample;
      
      if(j==(P-1)){
        log_Q_p(s) = star(sample);
        sample_mat.row(s) = cov_mem;
        if(current==sample){
          Theta_mu_lt.row(s) = Theta_mu_c;
          B_inv_lt.push_back(B_inv_c);
          inv_sigma2_lt(s) = is2_a/(0.5*(yTy-CTBiC_c));
          is2_b(s) = 0.5*(yTy-CTBiC_c);
        }
        else{
          Theta_mu_lt.row(s) = Theta_mu_n;
          B_inv_lt.push_back(B_inv_n);
          inv_sigma2_lt(s) = is2_a/(0.5*(yTy-CTBiC_n));
          is2_b(s) = 0.5*(yTy-CTBiC_n);
        }
      }
    } // end j
    std::random_shuffle(rand_cov.begin(), rand_cov.end(), randWrapper);
  }
  
  List final = List::create(sample_mat,Theta_mu_lt,B_inv_lt,inv_sigma2_lt,log_Q_p,is2_b);
  return(final);
}
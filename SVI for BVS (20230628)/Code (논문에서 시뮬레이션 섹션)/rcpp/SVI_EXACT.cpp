#include <Rcpp.h>
#include <random>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]

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

void log_prob_comp(VectorXd& cov_mem, int& nzero_num, int& zero_num,
                   VectorXd& C, MatrixXd& B, VectorXd& inv_V,
                   int& P, double& log_det_B, double& CTBiC_r,
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
    
    Theta_mu = expand_vec(Theta_mu_r,nzero_idx,nzero_num,P);
    MatrixXd B_r_inv = B_r.inverse();
    B_inv = retrieve_mat(B_r_inv,inv_V,nzero_idx,zero_idx,nzero_num,zero_num,P);
  }else{
    // all the elements of cov_mem are 0
    // only inv_V is involved
    log_det_B = inv_V.array().log().sum();
    CTBiC_r = 0;
    Theta_mu.setZero();
    B_inv = retrieve_mat2(inv_V,P);
  }
}

// [[Rcpp::export]]
List SVI_EXACT(MatrixXd& Comb, MatrixXd& B, VectorXd& C, VectorXd& inv_V,
               double E_logit, double is2_a, double yTy){
  // Note that B = WTW + Inv_V
  int G = Comb.rows();
  int P = B.cols();
  // Storage
  NumericVector Comb_sum(G);
  MatrixXd Theta_mu_lt(G,P);
  List B_inv_lt(0);
  VectorXd inv_sigma2_lt(G);
  // for log_prob
  double log_det_B = 0;
  double CTBiC = 0;
  int nzero_num = 0;
  int zero_num = 0;;
  // for Theta_mu,cov
  VectorXd Theta_mu(P);
  MatrixXd B_inv(P,P);
  
  for(int g=0; g<G; g++){
    VectorXd cov_mem = Comb.row(g);
    
    nzero_num = cov_mem.sum();
    zero_num = P-nzero_num;
    log_prob_comp(cov_mem,nzero_num,zero_num,
                  C,B,inv_V,
                  P,log_det_B,CTBiC,
                  Theta_mu,B_inv);
    
    Comb_sum(g) = -0.5*log_det_B + E_logit*nzero_num - is2_a*log(0.5*(yTy-CTBiC));
    
    Theta_mu_lt.row(g) = Theta_mu;
    B_inv_lt.push_back(B_inv);
    inv_sigma2_lt(g) = is2_a/(0.5*(yTy-CTBiC));
  }
  
  List final = List::create(Comb_sum,Theta_mu_lt,B_inv_lt,inv_sigma2_lt);
  return(final);
}
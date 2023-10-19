rm(list=ls())
library(copula)
library(Matrix)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(mnormt) # for pd.solve
sourceCpp('MEAN_FIELD_ALGORITHM.cpp')
sourceCpp('STRUCTURED_GIBBS_ALGORITHM.cpp')
sourceCpp('CAL_OMEGA.cpp')
sourceCpp('STRUCTURED_SMCSopt_ALGORITHM.cpp')

################################################

# global functions
normalize <- function(x) { return ((x - min(x)) / (max(x) - min(x))) }
standardize <- function(x) { return ((x - mean(x)) / sd(x))}
tol <- 1e-200
normalizing <- function(mat) mat/rowSums(mat)
mat_trace <- function(A) sum(diag(A))
ar1_cor <- function(n, rho) rho^abs(matrix(1:n-1,nrow=n,ncol=n,byrow=TRUE)-(1:n-1))
logit <- function(m) log((m+tol)/(1+m+tol))
log_sum_exp <- function(log_vec) log(sum(exp(log_vec-max(log_vec)))+tol) + max(log_vec)

################################################
################################################

# Data generating functions
# sim_data_corr_all : generating data with all equally correlated covariates
# sim_data_corr_all : generating data with AR(1) covariates
# parameters :
# Q = number of covaraites
# phi = correlation
# N = number of observations
# s = variance of error
# r = random seed
# p = number of nonzero coefficients
sim_data_corr_all <- function(Q,phi,N,s,r,p){
  set.seed(r)
  
  ########### shuffering linear
  rand_idx_x <- sample(c(1:Q),p)
  true.coeff <- runif(p,1,2)
  
  minus_idx <- unique(sample.int(p,p,replace=T))
  true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[rand_idx_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[rand_idx_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=Q,ncol=Q)
  corr.x <- phi*J - phi*diag(Q) + diag(Q)
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff))
}

sim_data_corr_ar1 <- function(Q,phi,N,s,r,p){
  set.seed(r)
  
  ########### shuffering linear
  rand_idx_x <- sample(c(1:Q),p)
  true.coeff <- runif(p,1,2)
  
  minus_idx <- unique(sample.int(p,p,replace=T))
  true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[rand_idx_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[rand_idx_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=Q,ncol=Q)
  corr.x <- ar1_cor(Q,phi)
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff))
}

################################################
################################################

# Algorithms
# MFVI : Mean Field variational Inferernce with Deterministic Annealing (DA)
# SVIG : Structured variational Inference with Gibbs Sampling (GS) and DA
# SVIS : Structured Variational Inference with Sequential Monte Carlo Sampler (SMCS) with GS and DA
# SVIM : Structured Variational Inference without MCMC (developing)

MFVI <- function(y,X,mem_true){
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  WTW <- t(W)%*%W
  
  #========== Iinitial VI parameter values ==========#
  #===== error =====#
  inv_sigma2 <- 1/sd(y)
  
  #===== coefficients =====#
  Theta.mu <- matrix(0,nrow=Q,ncol=1)
  
  #===== rho =====#
  A0 <- 1
  B0 <- Q
  
  #===== spike and slab prior =====#
  gam_matrix <- matrix(1/2,Q,2)
  gam <- matrix(gam_matrix[,2],ncol=1)
  Gam <- diag(c(gam),ncol=Q)
  Gam_2 <- diag(c(gam^2),ncol=Q)
  Omega <- gam%*%t(gam) - Gam_2 + Gam
  WTWO <- WTW*Omega
  
  #===== Dispersion parameters =====#
  lamb2_2 <- 1
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #===== Intercept =====#
  a.mu <- mean(y)
  a.sigma2 <- 1
  j1 <- matrix(1,ncol=1,nrow=N)
  yTW <- t(y-a.mu)%*%W
  
  #===== Deterministic Annealing(if needed)=====#
  # xi <- 1
  xi <- 0.001
  
  #========== Optimization ==========#
  num_iters <- 100
  before <- rep(1/2,Q)
  after <- rep(1/2,Q)
  diff <- c()
  
  for(idx in 1:num_iters){
    #===== Intercept update =====#
    a.sigma2 <- (inv_sigma2*N)^(-1)
    a.mu <- ( a.sigma2*(inv_sigma2*t(j1)%*%(y-W%*%(gam*Theta.mu)) ) )[1]
    a.mu2 <- a.sigma2 + a.mu^2
    yTW <- t(y-a.mu)%*%W
    
    Theta.cov <- pd.solve(inv_sigma2*(WTWO+Inv_V))
    Theta.mu <- inv_sigma2*gam*t(yTW); Theta.mu <- Theta.cov%*%Theta.mu
    Theta_ThetaT <- Theta.cov + Theta.mu%*%t(Theta.mu)
    theta.p2 <- diag(Theta_ThetaT)
    
    #===== Dispersion paramter =====#
    b_tau <- inv_sigma2*theta.p2
    inv_tau2 <- (sqrt(a_tau)*besselK(sqrt(a_tau*b_tau),p_tau+1))/
                (sqrt(b_tau)*besselK(sqrt(a_tau*b_tau),p_tau)) - ((2*p_tau)/b_tau)
    Inv_V <- diag(inv_tau2,Q)
    
    #===== logVk update =====# ***********
    A <- A0+sum(gam_matrix[,2]); B <- B0+(Q-sum(gam_matrix[,2]))
    E_logit <- digamma(A)-digamma(B)
    
    #===== Inv_sigma2 update =====#
    TMP <- sum((WTWO+Inv_V)*Theta_ThetaT)
    is2.a <- (N+Q)/2
    is2.b <- c(t(y)%*%y-2*a.mu*t(y)%*%j1+N*a.mu2 -
                 2*yTW%*%(gam*Theta.mu) + TMP)/2
    inv_sigma2 <- is2.a/is2.b
    
    #===== Gam update =====#
    PSI <- xi*(-0.5*inv_sigma2*(WTW*Theta_ThetaT)); psi <- diag(PSI); diag(PSI) <- 0
    psi <- xi*(E_logit + inv_sigma2*(c(yTW)*c(Theta.mu))) + psi
    gam_matrix <- MEAN_FIELD_ALGORITHM(Q,psi,PSI,gam_matrix)
    gam <- matrix(gam_matrix[,2],ncol=1)
    Gam <- diag(c(gam),ncol=Q)
    Gam_2 <- diag(c(gam^2),ncol=Q)
    Omega <- gam%*%t(gam) - Gam_2 + Gam
    WTWO <- WTW*Omega
    
    #==== xi update ====#
    xi <- xi+0.1
    if(xi>=1) xi<-1
    
    mem_est <- max.col(gam_matrix,'first')
    prob<- apply(gam_matrix,1,max)
    diff[idx] <- sum(abs(mem_true-mem_est))
    
    after <- c(-gam*log(gam+1e-10)-(1-gam)*log(1-gam+1e-10))
    delta_H <- max(abs(after-before)) 
    if(delta_H<=1e-4) break
    before <- after
  }
  print(diff)
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  
  resid <- y-(a.mu+W%*%(gam*Theta.mu))
  mse <- sqrt(norm(resid,type='2')/N)
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1))
}

SVIG <- function(y,X,mem_true,S){
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  WTW <- t(W)%*%W
  
  #========== Iinitial VI parameter values ==========#
  #===== error =====#
  inv_sigma2 <- 1/sd(y)
  
  #===== coefficients =====#
  Theta.mu <- matrix(0,nrow=Q,ncol=1)
  
  #===== rho =====#
  A0 <- 1
  B0 <- Q
  
  #===== spike and slab prior =====#
  gam_matrix <- matrix(1/2,Q,2)
  gam <- matrix(gam_matrix[,2],ncol=1)
  Gam <- diag(c(gam),ncol=Q)
  Gam_2 <- diag(c(gam^2),ncol=Q)
  Omega <- gam%*%t(gam) - Gam_2 + Gam
  WTWO <- WTW*Omega
  
  #===== Dispersion parameters =====#
  lamb2_2 <- 1
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #===== Intercept =====#
  a.mu <- mean(y)
  a.sigma2 <- 1
  j1 <- matrix(1,ncol=1,nrow=N)
  yTW <- t(y-a.mu)%*%W
  
  #===== Deterministic Annealing(if needed)=====#
  # xi <- 1
  xi <- 0.001
  
  #===== For gibbs sampling =====#
  cov_mem <- as.numeric(sample.int(2,size=Q,replace=T,prob=c(1/2,1/2))==2)
  
  #========== Optimization ==========#
  num_iters <- 100
  before <- rep(1/2,Q)
  after <- rep(1/2,Q)
  diff <- c()
  
  for(idx in 1:num_iters){
    #===== Intercept update =====#
    a.sigma2 <- (inv_sigma2*N)^(-1)
    a.mu <- ( a.sigma2*(inv_sigma2*t(j1)%*%(y-W%*%(gam*Theta.mu)) ) )[1]
    a.mu2 <- a.sigma2 + a.mu^2
    yTW <- t(y-a.mu)%*%W
    
    Theta.cov <- pd.solve(inv_sigma2*(WTWO+Inv_V))
    Theta.mu <- inv_sigma2*gam*t(yTW); Theta.mu <- Theta.cov%*%Theta.mu
    Theta_ThetaT <- Theta.cov + Theta.mu%*%t(Theta.mu)
    theta.p2 <- diag(Theta_ThetaT)
    
    #===== Dispersion paramter =====#
    b_tau <- inv_sigma2*theta.p2
    inv_tau2 <- (sqrt(a_tau)*besselK(sqrt(a_tau*b_tau),p_tau+1))/
      (sqrt(b_tau)*besselK(sqrt(a_tau*b_tau),p_tau)) - ((2*p_tau)/b_tau)
    Inv_V <- diag(inv_tau2,Q)
    
    #===== logVk update =====# ***********
    A <- A0+sum(gam_matrix[,2]); B <- B0+(Q-sum(gam_matrix[,2]))
    E_logit <- digamma(A)-digamma(B)
    
    #===== Inv_sigma2 update =====#
    TMP <- sum((WTWO+Inv_V)*Theta_ThetaT)
    is2.a <- (N+Q)/2
    is2.b <- c(t(y)%*%y-2*a.mu*t(y)%*%j1+N*a.mu2 -
                 2*yTW%*%(gam*Theta.mu) + TMP)/2
    inv_sigma2 <- is2.a/is2.b
    
    #===== Gam update =====#
    # Note that PSI_star = psi + PSI
    PSI_star <- xi*(diag(E_logit + inv_sigma2*(c(yTW)*c(Theta.mu)))
                    -0.5*inv_sigma2*(WTW*Theta_ThetaT))
    # random scan gibbs sampling
    rand_cov_x <- sample(seq.int(from=1,to=Q),replace=F)-1
    sample_mat <- STRUCTURED_GIBBS_ALGORITHM(PSI_star, rand_cov_x, cov_mem, S*1.1, Q)
    sample_mat <- sample_mat[-c(1:(S*0.1)),]
    
    Omega <- CAL_OMEGA(sample_mat, rep(1/dim(sample_mat)[1],dim(sample_mat)[1]), Q)
    gam <- diag(Omega)
    gam_matrix[,2] <- gam; gam_matrix[,1] <- 1-gam
    cov_mem <- tail(sample_mat,1)
    WTWO <- WTW*Omega
    
    #==== xi update ====#
    xi <- xi*1.2
    if(xi>=1) xi<-1
    
    mem_est <- max.col(gam_matrix,'first')
    prob<- apply(gam_matrix,1,max)
    diff[idx] <- sum(abs(mem_true-mem_est))
    
    after <- c(-gam*log(gam+1e-10)-(1-gam)*log(1-gam+1e-10))
    delta_H <- max(abs(after-before)) 
    if(delta_H<=1e-4) break
    before <- after
  }
  print(diff)
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  
  resid <- y-(a.mu+W%*%(gam*Theta.mu))
  mse <- sqrt(norm(resid,type='2')/N)
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1))
}

SVIS <- function(y,X,mem_true,S,tt,sweeps,eta){
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  WTW <- t(W)%*%W
  
  #========== Iinitial VI parameter values ==========#
  #===== error =====#
  inv_sigma2 <- 1/sd(y)
  
  #===== coefficients =====#
  Theta.mu <- matrix(0,nrow=Q,ncol=1)
  
  #===== rho =====#
  A0 <- 1
  B0 <- Q
  
  #===== spike and slab prior =====#
  gam_matrix <- matrix(1/2,Q,2)
  gam <- matrix(gam_matrix[,2],ncol=1)
  Gam <- diag(c(gam),ncol=Q)
  Gam_2 <- diag(c(gam^2),ncol=Q)
  Omega <- gam%*%t(gam) - Gam_2 + Gam
  WTWO <- WTW*Omega
  
  #===== Dispersion parameters =====#
  lamb2_2 <- 1
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #===== Intercept =====#
  a.mu <- mean(y)
  a.sigma2 <- 1
  j1 <- matrix(1,ncol=1,nrow=N)
  yTW <- t(y-a.mu)%*%W
  
  #===== Deterministic Annealing(if needed)=====#
  # xi <- 1
  xi <- 0.001
  
  #===== For SMCS =====#
  PSI_star_bf <- diag(rep(log(1/2),Q),nrow=Q)
  sample_mat_init <- t(sapply(1:S, function(s) sample.int(2,size=Q,replace=T,prob=c(1/2,1/2))))-1
  SAMPLE_MAT <- matrix(NA,nrow=0,ncol=Q)
  log_Prop_weights <- matrix(NA,nrow=0,ncol=S)
  # cooling_schedule <- c(0,ppoints((tt-1),a=0),1)
  beta <- 1/tt
  NWeights <- rep(1,S)
  thin <- 1; B <- 0
  h <- c()
  TTT <- c()
  RRR <- c()
  
  #========== Optimization ==========#
  num_iters <- 100
  before <- rep(1/2,Q)
  after <- rep(1/2,Q)
  diff <- c()
  
  for(idx in 1:num_iters){
    #===== Intercept update =====#
    a.sigma2 <- (inv_sigma2*N)^(-1)
    a.mu <- ( a.sigma2*(inv_sigma2*t(j1)%*%(y-W%*%(gam*Theta.mu)) ) )[1]
    a.mu2 <- a.sigma2 + a.mu^2
    yTW <- t(y-a.mu)%*%W
    
    Theta.cov <- pd.solve(inv_sigma2*(WTWO+Inv_V))
    Theta.mu <- inv_sigma2*gam*t(yTW); Theta.mu <- Theta.cov%*%Theta.mu
    Theta_ThetaT <- Theta.cov + Theta.mu%*%t(Theta.mu)
    theta.p2 <- diag(Theta_ThetaT)
    
    #===== Dispersion paramter =====#
    b_tau <- inv_sigma2*theta.p2
    inv_tau2 <- (sqrt(a_tau)*besselK(sqrt(a_tau*b_tau),p_tau+1))/
      (sqrt(b_tau)*besselK(sqrt(a_tau*b_tau),p_tau)) - ((2*p_tau)/b_tau)
    Inv_V <- diag(inv_tau2,Q)
    
    #===== logVk update =====# ***********
    A <- A0+sum(gam_matrix[,2]); B <- B0+(Q-sum(gam_matrix[,2]))
    E_logit <- digamma(A)-digamma(B)
    
    #===== Inv_sigma2 update =====#
    TMP <- sum((WTWO+Inv_V)*Theta_ThetaT)
    is2.a <- (N+Q)/2
    is2.b <- c(t(y)%*%y-2*a.mu*t(y)%*%j1+N*a.mu2 -
                 2*yTW%*%(gam*Theta.mu) + TMP)/2
    inv_sigma2 <- is2.a/is2.b
    
    #===== Gam update =====#
    PSI_star <- xi*(diag(E_logit + inv_sigma2*(c(yTW)*c(Theta.mu)))
                    -0.5*inv_sigma2*(WTW*Theta_ThetaT))

    # random scan gibbs sampling
    rand_cov_x <- sample(seq.int(from=1,to=Q),replace=F)-1
    LIST <- STRUCTURED_SMCSopt_ALGORITHM(rand_cov_x, sample_mat_init, NWeights,
                                         tt, S, sweeps, Q, eta, PSI_star, PSI_star_bf)
    ttt <- LIST[[5]]; TTT[idx] <- ttt
    R <- LIST[[6]]; RRR[idx] <- R
    sample_mat_init <- LIST[[1]]
    NWeights <- LIST[[2]]
    SAMPLE_MAT <- LIST[[3]]; log_wESS <- LIST[[4]]; 
    SAMPLE_MAT <- SAMPLE_MAT[1:(ttt*S),]; log_wESS <- log_wESS[1:ttt,]
    
    if(ttt==1){
      wESS <- exp(log_wESS-max(log_wESS))
      nwESS <- wESS/sum(wESS)
      h <- nwESS
    }else{
      log_WESSt <- sapply(1:ttt, function(t) log_sum_exp(log_wESS[t,]))
      log_wESS2 <- sapply(1:ttt, function(t) log_sum_exp(2*log_wESS[t,]))
      log_lt <- 2*log_WESSt - log_wESS2
      log_lambda_star <- log_lt - log_sum_exp(log_lt)
      lambda_star <- exp(log_lambda_star)
      
      wESS <- t(sapply(1:ttt, function(t) exp(log_wESS[t,]-max(log_wESS[t,]))))
      nwESS <- t(sapply(1:ttt, function(t) wESS[t,]/sum(wESS[t,])))
      h <- t(sapply(1:ttt, function(t) lambda_star[t]*nwESS[t,]))
    }
    
    # plot(c(t(h)))
    Omega <- CAL_OMEGA(SAMPLE_MAT, c(t(h)), Q)
    gam <- diag(Omega)
    gam_matrix[,2] <- gam; gam_matrix[,1] <- 1-gam
    WTWO <- WTW*Omega
    PSI_star_bf <- PSI_star
    
    #==== xi update ====#
    xi <- xi*1.2
    if(xi>=1) xi<-1
    
    mem_est <- max.col(gam_matrix,'first')
    prob<- apply(gam_matrix,1,max)
    diff[idx] <- sum(abs(mem_true-mem_est))
    
    after <- c(-gam*log(gam+1e-10)-(1-gam)*log(1-gam+1e-10))
    delta_H <- max(abs(after-before)) 
    if(delta_H<=1e-4) break
    before <- after
  }
  print(rbind(rbind(diff,TTT),RRR))
  # print(diff)
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  
  resid <- y-(a.mu+W%*%(gam*Theta.mu))
  mse <- sqrt(norm(resid,type='2')/N)
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,TTT))
}

SVIM <- function(y,X,mem_true,S,tt){
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  WTW <- t(W)%*%W
  
  #========== Iinitial VI parameter values ==========#
  #===== error =====#
  inv_sigma2 <- 1/sd(y)
  
  #===== coefficients =====#
  Theta.mu <- matrix(0,nrow=Q,ncol=1)
  
  #===== rho =====#
  A0 <- 1
  B0 <- Q
  
  #===== spike and slab prior =====#
  gam_matrix <- matrix(1/2,Q,2)
  gam <- matrix(gam_matrix[,2],ncol=1)
  Gam <- diag(c(gam),ncol=Q)
  Gam_2 <- diag(c(gam^2),ncol=Q)
  Omega <- gam%*%t(gam) - Gam_2 + Gam
  WTWO <- WTW*Omega
  
  #===== Dispersion parameters =====#
  lamb2_2 <- 1
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #===== Intercept =====#
  a.mu <- mean(y)
  a.sigma2 <- 1
  j1 <- matrix(1,ncol=1,nrow=N)
  yTW <- t(y-a.mu)%*%W
  
  #===== Deterministic Annealing(if needed)=====#
  # xi <- 1
  xi <- 0.001
  
  #===== For Monte Carlo sampling =====#
  sample_mat <- t(sapply(1:S, function(s) as.numeric(sample.int(2,size=Q,replace=T,prob=c(1/2,1/2))==2)))
  schedule <- c(rep(0,5),ppoints((tt-1),a=0),1)
  unique_num_arr <- c()
  
  #========== Optimization ==========#
  num_iters <- 100
  before <- rep(1/2,Q)
  after <- rep(1/2,Q)
  diff <- c()
  
  for(idx in 1:num_iters){
    #===== Intercept update =====#
    a.sigma2 <- (inv_sigma2*N)^(-1)
    a.mu <- ( a.sigma2*(inv_sigma2*t(j1)%*%(y-W%*%(gam*Theta.mu)) ) )[1]
    a.mu2 <- a.sigma2 + a.mu^2
    yTW <- t(y-a.mu)%*%W
    
    Theta.cov <- pd.solve(inv_sigma2*(WTWO+Inv_V))
    Theta.mu <- inv_sigma2*gam*t(yTW); Theta.mu <- Theta.cov%*%Theta.mu
    Theta_ThetaT <- Theta.cov + Theta.mu%*%t(Theta.mu)
    theta.p2 <- diag(Theta_ThetaT)
    
    #===== Dispersion paramter =====#
    b_tau <- inv_sigma2*theta.p2
    inv_tau2 <- (sqrt(a_tau)*besselK(sqrt(a_tau*b_tau),p_tau+1))/
      (sqrt(b_tau)*besselK(sqrt(a_tau*b_tau),p_tau)) - ((2*p_tau)/b_tau)
    Inv_V <- diag(inv_tau2,Q)
    
    #===== logVk update =====#
    A <- A0+sum(gam_matrix[,2]); B <- B0+(Q-sum(gam_matrix[,2]))
    E_logit <- digamma(A)-digamma(B)
    
    #===== Inv_sigma2 update =====#
    TMP <- sum((WTWO+Inv_V)*Theta_ThetaT)
    is2.a <- (N+Q)/2
    is2.b <- c(t(y)%*%y-2*a.mu*t(y)%*%j1+N*a.mu2 -
                 2*yTW%*%(gam*Theta.mu) + TMP)/2
    inv_sigma2 <- is2.a/is2.b
    
    #===== Gam update =====#
    # Note that PSI_star = psi + PSI
    PSI_star <- xi*(diag(E_logit + inv_sigma2*(c(yTW)*c(Theta.mu)))
                    -0.5*inv_sigma2*(WTW*Theta_ThetaT))
    
    if(idx==1){
      log_weights_t <- sapply(1:dim(sample_mat)[1], function(s){
        t(sample_mat[s,])%*%PSI_star%*%sample_mat[s,]
      })
      weights_t <- exp(log_weights_t-max(log_weights_t))
      nweights_t <- weights_t/sum(weights_t)
      sample_mat <- sample_mat[sample.int(dim(sample_mat)[1],dim(sample_mat)[1],
                                          replace=T,prob=nweights_t),]
    }
    
    for(t in 1:(length(schedule))){
      sample_mat <- t(sapply(1:S, function(s) abs(sample_mat[s,] - rbinom(n=Q,size=1,prob=1/Q))))
      
      log_weights_t <- sapply(1:dim(sample_mat)[1], function(s){
        t(sample_mat[s,])%*%PSI_star%*%sample_mat[s,]
      })
      weights_t <- exp(log_weights_t-max(log_weights_t))
      nweights_t <- weights_t/sum(weights_t)
      sample_mat <- sample_mat[sample.int(dim(sample_mat)[1],dim(sample_mat)[1],
                                          replace=T,prob=nweights_t),]
    }
    
    # filtering unique samples
    sample_mat_unique <- unique(sample_mat)
    unique_num <- dim(sample_mat_unique)[1]
    unique_count <- c()
    for(j in 1:unique_num){
      vec <- sample_mat_unique[j,]
      c <- 0
      for(i in 1:dim(sample_mat)[1]){
        vec_comp <- sample_mat[i,]
        if(sum(vec_comp==vec)==Q){
          c <- c+1
        }
      }
      unique_count[j] <- c
    }
    
    g_weights <- log(unique_count/unique_num)
    log_weights <- sapply(1:unique_num, function(s){
      t(sample_mat_unique[s,])%*%PSI_star%*%sample_mat_unique[s,] - g_weights[s]
    })
    weights <- exp(log_weights-max(log_weights))
    nweights <- weights/sum(weights)
    mat <- cbind(sample_mat_unique,nweights)
    Om22 <- matrix(0,ncol=Q,nrow=Q)
    for(m in 1:dim(mat)[1]){
      mem2 <- which(mat[m,1:Q]==1)
      if(length(mem2)>0){
        Om22[mem2,mem2] <- Om22[mem2,mem2] + mat[m,(Q+1)]
      }
    }
    Omega <- Om22
    gam <- diag(Om22)
    gam_matrix[,2] <- gam; gam_matrix[,1] <- 1-gam
    WTWO <- WTW*Omega
    
    #==== xi update ====#
    xi <- xi*1.2
    if(xi>=1) xi<-1
    
    mem_est <- max.col(gam_matrix,'first')
    prob<- apply(gam_matrix,1,max)
    diff[idx] <- sum(abs(mem_true-mem_est))
    
    after <- c(-gam*log(gam+1e-10)-(1-gam)*log(1-gam+1e-10))
    delta_H <- max(abs(after-before)) 
    if(delta_H<=1e-4) break
    before <- after
  }
  print(diff)
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  
  resid <- y-(a.mu+W%*%(gam*Theta.mu))
  mse <- sqrt(norm(resid,type='2')/N)
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1))
}

################################################
################################################

MISS1 <- matrix(0,nrow=0,ncol=3)
MISS2 <- matrix(0,nrow=0,ncol=3)
MISS3 <- matrix(0,nrow=0,ncol=3)
MISS4 <- matrix(0,nrow=0,ncol=3)

time1 <- matrix(0,nrow=0,ncol=1)
time2 <- matrix(0,nrow=0,ncol=1)
time3 <- matrix(0,nrow=0,ncol=1)
time4 <- matrix(0,nrow=0,ncol=1)

rho<-0.8; Q<-200; N<-100; s<-1; p<-20
n <- 5

SIM <- 1

for(r in 1:SIM){
  print('SIMULATION : '); print(r)
  List <- sim_data_corr_ar1(Q,rho,N,s,r,p)
  # List <- sim_data_corr_all(Q,rho,N,s,r,p)
  y <- List[[1]]; X <- List[[2]]
  data <- List[[3]]; mem_true <- List[[4]]; coeff <- List[[5]]
  
  print('MFVI'); pt<-proc.time(); ls1 <- MFVI(y,X,mem_true); ptFinal<-proc.time()-pt
  time1 <- rbind(time1,c(ptFinal[3])); print(c(ptFinal[3]))
  
  S <- 100
  print('SVIG'); pt<-proc.time(); ls2 <- SVIG(y,X,mem_true,S); ptFinal<-proc.time()-pt
  time2 <- rbind(time2,c(ptFinal[3])); print(c(ptFinal[3]))

  S <- 10; tt <- 200; sweeps <- 1
  eta <-1-(1e-3)
  print('SVIS'); pt<-proc.time(); ls3 <- SVIS(y,X,mem_true,S,tt,sweeps,eta); ptFinal<-proc.time()-pt
  time3 <- rbind(time3,c(ptFinal[3])); print(c(ptFinal[3]))
  
  S <- 100; tt <- 10
  print('SVIM'); pt<-proc.time(); ls4 <- SVIM(y,X,mem_true,S,tt); ptFinal<-proc.time()-pt
  time4 <- rbind(time4,c(ptFinal[3])); print(c(ptFinal[3]))

  MISS1 <- rbind(MISS1,ls1[[1]])
  MISS2 <- rbind(MISS2,ls2[[1]])
  MISS3 <- rbind(MISS3,ls3[[1]])
  MISS4 <- rbind(MISS4,ls4[[1]])
  
  print('Current : ')
  print(round(tail(MISS1,1),4));
  print(round(tail(MISS2,1),4));
  print(round(tail(MISS3,1),4));
  print(round(tail(MISS4,1),4));
  
  print('Average :')
  print(round(colMeans(MISS1),4));
  print(round(colMeans(MISS2),4));
  print(round(colMeans(MISS3),4));
  print(round(colMeans(MISS4),4));
}

par(mfrow=c(1,3))
boxplot(MISS1[,1],MISS2[,1],MISS4[,1])
boxplot(MISS1[,2],MISS2[,2],MISS4[,2])
boxplot(MISS1[,3],MISS2[,3],MISS4[,3])





rm(list=ls())
############################ Library we use ####################################
library(mvtnorm) # for multivariate gaussian sampling
library(mnormt) # for fast matrix inverse compute : pd.solve (Used in MFVI)
library(rmutil) # for inverse gaussian sampling (Used in MCMC)
library(coda) # for geweke diagnostics (Used in MCMC)
library(Matrix) # for making block diagonal matrix
# Rcpp
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
sourceCpp('MEAN_FIELD_ALGORITHM.cpp') # for MFVI
sourceCpp('SVI_GIBBS.cpp') # Gibbs Sampling within SVI
# Other Variable Selection Methods :
library(sparsevb) # 2022_VI for BVS
library(SSLASSO) # The Spike and Slab LASSO
library(varbvs) # 2012_VI for BVS
library(EMVS) # EMVS
################################################################################


############################ Global Helper ftns ################################
tol <- 1e-200 # preventing from overflow {ex. log(0+tol), a/(0+tol), etc}
normalize <- function(x) { return ((x - min(x)) / (max(x) - min(x))) }
standardize <- function(x) { return ((x - mean(x)) / sd(x))}
# normalizing <- function(mat) mat/rowSums(mat)
mat_trace <- function(A) sum(diag(A))
ar1_cor <- function(n, rho) rho^abs(matrix(1:n-1,nrow=n,ncol=n,byrow=TRUE)-(1:n-1))
logit <- function(m) log((m+tol)/(1+m+tol))
log_sum_exp <- function(log_vec) log(sum(exp(log_vec-max(log_vec)))+tol) + max(log_vec)
################################################################################


########################## Sample Generating ftns ##############################
# Strength of coeffs : With arithmaetic sequence (with alternating sign)
# ex. true coeff = (0.1, -0.2, 0.3, -0.4, ... , -2.6, 2.8, -3.0)
# pairwise correlated structure 
sim_data_cor_all <- function(Q,corr,N,s,r,p){
  # Q : total number of covariates
  # corr : strength of correlation
  # N : total number of samples
  # s : strength of noise standard deviation
  # r : random seed
  # p : number of valid covariates among Q
  set.seed(r)
  
  ########### Generate coeffs
  int <- round(Q/p,0)
  idx_valid_x <- seq(1,Q,int)[1:p]
  true.coeff <- seq(0.1,3,length.out=p)
  
  minus_idx <- which(seq_len(p)%%2==1)
  true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[idx_valid_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[idx_valid_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=Q,ncol=Q)
  corr.x <- corr*J - corr*diag(Q) + diag(Q)
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff))
}
# AR1 correlated structure
sim_data_cor_ar1 <- function(Q,corr,N,s,r,p){
  # Q : total number of covariates
  # corr : strength of correlation
  # N : total number of samples
  # s : strength of noise standard deviation
  # r : random seed
  # p : number of valid covariates among Q
  set.seed(r)
  
  ########### Generate coeffs
  int <- round(Q/p,0)
  idx_valid_x <- seq(1,Q,int)[1:p]
  true.coeff <- seq(0.1,3,length.out=p)
  
  minus_idx <- which(seq_len(p)%%2==1)
  true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[idx_valid_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[idx_valid_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=Q,ncol=Q)
  corr.x <- ar1_cor(Q,corr)
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff))
}
# blockwise correlated structure
sim_data_cor_block <- function(Q,corr,N,s,r,p){
  set.seed(r)
  
  ########### Generate coeffs their index
  int <- round(Q/p,0)
  idx_valid_x <- seq(1,Q,int)[1:p]
  true.coeff <- seq(0.1,3,length.out=p)
  
  minus_idx <- which(seq_len(p)%%2==1)
  true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[idx_valid_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[idx_valid_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=int,ncol=int)
  corr.x.blk <- corr*J - corr*diag(1,nrow=int) + diag(1,nrow=int)
  corr.x <- as.matrix(bdiag(rep(list(corr.x.blk),(p))))
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff))
}

# Random coeff
# pairwise correlated structure 
sim_data_cor_all_random <- function(Q,corr,N,s,r,p){
  set.seed(r)
  
  ########### Generate coeffs
  int <- round(Q/p,0)
  idx_valid_x <- seq(1,Q,int)[1:p]
  true.coeff <- runif(p,0,3)
  
  minus_idx <- unique(sample.int(p,p,replace=T))
  true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[idx_valid_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[idx_valid_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=Q,ncol=Q)
  corr.x <- corr*J - corr*diag(Q) + diag(Q)
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff))
}
# AR1 correlated structure
sim_data_cor_ar1_random <- function(Q,corr,N,s,r,p){
  set.seed(r)
  
  ########### Generate coeffs
  int <- round(Q/p,0)
  idx_valid_x <- seq(1,Q,int)[1:p]
  true.coeff <- runif(p,0,3)
  
  minus_idx <- unique(sample.int(p,p,replace=T))
  true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[idx_valid_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[idx_valid_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=Q,ncol=Q)
  corr.x <- ar1_cor(Q,corr)
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff))
}
# blockwise correlated structure
sim_data_cor_block_random <- function(Q,corr,N,s,r,p){
  set.seed(r)
  
  ########### Generate coeffs
  int <- round(Q/p,0)
  idx_valid_x <- seq(1,Q,int)[1:p]
  true.coeff <- runif(p,0,3)
  
  minus_idx <- unique(sample.int(p,p,replace=T))
  true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[idx_valid_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[idx_valid_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=int,ncol=int)
  corr.x.blk <- corr*J - corr*diag(1,nrow=int) + diag(1,nrow=int)
  corr.x <- as.matrix(bdiag(rep(list(corr.x.blk),(p))))
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff))
}
########################## Sample Generating ftns ##############################


################################# Algorithms ###################################
# Structured VI (SVI) approximated with Gibbs Sampling
SVIG <- function(y,X,mem_true,coeff,lamb,S){
  # y : dependent var
  # X : design matrix
  # mem_true : true membership(nonzero or zero coeff) of each covariate
  # lamb : scale parameter of laplace distribution
  # S : number of samples for MCMC
  pt<-proc.time(); # for recording time
  burn <- round(S/10,1) # burn-in number for MCMC
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  
  N_train <- round(N*(4/5),0)
  N_test <- N-N_train
  W_test <- W[(N_train+1):N,]
  W <- W[1:N_train,]
  y <- y[1:N_train]
  y <- y-mean(y)
  
  # We calculate the below statistics (which are used for every iteration) 
  # in advance to reduce the computational burden :
  WTW <- t(W)%*%W
  yTy <- sum(y*y)
  yT1 <- sum(y)
  WTy <- t(W)%*%y
  
  #========== Iinitialization of  VI parameter values ==========#
  #===== error =====#
  inv_sigma2 <- 1/c(var(y)) # E_q[1/sigma^2]
  inv_sigma2.lt <- c() # list of E_q[1/sigma^2 | gamma], there are total of S samples
  # 1/sigma2|gamma ~ Inverse_Gamma(is2.a, is2.b|gamma)
  is2.a <- N_train/2
  is2.b.lt <- c() # list of is2.b|gamma
  
  #===== coefficients =====#
  Theta.mu <- matrix(0,nrow=Q,ncol=1) # E_q[Theta]
  Theta.mu.lt <- list() # list of E_q[Theta | sigma2, gamma], there are total of S samples
  Theta_ThetaT <- c() # E_q[Theta%*%t(Theta)]
  B.inv.lt <- list() # list of inverse of B, inverse of B : refer to paper
  
  #===== spike and slab prior =====#
  # gam : E_q[gamma]
  # gam_matrix : 2nd column : gam, 1st column : 1-gam
  gam_matrix <- matrix(1/2,Q,2)
  gam <- matrix(gam_matrix[,2],ncol=1)
  # cov_mem : initial gamma value for MCMC
  cov_mem <- as.numeric(sample.int(2,size=Q,replace=T,prob=c((Q-1)/Q,1/Q))==2)
  sample_mat <- matrix(NA,nrow=S,ncol=Q) # gamma samples, there are S number of gamma samples
  
  #===== rho (hyperparameter for spike and slab prior) =====#
  # rho ~ Beta(a0,b0)
  a0 <- 1
  b0 <- Q
  E_logit <- logit((a0/b0)) # E_q[logit(rho)]
  
  #===== Dispersion parameters for coefficients =====#
  lamb2_2 <- lamb^2/2
  p_tau <- 1/2; a_tau <- 2*lamb2_2 # parameters for Generalized Inverse Gaussian
  inv_tau2 <- rep(1/1,Q) # E_q[1/tau^2]
  Inv_V <- diag(inv_tau2,Q)
  
  #========== Optimization ==========#
  num_iters <- 20 # Maximum number of iterations for VI algorithm
  pred_hist <- c() # Storage of ELBO
  delta_H <- 1000 # Initial threshold (we stop the algorithm when delta_H < 1e-2)
  diff <- c() # number of misclassified memberships of covariates, store for every iteration
  Delta_H <- c() # store delta_H for every iteration
  
  for(idx in 1:num_iters){
    #===== gam update =====#
    # MCMC : randomized Gibbs sampling (GS)
    # rand_cov : for randomizing the order of covariates to GS update
    rand_cov <- sample(seq.int(from=1,to=Q),replace=F)-1
    LIST <- SVI_GIBBS(rand_cov, WTW+Inv_V, WTy, diag(Inv_V),
                      E_logit, is2.a, yTy,
                      cov_mem, S+burn)
    sample_mat <- LIST[[1]][(burn+1):(S+burn),]
    cov_mem <- sample_mat[S,] # This is the initial gamma value for the next VI iter MCMC
    Theta.mu.lt <- LIST[[2]][(burn+1):(S+burn),]
    B.inv.lt <- LIST[[3]][(burn+1):(S+burn)]
    inv_sigma2.lt <- LIST[[4]][(burn+1):(S+burn)]
    log_Q <- LIST[[5]][(burn+1):(S+burn)] # We will use log_Q to calculate ELBO
    is2.b.lt <- LIST[[6]][(burn+1):(S+burn)]
    
    inv_sigma2 <- sum(inv_sigma2.lt)/S
    log_sigma2 <- mean(log(is2.b.lt+tol)-digamma(is2.a))
    gam <- colMeans(sample_mat)
    gam_matrix[,2] <- gam; gam_matrix[,1] <- 1-gam
    
    #===== Theta update =====#
    Theta.cov.lt <- lapply(1:S, function(s){
      B.inv.lt[[s]]/inv_sigma2.lt[s]
    })
    Theta_ThetaT.lt <- lapply(1:S, function(s){
      Theta.cov.lt[[s]] + Theta.mu.lt[s,]%*%t(Theta.mu.lt[s,])
    })
    is2_Theta_ThetaT.lt <- lapply(1:S, function(s){
      inv_sigma2.lt[s]*Theta_ThetaT.lt[[s]]
    })
    
    Theta.mu.mt <- matrix(0,nrow=S,ncol=Q)
    is2_Theta_ThetaT <- matrix(0,nrow=Q,ncol=Q)
    
    for(s in 1:S){
      is2_Theta_ThetaT <- is2_Theta_ThetaT + is2_Theta_ThetaT.lt[[s]]
    }
    Theta.mu <- colMeans(Theta.mu.lt)
    is2_Theta_ThetaT <- is2_Theta_ThetaT/S
    
    # O_Theta_ThetaT : E_q[gamma%*%t(gamma)]%*%Theta_ThetaT : Used for ELBO
    O_Theta_ThetaT <- matrix(0,nrow=Q,ncol=Q)
    for(s in 1:S){
      O_Theta_ThetaT <- O_Theta_ThetaT + (sample_mat[s,]%*%t(sample_mat[s,]))*is2_Theta_ThetaT.lt[[s]]
    }
    O_Theta_ThetaT <- O_Theta_ThetaT/S
    
    #===== rho update =====#
    a <- a0+sum(gam_matrix[,2]); b <- b0+(Q-sum(gam_matrix[,2]))
    E_logit <- digamma(a)-digamma(b)
    
    #===== Dispersion update =====#
    b_tau <- diag(is2_Theta_ThetaT)
    inv_tau2 <- (sqrt(a_tau)*besselK(sqrt(a_tau*b_tau),p_tau+1))/
      (sqrt(b_tau)*besselK(sqrt(a_tau*b_tau),p_tau)) - ((2*p_tau)/b_tau)
    Inv_V <- diag(inv_tau2,Q)
    
    #===== ELBO update =====#
    tau2 <- (sqrt(b_tau)*besselK(sqrt(a_tau*b_tau),p_tau+1))/
      (sqrt(a_tau)*besselK(sqrt(a_tau*b_tau),p_tau))
    var_tau2 <- (b_tau/a_tau)*
      (besselK(sqrt(a_tau*b_tau),p_tau+2)/besselK(sqrt(a_tau*b_tau),p_tau) -
         (besselK(sqrt(a_tau*b_tau),p_tau+1)/besselK(sqrt(a_tau*b_tau),p_tau))^2)
    is2_Theta.mu <- colMeans(inv_sigma2.lt*Theta.mu.lt)
    TMP <- sum(WTW*O_Theta_ThetaT) + sum(Inv_V*is2_Theta_ThetaT)
    log_1_rho <- digamma(b)-digamma(a+b)
    log_Z <- log_sum_exp(log_Q)
    log_tau2 <- log(tau2) - var_tau2/(2*tau2^2)
    
    # Store ELBO value for every iteration
    pred_hist[idx] <-
      sum((1/4-lamb2_2)*tau2) + lamb2_2*sum(inv_tau2) + sum((1/2-is2_Theta.mu)*log_tau2) +
      -0.5*inv_sigma2*yTy + t(WTy)%*%is2_Theta.mu - 0.5*TMP +
      Q*log_1_rho + log_Z +
      log(2*lamb2_2)/2*sum(b_tau) + sum(log(besselK(sqrt(lamb2_2),b_tau)))
    
    # For plotting gam, sample_mat, and pred_hist for every iteration
    
    # par(mfrow=c(3,1))
    # plot(c(gam))
    # M1 <- t(sample_mat)
    # image(1:nrow(M1),1:ncol(M1),M1,col=min(M1):max(M1))
    # plot(pred_hist,type='l')
    
    if(idx>3){
      delta_H <- abs((tol+pred_hist[idx]+pred_hist[idx-1]) - 
                       (tol+pred_hist[idx-1]+pred_hist[idx-2]))/
        abs(tol+pred_hist[idx-1]+pred_hist[idx-2])
      Delta_H[idx-3] <- delta_H
    }
    if(delta_H<=1e-2){
      rand_cov <- sample(seq.int(from=1,to=Q),replace=F)-1
      LIST <- SVI_GIBBS(rand_cov, WTW+Inv_V, WTy, diag(Inv_V),
                        E_logit, is2.a, yTy,
                        cov_mem, S+burn)
      sample_mat <- LIST[[1]][(burn+1):(S+burn),]
      Theta.mu.lt <- LIST[[2]][(burn+1):(S+burn),]
      log_Q <- LIST[[5]][(burn+1):(S+burn)]
      log_Z <- log_sum_exp(log_Q)
      
      gam <- colMeans(sample_mat)
      gam_matrix[,2] <- gam; gam_matrix[,1] <- 1-gam
      Theta.mu <- colMeans(Theta.mu.lt)
      
      mem_est <- max.col(gam_matrix,'first')
      diff[idx] <- sum(abs(mem_true-mem_est))
      ptFinal<-proc.time()-pt
      break
    }
    mem_est <- max.col(gam_matrix,'first')
    diff[idx] <- sum(abs(mem_true-mem_est))
    ptFinal<-proc.time()-pt
  }
  print(diff)
  resid <- W_test%*%coeff-(W_test%*%(Theta.mu))
  mse <- norm(resid,type='2')/sqrt(N_test)
  
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
}
# Mean Field VI
MFVI <- function(y,X,mem_true,coeff,lamb){
  pt<-proc.time()
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  
  N_train <- round(N*(4/5),0)
  N_test <- N-N_train
  W_test <- W[(N_train+1):N,]
  W <- W[1:N_train,]
  y <- y[1:N_train]
  
  WTW <- t(W)%*%W
  y <- y-mean(y)
  WTy <- t(W)%*%y
  
  #========== Iinitial VI parameter values ==========#
  # ===== error =====#
  inv_sigma2 <- 1/c(var(y))
  
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
  lamb2_2 <- lamb^2/2
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1/100,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #===== Deterministic Annealing(if needed)=====#
  xi <- 0.001
  
  #========== Optimization ==========#
  j1 <- matrix(1,ncol=1,nrow=N)
  num_iters <- 100
  pred_hist <- c()
  delta_H <- 1000
  diff <- c()
  
  for(idx in 1:num_iters){
    #===== logVk update =====# ***********
    A <- A0+sum(gam_matrix[,2]); B <- B0+(Q-sum(gam_matrix[,2]))
    E_logit <- digamma(A)-digamma(B)
    
    #===== Theta update =====#
    Theta.cov <- pd.solve(inv_sigma2*(WTWO+Inv_V))
    Theta.mu <- inv_sigma2*gam*WTy; Theta.mu <- Theta.cov%*%Theta.mu
    Theta_ThetaT <- Theta.cov + Theta.mu%*%t(Theta.mu)
    theta.p2 <- diag(Theta_ThetaT)
    
    #===== Dispersion paramter =====#
    b_tau <- inv_sigma2*theta.p2
    inv_tau2 <- (sqrt(a_tau)*besselK(sqrt(a_tau*b_tau),p_tau+1))/
      (sqrt(b_tau)*besselK(sqrt(a_tau*b_tau),p_tau)) - ((2*p_tau)/b_tau)
    Inv_V <- diag(inv_tau2,Q)
    
    #===== Inv_sigma2 update =====#
    TMP <- sum((WTWO+Inv_V)*Theta_ThetaT)
    is2.a <- (N_train+Q-1)/2
    is2.b <- c(t(y)%*%y - 2*t(WTy)%*%(gam*Theta.mu) + TMP)/2
    inv_sigma2 <- is2.a/is2.b
    
    #===== Gam update =====#
    PSI <- xi*(-0.5*inv_sigma2*(WTW*Theta_ThetaT)); psi <- diag(PSI); diag(PSI) <- 0
    psi <- xi*(E_logit + inv_sigma2*(c(t(WTy))*c(Theta.mu))) + psi
    gam_matrix <- MEAN_FIELD_ALGORITHM(Q,psi,PSI,gam_matrix)
    gam <- matrix(gam_matrix[,2],ncol=1)
    Gam <- diag(c(gam),ncol=Q)
    Gam_2 <- diag(c(gam^2),ncol=Q)
    Omega <- gam%*%t(gam) - Gam_2 + Gam
    WTWO <- WTW*Omega
    
    #==== xi update ====#
    xi <- xi*1.2
    if(xi>=1) xi<-1
    
    mem_est <- max.col(gam_matrix,'first')
    prob<- apply(gam_matrix,1,max)
    diff[idx] <- sum(abs(mem_true-mem_est))
    
    # I did a mistake...! I didn't use ELBO as a stopping criterion...
    pred_hist[idx] <- norm(W%*%(gam*Theta.mu),'2')/sqrt(N_train)
    if(idx>1){
      delta_H <- abs(pred_hist[idx]-pred_hist[idx-1])/abs(pred_hist[idx-1])
    }
    if(delta_H<=1e-5) break
  }
  ptFinal<-proc.time()-pt
  print(diff)
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  
  resid <- W_test%*%coeff-(W_test%*%(gam*Theta.mu))
  mse <- norm(resid,type='2')/sqrt(N_test)
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
}
# MCMC
MCMC <- function(y,X,mem_true,coeff,total){
  # total : limit the number of MCMC
  pt<-proc.time()
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  
  N_train <- round(N*(4/5),0)
  N_test <- N-N_train
  W_test <- W[(N_train+1):N,]
  W <- W[1:N_train,]
  y <- y[1:N_train]
  
  WTW <- t(W)%*%W
  y <- y-mean(y)
  yTy <- sum(y*y)
  yT1 <- sum(y)
  WTy <- t(W)%*%y
  
  #========== Iinitial VI parameter values ==========#
  #===== error =====#
  inv_Sigma2 <- 1/c(var(y))
  inv_sigma2.lt <- c()
  is2.b.lt <- c()
  is2.a <- N_train/2
  
  #===== coefficients =====#
  Theta.mu <- c()
  Theta_ThetaT <- c()
  Theta.lt <- matrix(0,nrow=0,ncol=Q)
  Theta_ThetaT.lt <- list()
  B.inv.lt <- list()
  
  #===== spike and slab prior =====#
  gam <- rep(0,Q)
  gam_matrix <- matrix(0,nrow=0,ncol=Q)
  sample_mat <- matrix(0,nrow=0,ncol=Q)
  
  #===== rho =====#
  a0 <- 1
  b0 <- Q
  rho <- a0/b0
  rho.lt <- c()
  logit_rho <- logit(rho)
  
  #===== Dispersion parameters =====#
  lamb <- 1/1; lamb2_2 <- lamb^2/2
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2.lt <- matrix(0,nrow=0,ncol=Q)
  inv_tau2 <- rep(1/1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #========== Optimization ==========#
  diff <- c()
  under_p_val <- Q # for stopping criterion
  
  for(idx in 1:total){
    #===== gam sampling =====#
    rand_cov <- sample(seq.int(from=1,to=Q),replace=F)-1
    LIST <- SVI_GIBBS(rand_cov, WTW+Inv_V, WTy, diag(Inv_V),
                      logit_rho, is2.a, yTy,
                      gam, 1)
    sample_mat <- rbind(sample_mat,LIST[[1]])
    gam <- tail(sample_mat,1)
    Theta.mu <- LIST[[2]]
    B.inv <- LIST[[3]][[1]]
    is2.b <- LIST[[6]]
    
    #===== rho sampling =====#
    a <- a0+sum(gam); b <- b0+(Q-sum(gam))
    rho <- rbeta(n=1,a,b)
    rho.lt <- c(rho.lt,rho)
    logit_rho <- logit(rho)
    
    #===== sigma2 sampling =====#
    inv_sigma2.lt <- c(inv_sigma2.lt,rgamma(1,shape=is2.a,scale=1/is2.b))
    inv_sigma2 <- tail(inv_sigma2.lt,1)
    
    #===== Theta sampling =====#
    Theta.cov <- 1/inv_sigma2*B.inv
    Theta.lt <- rbind(Theta.lt,rmvnorm(n=1,mean=Theta.mu,sigma=Theta.cov))
    Theta <- tail(Theta.lt,1)
    
    #===== tau2 sampling =====#
    b_tau <- inv_sigma2*Theta^2
    for(j in 1:Q){
      inv_tau2[j] <- rinvgauss(n=1,
                               m=sqrt((a_tau+tol)/(b_tau[j]+tol)),
                               s=a_tau)
    }
    inv_tau2.lt <- rbind(inv_tau2.lt,inv_tau2)
    Inv_V <- diag(inv_tau2,Q)
    
    mem_est <- as.integer(colMeans(sample_mat)>0.5)+1
    diff[idx] <- sum(abs(mem_true-mem_est))
    
    # Check the convergence for every 50 iteration
    # because this takes for a while for large number of covariates
    # which means if we check the convergence for every iteraion, MCMC will take much more time than we expect
    # but you can freely modify this criterion
    if(idx%%50==0){
      p_val <- pnorm(abs(geweke.diag(Theta.lt)$z),lower.tail=FALSE)*2
      under_p_val <- sum(p_val < 0.001) # I recommend to set p_val threshold as small as possible
      # because the MCMC algorithm won't converge if it is set higher
      cat(idx,diff[idx-1],under_p_val,'\n')
    }
    if(under_p_val==0) break
  }
  ptFinal<-proc.time()-pt
  print(ptFinal[3])
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  
  Theta_hat <- sample_mat*Theta.lt
  L <- dim(Theta_hat)[1]; L_eff <- round(L*0.5,0) # only use the last L_eff number of samples
  resid <- W_test%*%coeff-W_test%*%colMeans(tail(Theta_hat,L_eff))
  mse <- norm(resid,type='2')/sqrt(N_test)
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
}
# 2022_VI for BVS
SPVB <- function(y,X,mem_true,coeff){
  tryCatch(
    {
      pt<-proc.time()
      Q <- dim(X)[2]; N <- dim(X)[1]
      y <- y-mean(y)
      W <- apply(X,2,standardize)
      N_train <- round(N*(4/5),0)
      N_test <- N-N_train
      W_test <- W[(N_train+1):N,]
      W <- W[1:N_train,]
      y <- y[1:N_train]
      
      fit <- svb.fit(W,y,family='linear',intercept=F,slab='laplace')
      ptFinal<-proc.time()-pt
      message('succeed')
      mem_est <- as.numeric(fit$gamma>0.5)+1
      mem_true.x <- mem_true[1:Q]
      mem2.x <- mem_est[which(mem_true.x==2)]
      mem1.x <- mem_est[which(mem_true.x==1)]
      len2.x <- sum(mem_true.x==2)
      miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                       sum(mem1.x != rep(1,Q-len2.x)),
                       sum(mem_true.x != mem_est))
      gam <- fit$gamma
      Theta.mu <- fit$mu
      
      resid <- W_test%*%coeff-(W_test%*%(gam*Theta.mu))
      mse <- norm(resid,type='2')/sqrt(N_test)
      
      # F1 score
      TP <- sum((mem_est+mem_true)==4)
      FP <- sum((mem_est-mem_true)==1)
      FN <- sum((mem_est-mem_true)==-1)
      precision <- TP/(TP+FP)
      recall <- TP/(TP+FN)
      F1 <- (2*precision*recall)/(precision+recall)
      return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
    },
    error = function(e){
      message('Caught an error!')
      print(e)
      return(list(c(NaN,NaN,NaN),NaN,NaN,NaN,NaN))
    }
  )    
}
# 2012_VI for BVS
varbvss <- function(y,X,mem_true,coeff){
  Q <- dim(X)[2]; N <- dim(X)[1]
  y <- y-mean(y)
  W <- apply(X,2,standardize)
  N_train <- round(N*(4/5),0)
  N_test <- N-N_train
  W_test <- W[(N_train+1):N,]
  W <- W[1:N_train,]
  y <- y[1:N_train]
  
  pt<-proc.time()
  fit <- varbvs(X=W,Z=NULL,y=y,family='gaussian',tol=1e-4,maxiter=1e4,verbose=F)
  ptFinal<-proc.time()-pt
  mem_est <- as.integer(fit$pip>0.5)+1
  
  gam <- as.integer(mem_est==2)
  
  y_hat <- predict(fit, X=W_test, Z = NULL,family='gaussian')
  resid <- W_test%*%coeff-y_hat
  mse <- norm(resid,type='2')/sqrt(N_test)
  
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
}
# The Spike and Slab LASSO
SSLASSOO <- function(y,X,mem_true,coeff,X_test){
  Q <- dim(X)[2]; N <- dim(X)[1]
  y <- y-mean(y)
  W <- apply(X,2,standardize)
  N_train <- round(N*(4/5),0)
  N_test <- N-N_train
  W_test <- W[(N_train+1):N,]
  W <- W[1:N_train,]
  y <- y[1:N_train]
  
  lambda0 <- seq(0.1,Q,1)
  pt<-proc.time()
  fit <- SSLASSO(W,y,penalty='adaptive',variance='unknown',
                 lambda1=0.1,lambda0=lambda0,theta=1/Q,
                 a=1,b=Q)
  ptFinal<-proc.time()-pt
  linear_idx <- fit$model
  mse_lt <- rep(10000,length(lambda0))
  for(idx in 1:length(lambda0)){
    check_sum <- length(linear_idx)
    linear_idx_tmp <- which(fit$select[,idx]==1)
    if(sum(linear_idx==linear_idx_tmp) == check_sum){
      Theta.mu <- fit$beta[,idx]; intercept <- fit$intercept[,idx]
      resid <- W_test%*%coeff-(intercept + W_test%*%(Theta.mu))
      mse <- norm(resid,type='2')/sqrt(N_test)
      mse_lt[idx] <- mse
    }
  }
  min_idx <- which.min(mse_lt)
  mse <- mse_lt[min_idx]
  
  Theta.mu <- fit$beta[,min_idx]; intercept <- fit$intercept[,min_idx]
  mem_est <- rep(1,Q); mem_est[linear_idx] <- 2
  gam <- as.integer(mem_est==2)
  
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
}
# EMVS
EMVSS <- function(y,X,mem_true,coeff,X_test){
  Q <- dim(X)[2]; N <- dim(X)[1]
  y <- y-mean(y)
  W <- apply(X,2,standardize)
  N_train <- round(N*(4/5),0)
  N_test <- N-N_train
  W_test <- W[(N_train+1):N,]
  W <- W[1:N_train,]
  y <- y[1:N_train]
  
  v0 = seq(0.01, 0.1, length.out = 10)
  v1 = 1
  a = 1
  b = Q
  sigma_init = 1
  
  pt<-proc.time()
  invisible(capture.output(fit <- EMVS(y,W,v0=v0,v1=v1,type=c('betabinomial'),independent=F,
                                       sigma_init=sigma_init,a=a,b=b)))
  ptFinal<-proc.time()-pt
  idxes <- which(c(EMVSbest(fit)$log_g_function==fit$log_g_function)==T)
  
  mse_lt <- c()
  for(idx in 1:length(idxes)){
    i <- idxes[idx]
    Theta.mu <- fit$betas[i,]
    resid <- W_test%*%coeff-(W_test%*%(Theta.mu))
    mse <- norm(resid,type='2')/sqrt(N_test)
    mse_lt[idx] <- mse
  }
  min_idx <- which.min(mse_lt)
  mse <- mse_lt[min_idx]
  gam <- fit$prob_inclusion[min_idx,]
  
  mem_est <- as.numeric(gam>0.5)+1
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
}
################################################################################


############################ Experiment settings ###############################
# simulation setting with 'Corr_struct', 'SIM', 'Q', 'corr', 'N', 'sigma' ,'p', 'lamb', 'S'
# Corr_struct : correlation structure (all, ar1, block)
# SIM : total number of simulations
# corr : strength of correlation
# N : number of data
# sigma : strenght of sqrt(noise variance)
# p : number of valid covariates
# lamb : strength of laplace distribution (only for MFVI and SVIG)
# S : number of MCMC samples used in SVIG

# We do simulation with random seed from 'start' to 'SIM'
# If we are to store the results for the first time, set 'load'==FALSE
# If we already have the simulation result and we want to do more simulation with the same setting,
# set 'load'==TRUE to load the simulation result and this will append the additional simulation results
simulation <- function(Corr_struct,SIM,Q,corr,N,sigma,p,lamb,S,load,start){
  # storage
  MISS1 <- matrix(0,nrow=0,ncol=3)
  MISS2 <- matrix(0,nrow=0,ncol=3)
  MISS3 <- matrix(0,nrow=0,ncol=3)
  MISS4 <- matrix(0,nrow=0,ncol=3)
  MISS5 <- matrix(0,nrow=0,ncol=3)
  MISS6 <- matrix(0,nrow=0,ncol=3)
  MISS7 <- matrix(0,nrow=0,ncol=3)
  
  time1 <- time2 <- time3 <- time4 <- time5 <- time6 <- time7 <- c()
  mse1 <- mse2 <- mse3 <- mse4 <- mse5 <- mse6 <- mse7 <- c()
  F1 <- F2 <- F3 <- F4 <- F5 <- F6 <- F7 <- c()
  
  List <- list()
  
  main_text <- paste0('P',Q,', ',
                      'N',N,', ',
                      'p',p,', ',
                      Corr_struct,' ',corr)
  i <- 1
  if(load==T){
    i <- start
  }
  
  for(r in i:SIM){
    print('SIMULATION : '); print(r)
    ################################################################
    ######################## generate sample #######################
    if(Corr_struct=='all') List <- sim_data_cor_all(Q,corr,N,sigma,r,p)
    if(Corr_struct=='all_random') List <- sim_data_cor_all_random(Q,corr,N,sigma,r,p)
    if(Corr_struct=='block') List <- sim_data_cor_block(Q,corr,N,sigma,r,p)
    if(Corr_struct=='block_random') List <- sim_data_cor_block_random(Q,corr,N,sigma,r,p)
    if(Corr_struct=='ar1') List <- sim_data_cor_ar1(Q,corr,N,sigma,r,p)
    if(Corr_struct=='ar1_random') List <- sim_data_cor_ar1_random(Q,corr,N,sigma,r,p)
    
    y <- List[[1]]; X <- List[[2]]
    data <- List[[3]]; mem_true <- List[[4]]; coeff <- List[[5]]
    ################################################################
    ################################################################
    
    ################################################################
    ########################## Simulations #########################
    print('SVIG'); ls1 <- SVIG(y,X,mem_true,coeff,lamb,S)
    time1 <- c(time1,ls1[[5]])
    mse1 <- c(mse1,ls1[[3]])
    F1 <- c(F1,ls1[[4]])
    MISS1 <- rbind(MISS1,ls1[[1]])

    print('MFVI'); ls2 <- MFVI(y,X,mem_true,coeff,lamb)
    time2 <- c(time2,ls2[[5]])
    mse2 <- c(mse2,ls2[[3]])
    F2 <- c(F2,ls2[[4]])
    MISS2 <- rbind(MISS2,ls2[[1]])
    
    print('EMVS'); ls3 <- EMVSS(y,X,mem_true,coeff)
    time3 <- c(time3,ls3[[5]])
    mse3 <- c(mse3,ls3[[3]])
    F3 <- c(F3,ls3[[4]])
    MISS3 <- rbind(MISS3,ls3[[1]])
    
    print('varbvs'); ls4 <- varbvss(y,X,mem_true,coeff)
    time4 <- c(time4,ls4[[5]])
    mse4 <- c(mse4,ls4[[3]])
    F4 <- c(F4,ls4[[4]])
    MISS4 <- rbind(MISS4,ls4[[1]])

    print('SSLASSO'); ls5 <- SSLASSOO(y,X,mem_true,coeff)
    time5 <- c(time5,ls5[[5]])
    mse5 <- c(mse5,ls5[[3]])
    F5 <- c(F5,ls5[[4]])
    MISS5 <- rbind(MISS5,ls5[[1]])

    print('SPVB'); ls6 <- SPVB(y,X,mem_true,coeff)
    time6 <- c(time6,ls6[[5]])
    mse6 <- c(mse6,ls6[[3]])
    F6 <- c(F6,ls6[[4]])
    MISS6 <- rbind(MISS6,ls6[[1]])

    print('MCMC'); ls7<- MCMC(y,X,mem_true,coeff,1000)
    time7 <- c(time7,ls7[[5]])
    mse7 <- c(mse7,ls7[[3]])
    F7 <- c(F7,ls7[[4]])
    MISS7 <- rbind(MISS7,ls7[[1]])

    
    ################################################################
    ####################### Printout results########################
    print('Current FP, FN, and total missclassification : ')
    print(round(tail(MISS1,1),3));
    print(round(tail(MISS2,1),3));
    print(round(tail(MISS3,1),3));
    print(round(tail(MISS4,1),3));
    print(round(tail(MISS5,1),3));
    print(round(tail(MISS6,1),3));
    print(round(tail(MISS7,1),3));
    
    print('Average FP, FN, and total missclassification : ')
    print(round(colMeans(MISS1,na.rm=T),3));
    print(round(colMeans(MISS2,na.rm=T),3));
    print(round(colMeans(MISS3,na.rm=T),3));
    print(round(colMeans(MISS4,na.rm=T),3));
    print(round(colMeans(MISS5,na.rm=T),3));
    print(round(colMeans(MISS6,na.rm=T),3));
    print(round(colMeans(MISS7,na.rm=T),3));
    
    print('Average rmse : ')
    print(round(mean(mse1,na.rm=T),3));
    print(round(mean(mse2,na.rm=T),3));
    print(round(mean(mse3,na.rm=T),3));
    print(round(mean(mse4,na.rm=T),3));
    print(round(mean(mse5,na.rm=T),3));
    print(round(mean(mse6,na.rm=T),3));
    print(round(mean(mse7,na.rm=T),3));
    ################################################################
    ################################################################
    par(mfrow=c(3,1))
    boxplot(F2,F3,F4,F5,F6,F1,F7,
            main=main_text,
            names=c('MFVI','EMVS','varbvs','SSLASSO','SPVB','SVIG','MCMC'),
            ylab='F1-score',cex.lab=1.5)
    boxplot(mse2,mse3,mse4,mse5,mse6,mse1,mse7,
            names=c('MFVI','EMVS','varbvs','SSLASSO','SPVB','SVIG','MCMC'),
            ylab='RMSE',cex.lab=1.5)
    time_box <- boxplot(time2,time3,time4,time5,time6,time1,time7,
                        main='Time',
                        names=c('MFVI','EMVS','varbvs','SSLASSO','SPVB','SVIG','MCMC'))
  }
  
  F1_LIST = list(F1,F2,F3,F4,F5,F6,F7)
  RMSE_LIST = list(mse1,mse2,mse3,mse4,mse5,mse6,mse7)
  TIME_LIST = list(time1,time2,time3,time4,time5,time6,time7)
  
  storage_text <- paste0('P',Q,'_',
                         'N',N,'_',
                         'p',p,'_',
                         Corr_struct,'_',corr)
  storage_place <- paste0("C:\\Users\\Public\\",storage_text,'.RData')
  
  if(load==F){
    save_file <- list(F1_list = F1_LIST,
                      RMSE_list = RMSE_LIST,
                      TIME_list = TIME_LIST)
    save(save_file, file=storage_place)
  }else{
    load(file=storage_place)
    M <- length(save_file$F1_list)
    for(m in 1:M){
      save_file$F1_list[[m]] <- c(save_file$F1_list[[m]],F1_LIST[[m]])
      save_file$RMSE_list[[m]] <- c(save_file$RMSE_list[[m]],RMSE_LIST[[m]])
      save_file$TIME_list[[m]] <- c(save_file$TIME_list[[m]],TIME_LIST[[m]])
    }
    save(save_file, file=storage_place)
  }
  
  store_image <- paste0(storage_text,'.png')
  png(store_image, width = 600, height = 800)
  par(mfrow=c(3,1))
  boxplot(F2,F3,F4,F5,F6,F1,F7,
          main=main_text,
          names=c('MFVI','EMVS','varbvs','SSLASSO','SPVB','SVIG','MCMC'),
          ylab='F1-score',cex.lab=1.5)
  boxplot(mse2,mse3,mse4,mse5,mse6,mse1,mse7,
          names=c('MFVI','EMVS','varbvs','SSLASSO','SPVB','SVIG','MCMC'),
          ylab='RMSE',cex.lab=1.5)
  time_box <- boxplot(time2,time3,time4,time5,time6,time1,time7,
                      main='Time',
                      names=c('MFVI','EMVS','varbvs','SSLASSO','SPVB','SVIG','MCMC'))
  # dev.copy(png, filename=store_image)
  # dev.print(width = 6, height = 6)
  dev.off()
  # graphics.off()
}
################################################################################


############################### Do Simulations! ################################
N_train <- 70; N <- round(N_train*(5/4),0)
Q<-200; p<-10; sigma<-sqrt(1)
SIM<-10
Corr_struct_list <- c('ar1','block','all')
# Corr_struct_list <- c('ar1_random','block_random','all_random')
S <- 10
lamb <- 1
for(c in 2:1){
  corr <- 0.4*c
  print(corr)
  for(cs in c(2,1,3)){
    Corr_struct <- Corr_struct_list[cs]
    print(Corr_struct)
    simulation(Corr_struct,SIM,Q,corr,N,sigma,p,lamb,S,F,1)
  }
}
################################################################################













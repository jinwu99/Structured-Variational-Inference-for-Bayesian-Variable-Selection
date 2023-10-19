rm(list=ls())
library(Matrix)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(mnormt) # for pd.solve
library(sparsevb)
sourceCpp('MEAN_FIELD_ALGORITHM.cpp')
sourceCpp('SVI_GIBBS.cpp')
sourceCpp('SVI_EXACT.cpp')

################################################
############## global functions ################
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

################################################
############ Sample generating ftns ############
# all the pairwise correlations are corr
sim_data_cor_all <- function(Q,corr,N,s,r,p){
  set.seed(r)
  
  ########### Generate coeffs their index
  int <- round(Q/p,0)
  idx_valid_x <- seq(1,Q,int)[1:p]
  true.coeff <- seq(2,2,length.out=p)
  
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
# corr with AR(1)
sim_data_cor_ar1 <- function(Q,corr,N,s,r,p){
  set.seed(r)
  
  ########### Generate coeffs their index
  int <- round(Q/p,0)
  idx_valid_x <- seq(1,Q,int)[1:p]
  true.coeff <- seq(1,2,length.out=p)
  
  # minus_idx <- unique(sample.int(p,p,replace=T))
  # true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[idx_valid_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[idx_valid_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=Q,ncol=Q)
  corr.x <- ar1_cor(Q,corr)
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  X_test <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff,X_test))
}
# blockwise correlation corr
sim_data_cor_block <- function(Q,corr,N,s,r,p){
  set.seed(r)
  
  ########### Generate coeffs their index
  int <- round(Q/p,0)
  idx_valid_x <- seq(1,Q,int)[1:p]
  true.coeff <- seq(1,2,length.out=p)
  
  # minus_idx <- unique(sample.int(p,p,replace=T))
  # true.coeff[minus_idx] <- -1*true.coeff[minus_idx]
  
  coeff <- rep(0,Q)
  coeff[idx_valid_x] <- true.coeff
  
  mem_true_x <- rep(1,Q)
  mem_true_x[idx_valid_x] <- 2
  
  ########### Covariance structure
  J <- matrix(1,nrow=int,ncol=int)
  corr.x.blk <- corr*J - corr*diag(1,nrow=int) + diag(1,nrow=int)
  corr.x <- as.matrix(bdiag(rep(list(corr.x.blk),(p))))
  X <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  X_test <- rmvnorm(n=N, mean=rep(0,Q), sigma=corr.x)
  
  ########### making data
  intercept <- 1
  y <- intercept + X%*%coeff + rnorm(N,mean=0,sd=s)
  
  data <- data.frame(y,X)
  return(list(y,X,data,mem_true_x,coeff,X_test))
}
################################################
################################################

################################################
################## Algorithms ##################
# Mean Field VI
MFVI <- function(y,X,mem_true,coeff,X_test,xi){
  pt<-proc.time()
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  W_test <- apply(X_test,2,standardize)
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
  lamb <- 1/1; lamb2_2 <- lamb^2/2
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1/1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #===== Deterministic Annealing(if needed)=====#
  # xi <- 1
  # xi <- 0.01
  
  #========== Optimization ==========#
  j1 <- matrix(1,ncol=1,nrow=N)
  num_iters <- 100
  pred_hist <- c()
  delta_H <- 1000
  diff <- c()
  
  for(idx in 1:num_iters){
    Theta.cov <- pd.solve(inv_sigma2*(WTWO+Inv_V))
    Theta.mu <- inv_sigma2*gam*WTy; Theta.mu <- Theta.cov%*%Theta.mu
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

    #===== Gam update =====#
    PSI <- xi*(-0.5*inv_sigma2*(WTW*Theta_ThetaT)); psi <- diag(PSI); diag(PSI) <- 0
    psi <- xi*(E_logit + inv_sigma2*(c(t(WTy))*c(Theta.mu))) + psi
    gam_matrix <- MEAN_FIELD_ALGORITHM(Q,psi,PSI,gam_matrix)
    gam <- matrix(gam_matrix[,2],ncol=1)
    Gam <- diag(c(gam),ncol=Q)
    Gam_2 <- diag(c(gam^2),ncol=Q)
    Omega <- gam%*%t(gam) - Gam_2 + Gam
    WTWO <- WTW*Omega
    
    #===== Inv_sigma2 update =====#
    TMP <- sum((WTWO+Inv_V)*Theta_ThetaT)
    is2.a <- (N+Q-1)/2
    is2.b <- c(t(y)%*%y - 2*t(WTy)%*%(gam*Theta.mu) + TMP)/2
    inv_sigma2 <- is2.a/is2.b
    
    #==== xi update ====#
    xi <- xi*1.2
    if(xi>=1) xi<-1

    mem_est <- max.col(gam_matrix,'first')
    prob<- apply(gam_matrix,1,max)
    diff[idx] <- sum(abs(mem_true-mem_est))
    
    # I did a mistake...! I didn't use ELBO as a stopping criterion...
    pred_hist[idx] <- norm(W%*%(gam*Theta.mu),'2')/sqrt(N)
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
  mse <- norm(resid,type='2')/sqrt(N)
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
}
# Structured VI with Exact computation
SVIE <- function(y,X,mem_true,coeff,X_test){
  pt<-proc.time();
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  W_test <- apply(X_test,2,standardize)
  WTW <- t(W)%*%W
  y <- y-mean(y)
  yTy <- sum(y*y)
  yT1 <- sum(y)
  WTy <- t(W)%*%y
  
  #========== Iinitial VI parameter values ==========#
  #===== error =====#
  inv_sigma2.lt <- c()
  inv_sigma2 <- 1/c(var(y))
  is2.b.lt <- c()
  is2.a <- N/2
  
  #===== coefficients =====#
  Theta.mu <- matrix(0,nrow=Q,ncol=1)
  Theta.mu.lt <- list()
  Theta_ThetaT <- c()
  B.inv.lt <- list()
  
  #===== spike and slab prior =====#
  gam_matrix <- matrix(1/2,Q,2)
  gam <- matrix(gam_matrix[,2],ncol=1)
  
  #===== rho =====#
  a0 <- 1
  b0 <- Q
  E_logit <- logit((a0/b0))
  
  #===== Dispersion parameters =====#
  lamb <- 1/1; lamb2_2 <- lamb^2/2
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1/1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #========== Optimization ==========#
  num_iters <- 20
  pred_hist <- c()
  delta_H <- 1000
  diff <- c()
  Delta_H <- c()
  
  # Comb : all the 2^p combination that gamma can have
  Comb <- as.matrix(expand.grid(lapply(1:Q,function(x) c(0,1))))
  G <- dim(Comb)[1]

  for(idx in 1:num_iters){
    #===== gam update =====#
    LIST <- SVI_EXACT(Comb, WTW+Inv_V, WTy, diag(Inv_V),
                      E_logit, is2.a, yTy)
    log_Q <- LIST[[1]]
    Theta.mu.lt <- LIST[[2]]
    B.inv.lt <- LIST[[3]]
    inv_sigma2.lt <- LIST[[4]]
    inv_sigma2 <- mean(inv_sigma2.lt)
    
    Q_gam <- exp(log_Q - max(log_Q))
    q_gam <- Q_gam/sum(Q_gam)
    gam <- colSums(q_gam*Comb)
    gam_matrix[,2] <- gam; gam_matrix[,1] <- 1-gam
    
    #===== Theta update =====#
    Theta.cov.lt <- lapply(1:G, function(s){
      B.inv.lt[[s]]/inv_sigma2.lt[s]
    })
    Theta_ThetaT.lt <- lapply(1:G, function(s){
      Theta.cov.lt[[s]] + Theta.mu.lt[s,]%*%t(Theta.mu.lt[s,])
    })
    is2_Theta_ThetaT.lt <- lapply(1:G, function(s){
      inv_sigma2.lt[s]*Theta_ThetaT.lt[[s]]
    })
    
    Theta.mu.mt <- matrix(0,nrow=G,ncol=Q)
    is2_Theta_ThetaT <- matrix(0,nrow=Q,ncol=Q)
    
    for(s in 1:G){
      is2_Theta_ThetaT <- is2_Theta_ThetaT + q_gam[s]*is2_Theta_ThetaT.lt[[s]]
    }
    Theta.mu <- colSums(q_gam*Theta.mu.lt)
    
    O_Theta_ThetaT <- matrix(0,nrow=Q,ncol=Q)
    for(s in 1:G){
      O_Theta_ThetaT <- O_Theta_ThetaT + q_gam[s]*(Comb[s,]%*%t(Comb[s,]))*is2_Theta_ThetaT.lt[[s]]
    }
    
    #===== logVk update =====#
    a <- a0+sum(gam_matrix[,2]); b <- b0+(Q-sum(gam_matrix[,2]))
    E_logit <- digamma(a)-digamma(b)
    
    #===== Dispersion paramter =====#
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
    is2_Theta.mu <- colMeans(q_gam*inv_sigma2.lt*Theta.mu.lt)
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
    
    if(idx>3){
      delta_H <- abs((tol+pred_hist[idx]+pred_hist[idx-1]) - 
                       (tol+pred_hist[idx-1]+pred_hist[idx-2]))/
        abs(tol+pred_hist[idx-1]+pred_hist[idx-2])
      Delta_H[idx-3] <- delta_H
    }
    if(delta_H<=1e-2) break
    mem_est <- max.col(gam_matrix,'first')
    diff[idx] <- sum(abs(mem_true-mem_est))
    ptFinal<-proc.time()-pt
  }
  print(diff)
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  
  resid <- W_test%*%coeff-(W_test%*%(Theta.mu))
  mse <- norm(resid,type='2')/sqrt(N)
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
}
# Structured VI with Gibbssampling
SVIG <- function(y,X,mem_true,coeff,X_test,S){
  pt<-proc.time();
  burn <- round(S/10,1)
  Q <- dim(X)[2]; N <- dim(X)[1]
  W <- apply(X,2,standardize)
  W_test <- apply(X_test,2,standardize)
  WTW <- t(W)%*%W
  y <- y-mean(y)
  yTy <- sum(y*y)
  yT1 <- sum(y)
  WTy <- t(W)%*%y
  
  #========== Iinitial VI parameter values ==========#
  #===== error =====#
  inv_sigma2.lt <- c()
  inv_sigma2 <- 1/c(var(y))
  is2.b.lt <- c()
  is2.a <- N/2
  
  #===== coefficients =====#
  Theta.mu <- matrix(0,nrow=Q,ncol=1)
  Theta.mu.lt <- list()
  Theta_ThetaT <- c()
  B.inv.lt <- list()
  
  #===== spike and slab prior =====#
  gam_matrix <- matrix(1/2,Q,2)
  gam <- matrix(gam_matrix[,2],ncol=1)
  cov_mem <- as.numeric(sample.int(2,size=Q,replace=T,prob=c((Q-1)/Q,1/Q))==2)
  sample_mat <- matrix(NA,nrow=S,ncol=Q)
  
  #===== rho =====#
  a0 <- 1
  b0 <- Q
  E_logit <- logit((a0/b0))
  
  #===== Dispersion parameters =====#
  lamb <- 1/1; lamb2_2 <- lamb^2/2
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1/1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #========== Optimization ==========#
  num_iters <- 20
  pred_hist <- c()
  delta_H <- 1000
  diff <- c()
  Delta_H <- c()
  
  for(idx in 1:num_iters){
    #===== gam update =====#
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
    
    # stopping criterion
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
  mem_true.x <- mem_true[1:Q]
  mem2.x <- mem_est[which(mem_true.x==2)]
  mem1.x <- mem_est[which(mem_true.x==1)]
  len2.x <- sum(mem_true.x==2)
  miss_rate.x <- c(sum(mem2.x != rep(2,len2.x)),
                   sum(mem1.x != rep(1,Q-len2.x)),
                   sum(mem_true.x != mem_est))
  
  resid <- W_test%*%coeff-(W_test%*%(Theta.mu))
  mse <- norm(resid,type='2')/sqrt(N)
  
  # F1 score
  TP <- sum((mem_est+mem_true)==4)
  FP <- sum((mem_est-mem_true)==1)
  FN <- sum((mem_est-mem_true)==-1)
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  F1 <- (2*precision*recall)/(precision+recall)
  
  return(list(c(FP,FN,miss_rate.x[3]),mem_est,mse,F1,ptFinal[3]))
}

################################################
################################################

################################################
################ Do Simulations ################
MISS1 <- matrix(0,nrow=0,ncol=3)
MISS2 <- matrix(0,nrow=0,ncol=3)
MISS3 <- matrix(0,nrow=0,ncol=3)

time1 <- time2 <- time3 <- c()
mse1 <- mse2 <- mse3 <- c()
F1 <- F2 <- F3 <- c()

corr<-0.8; P<-12; N<-30; valid_p<-3; sigma<-1

SIM <- 100

for(r in 1:SIM){
  print('SIMULATION : '); print(r)
  ################################################################
  ######################## generate sample #######################
  # List <- sim_data_cor_all(P,corr,N,sigma,r,valid_p)
  # List <- sim_data_cor_ar1(P,corr,N,sigma,r,valid_p)
  List <- sim_data_cor_block(P,corr,N,sigma,r,valid_p)
  y <- List[[1]]; X <- List[[2]]
  data <- List[[3]]; mem_true <- List[[4]]; coeff <- List[[5]]; X_test <- List[[6]]
  ################################################################
  ################################################################
  
  ################################################################
  ########################## Simulations #########################
  print('MFVI'); pt<-proc.time(); ls1 <- MFVI(y,X,mem_true,coeff,X_test,1); ptFinal<-proc.time()-pt
  time1 <- c(time1,c(ptFinal[3])); print(c(ptFinal[3]))
  
  print('SVIE'); pt<-proc.time(); ls2 <- SVIE(y,X,mem_true,coeff,X_test); ptFinal<-proc.time()-pt
  time2 <- c(time2,c(ptFinal[3])); print(c(ptFinal[3]))

  S <- 30
  print('SVIG'); pt<-proc.time(); ls3 <- SVIG(y,X,mem_true,coeff,X_test,S); ptFinal<-proc.time()-pt
  time3 <- c(time3,c(ptFinal[3])); print(c(ptFinal[3]))

  ################################################################
  ################################################################
  
  ################################################################
  ########################## store results #######################
  mse1 <- c(mse1,ls1[[3]])
  mse2 <- c(mse2,ls2[[3]])
  mse3 <- c(mse3,ls3[[3]])
  
  F1 <- c(F1,ls1[[4]])
  F2 <- c(F2,ls2[[4]])
  F3 <- c(F3,ls3[[4]])
  
  MISS1 <- rbind(MISS1,ls1[[1]])
  MISS2 <- rbind(MISS2,ls2[[1]])
  MISS3 <- rbind(MISS3,ls3[[1]])
  ################################################################
  ################################################################
  
  ################################################################
  ####################### Printout results########################
  print('Current FP, FN, and total missclassification : ')
  print(round(tail(MISS1,1),4));
  print(round(tail(MISS2,1),4));
  print(round(tail(MISS3,1),4));
  
  print('Average FP, FN, and total missclassification : ')
  print(round(colMeans(MISS1),4));
  print(round(colMeans(MISS2),4));
  print(round(colMeans(MISS3),4));
  ################################################################
  ################################################################
}


paste0('MVF : ',round(mean(F1,na.rm=T),2),' (',round(sd(F1,na.rm=T),2),')',' ',round(mean(mse1),2),' (',round(sd(mse1),2),')')
paste0('SVF : ',round(mean(F2,na.rm=T),2),' (',round(sd(F2,na.rm=T),2),')',' ',round(mean(mse2),2),' (',round(sd(mse2),2),')')
paste0('SVFmc : ',round(mean(F3,na.rm=T),2),' (',round(sd(F3,na.rm=T),2),')',' ',round(mean(mse3),2),' (',round(sd(mse3),2),')')


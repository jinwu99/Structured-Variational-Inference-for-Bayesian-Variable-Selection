rm(list=ls())
library(Matrix)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(mnormt)
library(sparsevb)
library(SSLASSO)
library(varbvs)
library(rmutil)
library(coda)
library(EMVS)
sourceCpp('MEAN_FIELD_ALGORITHM.cpp')
sourceCpp('SVI_GIBBS.cpp')
# sourceCpp('LOG_Q.cpp')

################################################
############## global functions ################
# normalize <- function(x) { return ((x - min(x)) / (max(x) - min(x))) }
center <- function(x) { return(x-mean(x))}
standardize <- function(x) { return ((x - mean(x)) / sd(x))}
tol <- 1e-200
normalizing <- function(mat) mat/rowSums(mat)
mat_trace <- function(A) sum(diag(A))
ar1_cor <- function(n, rho) rho^abs(matrix(1:n-1,nrow=n,ncol=n,byrow=TRUE)-(1:n-1))
logit <- function(m) log((m+tol)/(1+m+tol))
log_sum_exp <- function(log_vec) log(sum(exp(log_vec-max(log_vec)))+tol) + max(log_vec)
create_interaction_mat <- function(X,q){ # for adding interaction & squared terms
  p <- dim(X)[2]
  W <- apply(X,2,center)
  for(j1 in 1:(p-1)){
    for(j2 in (j1+1):p){
      W <- cbind(W,W[,j1]*W[,j2])
    }
  }
  for(j in 1:q){
    W <- cbind(W,W[,j]*W[,j])
  }
  print(dim(W)[2])
  return(W)
}
################################################
################################################

################################################
################## Algorithms ##################
# Mean Field VI
MFVI <- function(y,X,xi,train_idx,test_idx,lamb){
  pt<-proc.time()
  W <- apply(X,2,standardize)
  y <- y-mean(y)
  
  y_t <- y[test_idx]
  y <- y[train_idx]
  W_t <- W[test_idx,]; 
  W <- W[train_idx,]
  
  Q <- dim(W)[2]; N <- dim(W)[1]
  WTW <- t(W)%*%W
  WTy <- t(W)%*%y
  
  #========== Iinitial VI parameter values ==========#
  # ===== error =====#
  inv_sigma2 <- 1/c(var(y))
  # inv_sigma2 <- 1000
  
  #===== coefficients =====#
  Theta.mu <- matrix(0,nrow=Q,ncol=1)
  
  #===== rho =====#
  A0 <- 1
  B0 <- Q
  
  #===== spike and slab prior =====#
  gam_matrix <- matrix(1/2,Q,2)
  gam_matrix <- matrix(0,Q,2); gam_matrix[,2] <- 1
  gam <- matrix(gam_matrix[,2],ncol=1)
  Gam <- diag(c(gam),ncol=Q)
  Gam_2 <- diag(c(gam^2),ncol=Q)
  Omega <- gam%*%t(gam) - Gam_2 + Gam
  WTWO <- WTW*Omega
  
  #===== Dispersion parameters =====#
  lamb2_2 <- lamb^2/2
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1/1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #===== Deterministic Annealing(if needed)=====#
  # xi <- 1
  xi <- 0.001
  
  #========== Optimization ==========#
  j1 <- matrix(1,ncol=1,nrow=N)
  num_iters <- 100
  pred_hist <- c()
  delta_H <- 1000
  
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
    is2.a <- (N+Q-1)/2
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
    
    # I did a mistake...! I didn't use ELBO as a stopping criterion...
    pred_hist[idx] <- norm(W%*%(gam*Theta.mu),'2')/sqrt(N)
    if(idx>1){
      delta_H <- abs(pred_hist[idx]-pred_hist[idx-1])/abs(pred_hist[idx-1])
    }
    if(delta_H<=1e-5) break
  }
  ptFinal<-proc.time()-pt
  resid <- y_t-(W_t%*%(gam*Theta.mu))
  rmse <- sqrt(sum(resid^2)/length(test_idx))
  
  return(list(mem_est,rmse,gam,Theta.mu,ptFinal[3]))
}
# Structured VI with Gibbssampling
SVIG <- function(y,X,train_idx,test_idx,lamb,S){
  pt<-proc.time();
  burn <- round(S/10,1)
  W <- apply(X,2,standardize)
  y <- y-mean(y)
  
  y_t <- y[test_idx]
  y <- y[train_idx]
  W_t <- W[test_idx,]; 
  W <- W[train_idx,]

  Q <- dim(W)[2]; N <- dim(W)[1] # Note that N is a number of training samples, not the whole data
  WTW <- t(W)%*%W
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
  lamb2_2 <- lamb^2/2
  p_tau <- 1/2; a_tau <- 2*lamb2_2
  inv_tau2 <- rep(1/1,Q)
  Inv_V <- diag(inv_tau2,Q)
  
  #========== Optimization ==========#
  num_iters <- 20
  pred_hist <- c()
  delta_H <- 1000
  Delta_H <- c()
  
  for(idx in 1:num_iters){
    #===== gam update =====#
    rand_cov <- sample(seq.int(from=1,to=Q),replace=F)-1
    LIST <- SVI_GIBBS(rand_cov, WTW+Inv_V, WTy, diag(Inv_V),
                      E_logit, is2.a, yTy,
                      cov_mem, S+burn)
    sample_mat <- LIST[[1]][(burn+1):(S+burn),]; cov_mem <- sample_mat[S,]
    Theta.mu.lt <- LIST[[2]][(burn+1):(S+burn),]
    B.inv.lt <- LIST[[3]][(burn+1):(S+burn)]
    inv_sigma2.lt <- LIST[[4]][(burn+1):(S+burn)]
    log_Q <- LIST[[5]][(burn+1):(S+burn)]
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
    
    #===== Dispersion paramter =====#
    b_tau <- diag(is2_Theta_ThetaT)
    inv_tau2 <- (sqrt(a_tau)*besselK(sqrt(a_tau*b_tau),p_tau+1))/
      (sqrt(b_tau)*besselK(sqrt(a_tau*b_tau),p_tau)) - ((2*p_tau)/b_tau)
    Inv_V <- diag(inv_tau2,Q)
    
    # stopping criterion
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
    
    pred_hist[idx] <-
      sum((1/4-lamb2_2)*tau2) + lamb2_2*sum(inv_tau2) + sum((1/2-is2_Theta.mu)*log_tau2) +
      -0.5*inv_sigma2*yTy + t(WTy)%*%is2_Theta.mu - 0.5*TMP +
      Q*log_1_rho + log_Z +
      log(2*lamb2_2)/2*sum(b_tau) + sum(log(besselK(sqrt(lamb2_2),b_tau)))

    # par(mfrow=c(2,1))
    # plot(c(gam))
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
      ptFinal<-proc.time()-pt
      break
    }
    mem_est <- max.col(gam_matrix,'first')
    ptFinal<-proc.time()-pt
  }
  resid <- y_t-(W_t%*%(Theta.mu))
  rmse <- sqrt(sum(resid^2)/length(test_idx))
  return(list(mem_est,rmse,gam,Theta.mu,ptFinal[3]))
}

SPVB <- function(y,X,train_idx,test_idx){
  W <- apply(X,2,standardize)
  y <- y-mean(y)
  
  y_t <- y[test_idx]
  y <- y[train_idx]
  W_t <- W[test_idx,]; 
  W <- W[train_idx,]
  
  Q <- dim(W)[2]; N <- dim(W)[1]
  
  pt<-proc.time()
  fit <- svb.fit(W,y, family ="linear",intercept=F,slab='laplace')
  ptFinal<-proc.time()-pt
  mem_est <- as.numeric(fit$gamma>0.5)+1
  gam <- fit$gamma
  Theta.mu <- fit$mu
  
  resid <- y_t-(W_t%*%(gam*Theta.mu))
  rmse <- sqrt(sum(resid^2)/length(test_idx))
  
  return(list(mem_est,rmse,gam,Theta.mu,ptFinal[3]))
}

varbvss <- function(y,X,train_idx,test_idx){
  W <- apply(X,2,standardize)
  y <- y-mean(y)
  
  y_t <- y[test_idx]
  y <- y[train_idx]
  W_t <- W[test_idx,]; 
  W <- W[train_idx,]
  
  Q <- dim(W)[2]; N <- dim(W)[1]
  pt<-proc.time()
  fit <- varbvs(X=W,Z=NULL,y=y,family='gaussian',tol=1e-4,maxiter=1e4,verbose=F)
  ptFinal<-proc.time()-pt
  mem_est <- as.integer(fit$pip>0.5)+1
  
  gam <- as.integer(mem_est==2)
  
  y_hat <- predict(fit, X=W_t, Z = NULL,family='gaussian')
  resid <- y_t-y_hat
  rmse <- sqrt(sum(resid^2)/length(test_idx))
  
  return(list(mem_est,rmse,gam,ptFinal[3]))
}

SSLASSOO <- function(y,X,train_idx,test_idx){
  W <- apply(X,2,standardize)
  y <- y-mean(y)
  
  y_t <- y[test_idx]
  y <- y[train_idx]
  W_t <- W[test_idx,]; 
  W <- W[train_idx,]
  
  Q <- dim(W)[2]; N <- dim(W)[1]
  
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
      resid <- y_t-(intercept + W_t%*%(Theta.mu))
      rmse <- sqrt(sum(resid^2)/length(test_idx))
      mse_lt[idx] <- rmse
    }
  }
  min_idx <- which.min(mse_lt)
  rmse <- mse_lt[min_idx]
  
  Theta.mu <- fit$beta[,min_idx]; intercept <- fit$intercept[,min_idx]
  mem_est <- rep(1,Q); mem_est[linear_idx] <- 2
  gam <- as.integer(mem_est==2)
  
  return(list(mem_est,rmse,gam,Theta.mu,ptFinal[3]))
}

EMVSS <- function(y,X,train_idx,test_idx){
  W <- apply(X,2,standardize)
  y <- y-mean(y)
  
  y_t <- y[test_idx]
  y <- y[train_idx]
  W_t <- W[test_idx,]; 
  W <- W[train_idx,]
  
  Q <- dim(W)[2]; N <- dim(W)[1]
  
  v0 = seq(0.01, 0.1, length.out = 10)
  v1 = 1
  a = 1
  b = Q
  
  pt<-proc.time()
  invisible(capture.output(fit <- EMVS(y,W,v0=v0,v1=v1,type=c('betabinomial'),independent=F,
                                       a=a,b=b)))
  ptFinal<-proc.time()-pt
  idxes <- which(c(EMVSbest(fit)$log_g_function==fit$log_g_function)==T)
  
  mse_lt <- c()
  for(idx in 1:length(idxes)){
    i <- idxes[idx]
    Theta.mu <- fit$betas[i,]
    resid <- y_t-(W_t%*%(Theta.mu))
    rmse <- sqrt(sum(resid^2)/length(test_idx))
    mse_lt[idx] <- rmse
  }
  min_idx <- which.min(mse_lt)
  rmse <- mse_lt[min_idx]
  gam <- fit$prob_inclusion[min_idx,]
  mem_est <- as.numeric(gam>0.5)+1
  
  return(list(mem_est,rmse,gam,Theta.mu,ptFinal[3]))
}

MCMC <- function(y,X,train_idx,test_idx,total){
  pt<-proc.time();
  W <- apply(X,2,standardize)
  y <- y-mean(y)
  
  y_t <- y[test_idx]
  y <- y[train_idx]
  W_t <- W[test_idx,]; 
  W <- W[train_idx,]
  
  Q <- dim(W)[2]; N <- dim(W)[1]
  WTW <- t(W)%*%W
  yTy <- sum(y*y)
  yT1 <- sum(y)
  WTy <- t(W)%*%y
  
  #========== Iinitial VI parameter values ==========#
  #===== error =====#
  inv_Sigma2 <- 1/c(var(y))
  inv_sigma2.lt <- c()
  is2.b.lt <- c()
  is2.a <- N/2
  
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
  under_p_val <- Q
  
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
    
    # if(idx%%5==0){ # for checking convergence in eye
    #   par(mfrow=c(2,1))
    #   plot(colMeans(sample_mat))
    #   plot(colMeans(Theta.lt))
    # }
    
    if(idx%%50==0){ # for checking convergence with geweke diag
      p_val <- pnorm(abs(geweke.diag(Theta.lt)$z),lower.tail=FALSE)*2
      under_p_val <- sum(p_val < 0.001)
      cat(idx,under_p_val,'\n')
    }
    if(under_p_val==0) break
  }
  ptFinal<-proc.time()-pt
  print(ptFinal[3])
  Theta_hat <- sample_mat*Theta.lt
  L <- dim(Theta_hat)[1]; L_eff <- round(L*0.5,0)
  
  resid <- y_t-(W_t%*%colMeans(tail(Theta_hat,L_eff)))
  rmse <- sqrt(sum(resid^2)/length(test_idx))
  return(list(mem_est,rmse,gam,Theta.mu,ptFinal[3]))
}
################################################
################################################

################################################
################ Read Data #####################
##### Mcdonald data
library('bestglm')
data(mcdonald)
y <- mcdonald$MORT
X <- mcdonald[,-ncol(mcdonald)]
X <- create_interaction_mat(X,dim(X)[2])

##### TopGear data
library('robustHD')
data("TopGear")
TopGear <- TopGear[complete.cases(TopGear),] 
y <- log(TopGear$MPG)
X <- TopGear[,-13]
using_vars <- colnames(X)[c(5,7,9,10,11,12,13,14,15,16,4,18,19,27,8,17,20,21,22,23,24,25,26,28,29,31)]
X <- X[,using_vars]
dim(X)[2]
library('fastDummies')
X <- dummy_cols(X,remove_first_dummy=T,remove_selected_columns=T)
# X[,c(1,2,3,4,6,7,8,9,10)] <- log(X[,c(1,2,3,4,6,7,8,9,10)])
X <- create_interaction_mat(X,10)

##### BostonHousing data
library('mlbench')
data(BostonHousing)
y <- BostonHousing$medv
X <- BostonHousing[,-dim(BostonHousing)[2]]
X$chas <- as.numeric(X$chas)-1

##### Ozone data
library('spikeslab')
data(ozoneI, package = "spikeslab")
y <- ozoneI$ozone
X <- as.matrix(ozoneI[,-1])

##### Eyedata
library('flare')
data(eyedata)
X <- x

##### Riboflavin data
library('hdi')
data(riboflavin)
dim(riboflavin)
y <- riboflavin$y
X <- as.matrix(riboflavin$x)
# Screening using simple linear regression 
n <- dim(X)[1]; p <- dim(X)[2]
storage <- matrix(NA,nrow=dim(X)[2],ncol=2)
colnames(storage) <- c('coeffs','p_values')
for(j in 1:dim(X)[2]){
  print(j)
  lm0 <- lm(y~X[,j])
  coeff <- summary(lm0)$coefficients[,1][2]
  p_value <- summary(lm0)$coefficients[,4][2]
  storage[j,1] <- coeff; storage[j,2] <- p_value
}
sort_idx <- sort(storage[,2],index.return=T,decreasing=F)$ix
X <- X[,sort_idx[1:500]]
################################################
################################################

################################################
################# Simulation ###################
n <- dim(X)[1]; p <- dim(X)[2]
fold <- 5
por <- 1/fold
n_test <- round(n*por,0); n_train <- n-n_test

all_idx <- 1:n
cov_list1 <- cov_list2 <- cov_list3 <- cov_list4 <- cov_list5 <- cov_list6 <- cov_list7 <- c()
mse_list1 <- mse_list2 <- mse_list3 <- mse_list4 <- mse_list5 <- mse_list6 <- mse_list7 <- c()
time1 <- time2 <- time3 <- time4 <- time5 <- time6 <- time7 <- c()

lamb <- 1; total <- 1000
for(i in 1:10){
  print(i)
  set.seed(i)

  test_idx <- sample.int(n,n_test,replace=F)
  train_idx <- all_idx[-test_idx]

  ls1 <- MFVI(y,X,1,train_idx,test_idx,lamb)
  mem_est1 <- ls1[[1]]; mse1 <- ls1[[2]]; gam1 <- ls1[[3]]; Theta.mu1 <- ls1[[4]]
  time1 <- c(ls1[[5]],time1)

  S <- 30
  ls2 <- SVIG(y,X,train_idx,test_idx,lamb,S)
  mem_est2 <- ls2[[1]]; mse2 <- ls2[[2]]; gam2 <- ls2[[3]]; Theta.mu2 <- ls2[[4]]
  time2 <- c(ls2[[5]],time2)

  ls3 <- MCMC(y,X,train_idx,test_idx,total)
  mem_est3 <- ls3[[1]]; mse3 <- ls3[[2]]; gam3 <- ls3[[3]]; Theta.mu3 <- ls3[[4]]
  time3 <- c(ls3[[5]],time3)
  
  ls4 <- SPVB(y,X,train_idx,test_idx)
  mem_est4 <- ls4[[1]]; mse4 <- ls4[[2]]; gam4 <- ls4[[3]]; Theta.mu4 <- ls4[[4]]
  time4 <- c(ls4[[5]],time4)

  ls5 <- varbvss(y,X,train_idx,test_idx)
  mem_est5 <- ls5[[1]]; mse5 <- ls5[[2]]; gam5 <- ls5[[3]]
  time5 <- c(ls5[[4]],time5)

  ls6 <- EMVSS(y,X,train_idx,test_idx)
  mem_est6 <- ls6[[1]]; mse6 <- ls6[[2]]; gam6 <- ls6[[3]]; Theta.mu6 <- ls6[[4]]
  time6 <- c(ls6[[5]],time6)

  ls7 <- SSLASSOO(y,X,train_idx,test_idx)
  mem_est7 <- ls7[[1]]; mse7 <- ls7[[2]]; gam7 <- ls7[[3]]; Theta.mu7 <- ls7[[4]]
  time7 <- c(ls7[[5]],time7)
  
  mse_list1[i] <- mse1
  cov_list1[i] <- sum(gam1>0.5)
  mse_list2[i] <- mse2
  cov_list2[i] <- sum(gam2>0.5)
  mse_list3[i] <- mse3
  cov_list3[i] <- sum(gam3>0.5)
  mse_list4[i] <- mse4
  cov_list4[i] <- sum(gam4>0.5)
  mse_list5[i] <- mse5
  cov_list5[i] <- sum(gam5>0.5)
  mse_list6[i] <- mse6
  cov_list6[i] <- sum(gam6>0.5)
  mse_list7[i] <- mse7
  cov_list7[i] <- sum(gam7>0.5)

  par(mfrow=c(3,1))
  boxplot(mse_list1,mse_list2,mse_list3,mse_list4,mse_list5,mse_list6,mse_list7)
  boxplot(cov_list1,cov_list2,cov_list3,cov_list4,cov_list5,cov_list6,cov_list7)
  boxplot(time1,time2,time3,time4,time5,time6,time7)
  cat(mean(mse_list1,na.rm=T),mean(mse_list2,na.rm=T),mean(mse_list3,na.rm=T),
      mean(mse_list4,na.rm=T),mean(mse_list5,na.rm=T),mean(mse_list6,na.rm=T),mean(mse_list7,na.rm=T),'\n')
  cat(mean(cov_list1,na.rm=T),mean(cov_list2,na.rm=T),
      mean(cov_list4,na.rm=T),mean(cov_list5,na.rm=T),mean(cov_list6,na.rm=T),mean(cov_list7,na.rm=T),'\n')
}

cat(mean(mse_list1),mean(mse_list2),mean(mse_list3),mean(mse_list4),mean(mse_list5),mean(mse_list6),mean(mse_list7),'\n')
cat(sd(mse_list1),sd(mse_list2),sd(mse_list3),sd(mse_list4),sd(mse_list5),sd(mse_list6),sd(mse_list7),'\n')
cat(mean(cov_list1),mean(cov_list2),mean(cov_list3),mean(cov_list4),mean(cov_list5),mean(cov_list6),mean(cov_list7),'\n')
cat(sd(cov_list1),sd(cov_list2),sd(cov_list3),sd(cov_list4),sd(cov_list5),sd(cov_list6),sd(cov_list7),'\n')
cat(mean(time1),mean(time2),mean(time3),mean(time4),mean(time5),mean(time6),mean(time7),'\n')
cat(sd(time1),sd(time2),sd(time3),sd(time4),sd(time5),sd(time6),sd(time7),'\n')



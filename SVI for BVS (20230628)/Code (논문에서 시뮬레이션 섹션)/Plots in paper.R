rm(list=ls())
library(FSA)
library(ggplot2)
library(dplyr)
library(ggpubr)

#####################################################################
# F1 boxplot
N_train <- 70; N <- round(N_train*(5/4),0)
Q<-200; p<-10; sigma<-sqrt(1)
SIM<-10
Corr_struct_list <- c('ar1','block','all')
Corr_struct_title_list <- c('AR(1)','Block','Pairwise')
Method.lt <- c('SVI+G','MFVI','EMVS','VARBVS','SSLASSO','SPVB','MCMC')
DATA_total <- list()
for(s in 1:3){
  Corr_struct <- Corr_struct_list[[s]]
  DATA.lt <- list()
  for(c in 1:2){
    corr <- 0.4*c
    storage_text <- paste0('P',Q,'_','N',N,'_','p',p,'_',Corr_struct,'_',corr)
    storage_place <- paste0("C:\\Users\\Public\\",storage_text,'.RData')
    load(file=storage_place)
    
    M <- length(save_file$F1_list)
    data.lt <- list()
    for(m in 1:M){
      SIM <- length(save_file$F1_list[[m]])
      data.lt[[m]] <- data.frame(F1=save_file$F1_list[[m]],
                                 Method=rep(Method.lt[[m]],SIM),
                                 Corr=rep(paste0(Corr_struct_title_list[s],' : ',corr),SIM))
    }
    DATA.lt[[c]] <- Reduce(rbind,data.lt)
  }
  DATA <- Reduce(rbind,DATA.lt)
  DATA$Struct <- Corr_struct
  DATA_total[[s]] <- DATA
}

pd = position_dodge(.5)
Sum.lt <- list()
gp.lt <- list()
# ymin <- c(0.4,0.05,0.05)
for(s in 1:3){
  Sum.lt[[s]] <- DATA_total[[s]] %>% 
    as.data.frame() %>%
    mutate(Method=factor(Method,levels=c('MFVI','EMVS','VARBVS','SSLASSO','SPVB','SVI+G','MCMC')))
  
  gp.lt[[s]] <- ggplot(Sum.lt[[s]],
                       aes(x = Method,
                           y = F1,
                           fill = Method)) +
    facet_wrap(vars(Corr), ncol=2) +
    geom_boxplot(outlier.shape = NA,
                 position = pd) +
    labs(y="F1-score",x='') +
    # ylim(ymin[s],0.98) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}
g1 <- ggarrange(gp.lt[[1]], gp.lt[[2]], gp.lt[[3]],
                ncol=3, nrow=1, common.legend = TRUE, legend="right")

#####################################################################
#####################################################################
# RMSE boxplot
DATA_total <- list()
for(s in 1:3){
  Corr_struct <- Corr_struct_list[[s]]
  DATA.lt <- list()
  for(c in 1:2){
    corr <- 0.4*c
    storage_text <- paste0('P',Q,'_','N',N,'_','p',p,'_',Corr_struct,'_',corr)
    storage_place <- paste0("C:\\Users\\Public\\",storage_text,'.RData')
    load(file=storage_place)
    
    M <- length(save_file$F1_list)
    data.lt <- list()
    for(m in 1:M){
      SIM <- length(save_file$F1_list[[m]])
      data.lt[[m]] <- data.frame(RMSE=save_file$RMSE_list[[m]],
                                 Method=rep(Method.lt[[m]],SIM),
                                 Corr=rep(paste0(Corr_struct_title_list[s],' : ',corr),SIM))
    }
    DATA.lt[[c]] <- Reduce(rbind,data.lt)
  }
  DATA <- Reduce(rbind,DATA.lt)
  DATA$Struct <- Corr_struct
  DATA_total[[s]] <- DATA
}

pd = position_dodge(.5)
Sum.lt <- list()
gp2.lt <- list()
ymax <- c(2.5,2,2.2)
for(s in 1:3){
  Sum.lt[[s]] <- DATA_total[[s]] %>% 
    as.data.frame() %>%
    mutate(Method=factor(Method,levels=c('MFVI','EMVS','VARBVS','SSLASSO','SPVB','SVI+G','MCMC')))
  
  gp2.lt[[s]] <- ggplot(Sum.lt[[s]],
                        aes(x = Method,
                            y = RMSE,
                            fill = Method)) +
    facet_wrap(vars(Corr), ncol=2) +
    geom_boxplot(outlier.shape = NA,
                 position = pd) +
    labs(y="RMSE",x='') +
    ylim(0.2,ymax[s]) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}
g2 <- ggarrange(gp2.lt[[1]], gp2.lt[[2]], gp2.lt[[3]],
                ncol=3, nrow=1, common.legend = TRUE, legend="right")

#####################################################################
#####################################################################
# Time boxplot
DATA_total <- list()
for(s in 1:3){
  Corr_struct <- Corr_struct_list[[s]]
  DATA.lt <- list()
  for(c in 1:2){
    corr <- 0.4*c
    storage_text <- paste0('P',Q,'_','N',N,'_','p',p,'_',Corr_struct,'_',corr)
    storage_place <- paste0("C:\\Users\\Public\\",storage_text,'.RData')
    load(file=storage_place)
    
    M <- length(save_file$TIME_list)
    data.lt <- list()
    for(m in 1:M){
      SIM <- length(save_file$TIME_list[[m]])
      data.lt[[m]] <- data.frame(TIME=log(save_file$TIME_list[[m]]),
                                 Method=rep(Method.lt[[m]],SIM),
                                 Corr=rep(paste0(Corr_struct_title_list[s],' : ',corr),SIM))
    }
    DATA.lt[[c]] <- Reduce(rbind,data.lt)
  }
  DATA <- Reduce(rbind,DATA.lt)
  DATA$Struct <- Corr_struct
  DATA_total[[s]] <- DATA
}

Sum.lt <- list()
gp3.lt <- list()
for(s in 1:3){
  Sum.lt[[s]] <- DATA_total[[s]] %>% 
    as.data.frame() %>%
    mutate(Method=factor(Method,levels=c('MFVI','EMVS','VARBVS','SSLASSO','SPVB','SVI+G','MCMC')))
  
  gp3.lt[[s]] <- ggplot(Sum.lt[[s]],
                        aes(x = Method,
                            y = TIME,
                            fill = Method)) +
    facet_wrap(vars(Corr), ncol=2) +
    geom_boxplot(outlier.shape = NA,
                 position = pd) +
    labs(y="log(TIME)",x='') +
    # ylim(ymin[s],0.98) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}
g3 <- ggarrange(gp3.lt[[1]], gp3.lt[[2]], gp3.lt[[3]],
                ncol=3, nrow=1, common.legend = TRUE, legend="right")

# Merge plots
g4 <- ggarrange(g1,g2,g3,nrow=3,common.legend=T,legend='right')





#Load packages needed
library(Rcpp)
library(RcppParallel)
library(matrixStats)
library(tidyverse)

########
## function to compute the mcmc output
para_summary <- function(mcmc,a,b,print){
  y <- matrix(NA,ncol(mcmc),4)
  for (i in 1:ncol(mcmc)){
    y[i,1:3] <- quantile(mcmc[,i],c(0.5,0.025,0.975))
    y[i,4] <- sum(diff(mcmc[,i])!=0)/nrow(mcmc)
  }
  layout(matrix(1:(a*b),nrow=a,byrow=T))
  par(mar=c(2,4,1,1))
  if (print==1){
    for (i in 1:ncol(mcmc)){
      plot(mcmc[,i],type="l")
    }
  }
  return(y)
}

#function to list possible infection time scenarios
#considering probability threshold
hh_sep <- function(hh, incu, thres){
  inf <- c()
  inf_time <- c()
  end_time <- hh[5]
  norm_prob <- 1
  for(m in c(0:(as.integer(length(hh)/10)-1))){
    if(hh[10*m+7]==1){
      inf <- c(inf, m)
      inf_time <- c(inf_time, hh[10*m+8])
    }
  }
  if(max(incu)^length(inf)<thres){
    print(paste0("Drop the household ", hh[1]))
    return(rep(NA,length(hh)))
  }
  out <- hh
  vec <- rep(1, length(inf))
  while(vec[length(inf)] <=length(incu)){
    prob <- 1
    for(i in c(1:length(inf))){
      prob <- prob*incu[vec[i]]
    }
    if(max(inf_time - vec) > end_time){
      norm_prob <- norm_prob - prob
    }
    if(prob >= thres & max(inf_time - vec) <= end_time){
      temp <- hh
      start <- temp[5]
      for(member in c(1:length(inf))){
        temp[10*inf[member]+9] <- temp[10*inf[member]+8] - vec[member]
        temp[10*inf[member]+16] <- 0
        if(start > temp[10*inf[member]+9]){
          start <- temp[10*inf[member]+9]
        }
      }
      temp_flag <- 0
      for(member in c(1:length(inf))){
        if(start == temp[10*inf[member]+9]){
          if(temp_flag==1){
            temp[10*inf[member]+16] <- 2
          }
          if(temp_flag==0){
            temp[10*inf[member]+16] <- 1
            temp_flag <- 1
          }
        }
      }
      temp[4] <- start
      temp[6] <- log(prob)
      out <- rbind(out, temp)
    }
    vec[1] <- vec[1]+1
    if(length(inf)>1){
      for(j in c(1:(length(inf)-1))){
        if(vec[j]==length(incu)+1){
          vec[j] <- 1
          vec[j+1] <- vec[j+1]+1
        }
      }
    }
    
  }
  if(!is.null(nrow(out))){
    out[2:nrow(out),6] <- out[2:nrow(out),6] -log(norm_prob)
    out <- out[2:nrow(out),]
    return(out)
  }
  if(is.null(nrow(out))){
    print(paste0("Drop the household ", hh[1]))
    return(rep(NA,length(hh)))
  }
}

#function to calculate DIC value
DIC_comp <- function(data1,para_est,inci,incu,shift,prob){
  inc <- 10000+1:10000*2
  data_sep <- hh_sep(data1[1,],incu,prob)
  for(i in c(2:nrow(data1))){
    sep <- hh_sep(data1[i,],incu,prob)
    data_sep <- rbind(data_sep, sep)
  }
  data_sep <- data_sep[!is.na(data_sep[,1]),]
  print("Data Seperation: DONE")
  
  mean_para <- colMeans(para_est)
  mean_loglik <- sum(loglik(data_sep,inci,mean_para,shift,6,10)[[1]])
  #print(paste0("Loglik Mean: ", mean_loglik))
  loglik_list <- rep(0,10000)
  for(i in c(1:10000)){
    loglik_list[i] <- sum(loglik(data_sep,inci,para_est[i,],shift,6,10)[[1]])
  }
  #print(paste0("Mean Loglik: ", mean(loglik_list)))
  dic <- -2*(2*mean(loglik_list) - mean_loglik)
  return(dic)
}

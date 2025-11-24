# Data for the paper Estimating pre-symptomatic transmission potential of influenza A and B viruses in household transmission studies
# Code for influenza B model

rm(list = ls())
source("helper_funcs.R")
sourceCpp("infB_model.cpp")

#incubation period
incu <- incubation(log(0.6),log(1.51))

inci <- rep(1,4000)
data <- read.csv("infB_hk_data.csv",as.is = T)
data1 <- as.matrix(data)

#List all possible infection time scenarios
data_sep <- hh_sep(data1[1,],incu,0.0001)
for(i in c(2:nrow(data1))){
  sep <- hh_sep(data1[i,],incu,0.0001)
  data_sep <- rbind(data_sep, sep)
}
data_sep <- data_sep[!is.na(data_sep[,1]),]
print("Listing scenarios: DONE")

#Initial values of model parameters
para <- c(log(0.6), log(1.51), 3, 1, 1, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
          0, 0, 0, 
          0, 0, 0, 
          0, 0, 0)

# input the initial s.d. for the metropolis hasting
sigma <- c(1,rep(0.1,30))
# specific if the parameter is sampled
move <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 
          1, 1, 0, 
          0, 0, 0, 
          0, 0, 0)

#Model inference
#shift = 1, infectiousness starting at symptom onset
aaaaa1 <- Sys.time()
tt <- mcmc(data_sep,inci,1,30000,para,move,sigma)
aaaaa2 <- Sys.time()
print(aaaaa2-aaaaa1)

inc <- 10000+1:10000*2
dic <- DIC_comp(data1,tt[[1]][inc,],inci,incu,1,0.0001)
tt[[3]] <- dic

tt_infB <- tt
save(tt_infB, file = "tt_infB.RData")

z1 <- para_summary((tt[[1]][inc,]),4,4,1)
write.csv(z1,paste("mcmc_summary_infB.csv",sep=""))


##########################################################################################
## for compute model adequency
## this will generate a file that can be coverted to the model adequency check in appendix
## create a function to generate a model adequency check
count <- function(temp){
  out <- matrix(NA,1,14*14)
  infn <- rowSums(temp[,0:7*10+7]==1)
  
  for (i in 0:13){
    for (j in 0:13){
      out[1,14*i+j+1] <- sum(temp[,2]==i & infn==j)  
    }  
  }
  
  out
}

# the first row is the data
adeq <- matrix(NA,10000,196)
real <- count(data1)


a1 <- tt[[1]][inc,]
#a2 <- tt[[3]][1:10000,]

for (i in 1:10000){
  if (i%%100==0){print(i)}
  uu <- sim_data(data1,inci,a1[i,],1,6,10)
  adeq[i,] <- count(uu[[1]])
}

z3 <- para_summary(adeq,4,3,0)
z3[,4] <- real

write.csv(z3,paste("adeq_infB.csv",sep=""))
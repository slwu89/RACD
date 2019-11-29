rm(list=ls());gc()
library(here)

source(here("/mosy-equilibrium.R"))
Rcpp::sourceCpp(here('modular/mosquito.cpp'))

parameters <- c(
  ## Mosquito life cycle parameters:
  beta = 21.19, # Number of eggs laid per day by female mosquito
  muEL = 0.034, # Early larval instar daily mortality
  muLL = 0.035, # Late larval instar daily mortality
  muPL = 0.25, # Pupal daily mortality
  durEL = 6.64, # Duration of early instar stage
  durLL = 3.72, # Duration of late instar stage
  durPL = 0.64, # Duration of pupal stage
  durEV = 10, # Duration of latent period in mosquito (days)
  gamma = 13.25, # Effect of density-dependence on late instars relative to early instars
  tau1 = 0.68, # Time spent foraging for a blood meal at 0% ITN coverage
  tau2 = 2.32, # Time spent resting and ovipositing by a mosquito

  ## Species-specific parameters:
  ## An. gambiae:
  muV = 1/7.6, # Adult mosquito daily mortality
  Q0 = 0.92, # Human blood index
  phiB = 0.89, # Proportion of bites on a person while they are in bed
  phiI = 0.97, # Proportion of bites on a person while they are indoors
  rITN = 0.56, # Probability of mosquito repeating a feeding attempt due to IRS
  sITN = 0.03, # Probability of mosquito feeding and surviving in presence of ITNs
  rIRS = 0.60, # Probability of mosquito repeating a feeding attempt due to IRS
  sIRS = 0, # Probability of mosquito feeding and surviving in presence of IRS

  ## An. arabiensis:
  # muV = 1/7.6, # Adult mosquito daily mortality
  # Q0 = 0.71, # Human blood index
  # phiB = 0.90, # Proportion of bites on a person while they are in bed
  # phiI = 0.96, # Proportion of bites on a person while they are indoors
  # rITN = 0.48, # Probability of mosquito repeating a feeding attempt due to IRS
  # sITN = 0.39, # Probability of mosquito feeding and surviving in presence of ITNs
  # rIRS = 0.60, # Probability of mosquito repeating a feeding attempt due to IRS
  # sIRS = 0, # Probability of mosquito feeding and surviving in presence of IRS

  ## An. funestus:
  # muV = 1/8.9, # Adult mosquito daily mortality
  # Q0 = 0.94, # Human blood index
  # phiB = 0.90, # Proportion of bites on a person while they are in bed
  # phiI = 0.98, # Proportion of bites on a person while they are indoors
  # rITN = 0.56, # Probability of mosquito repeating a feeding attempt due to IRS
  # sITN = 0.03, # Probability of mosquito feeding and surviving in presence of ITNs
  # rIRS = 0.63, # Probability of mosquito repeating a feeding attempt due to IRS
  # sIRS = 0, # Probability of mosquito feeding and surviving in presence of IRS

  ## Additional transmission parameters:
  f0 = 1/3 # Daily biting rate by mosquitoes on animals and humans

)

delta <- 1.0/(parameters["tau1"]+parameters["tau2"]) # Inverse of gonotrophic cycle without ITNs/IRS
eggOV <- parameters["beta"]*(exp(parameters["muV"]/delta)-1)/parameters["muV"] # Number of eggs per oviposition per mosquito
parameters["eggOV"] <- eggOV

a <- parameters["f0"]*parameters["Q0"]

set.seed(2134124325L)
n <- 500
hh <- sample(x = 1:100,size = n,replace = T)
psi <- rexp(length(unique(hh))); psi <- psi/sum(psi)


pi <- sapply(unique(hh),function(h){
  pi <- rlnorm(n = sum(hh == h))
  pi/sum(pi)
})

w <- sapply(pi,function(p){
  runif(length(p))
})
c <- sapply(pi,function(p){
  rbeta(n = length(p),shape1 = 5,shape2 = 25)
})

# house infectivity
h_cc <- mapply(FUN = function(c,p,w){
  sum(c*(p*w)/sum(p*w))
},c=c,p=pi,w=w)

lambda <- a * sum(psi * h_cc)


eq <- RACD_mosq_equilibrium(theta = parameters,dt = 1,IV = 500,lambdaV = lambda)

mosy_ptr <- init_mosquitos(EL_ = eq$EL_eq,LL_ = eq$LL_eq,PL_ = eq$PL_eq,
                           SV_ = eq$SV_eq,EV_ = eq$EV_eq,IV_ = eq$IV_eq,
                           K_ = eq$K_eq,lambdaV_ = lambda)

bite_probs <- rep(0,3)
names(bite_probs) <- c("WW","ZZ","CC")
bite_probs[1] <- 1
bite_probs[2] <- 0.001
bite_probs[3] <- sum(h_cc*psi)

out <- matrix(data = 0,nrow = 1e3,ncol = 3)
for(i in 1:1e3){
  feeding_cycle(mosy = mosy_ptr,bite_probs = bite_probs,psi = psi,parameters = parameters)
  euler_step(mosy = mosy_ptr,parameters = parameters)
  out[i,] <- track_mosquito(mosy = mosy_ptr)
}

matplot(out,type="l",lty=1)
